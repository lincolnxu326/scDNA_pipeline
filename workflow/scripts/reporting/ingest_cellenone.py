#!/usr/bin/env python3
"""
Parse a CellenONE `.Run` folder into per-well records the review report can use.

Parse-only and cheap: this reads the run's tables and indexes the image filenames,
but never opens an image. All pixel work lives in render_cellenone_images.py.

What the run folder actually contains
-------------------------------------
* `Reordered_*_isolated.xls` — despite the extension this is **tab-separated text**
  with CRLF line endings and a run of empty trailing fields. One row per
  (drop, channel): 384 Transmission + 384 of each fluorescence channel the run used.
  Which channels exist varies per run: plate21/22 recorded Blue + Orange (+ a few
  Red), plate24 recorded Green + Orange + Red. Carries
  the isolated cell's X/Y in the 952x471 camera frame, its diameter/elongation/
  circularity/intensity, and — importantly — the run's own detection and isolation
  criteria, which **vary per run** and must always be read from the file rather than
  hardcoded (plate21 isolates 12-25 um, plate22 9-25 um).
* `Reordered_*_geoprops.xls` — every detected object in every drop, printed or not.
  For a printed drop it holds exactly ONE Transmission row (the isolated cell), so it
  cannot answer "how many objects are in this well's image" — that is what the
  detection pass in render_cellenone_images.py is for. Used here only for run-level
  statistics (drops attempted, isolation frequency).
* `*.par` — the instrument parameter file; `[Cells] TParameters` restates the
  detection/isolation thresholds as JSON.
* `<drop>_Printed_..._({POS})_Trans_...png` and `({POS})Blue` / `({POS})Orange` —
  the raw grayscale camera images, matched to wells by the `(POS)` token in the
  filename. NEVER by count or index: this run has 385 Blue and 385 Orange PNGs but
  only 384 table rows each, and Red exists for just 27 of 384 wells.

The two ejection lines
----------------------
Objects to the RIGHT of the purple line are still up the capillary and will not be
ejected; the nozzle tip is at the LEFT and cells sediment toward it. The lines are
resolved in this order, and `run_meta.json` records which source won:

  1. hue auto-detect on an annotated `cellenREPORT/Images_iso_det/Trans_*.jpg`
  2. the `BackgroundEjZone` image width (purple only — it IS the zone's right edge)
  3. the config defaults (396 / 647)

Note the `EjBound` column in the xls is 330 and is **not** the green line.

Outputs
-------
    <plate>/cellenone/run_meta.json   lines + their source, thresholds, run stats
    <plate>/cellenone/wells_raw.tsv   one row per sequencing well
"""

import argparse
import csv
import json
import logging
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plate384_layout as L384

# `(A-1)` / `(P-24)` in an image filename, followed by the channel token:
#   …_(K-22)_Trans_K1563_plate_2_Run.png
#   …_(K-22)Blue_K1563_plate_2_Run.png     <- note: no separator before Blue/Orange
# The trailing guard is `(?![A-Za-z])`, NOT `\b`: the next character is `_`, which is
# itself a word character, so `\b` never fires there.
IMG_RE = re.compile(r"\(([A-Pa-p]-\d{1,2})\)_?(Trans|Blue|Green|Orange|Red)(?![A-Za-z])", re.I)
HYPERLINK_RE = re.compile(r'=HYPERLINK\("([^"]*)"\)', re.I)
CHANNEL_OF_TEG = {"transmission": "trans", "blue": "blue", "green": "green",
                  "orange": "orange", "red": "red"}


def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    return logging.getLogger("cellenone_ingest")


def parse_args():
    p = argparse.ArgumentParser(description="Ingest a CellenONE run folder.")
    p.add_argument("--run-dir", required=True, help="The `.Run` folder for this plate")
    p.add_argument("--plate", required=True, help="Plate name")
    p.add_argument("--config", required=True, help="Pipeline config.yaml")
    p.add_argument("--outdir", required=True, help="Output dir (<plate>/cellenone)")
    p.add_argument("--subplates", nargs="*", default=[],
                   help="Subplate names in SL order; empty => 96-well plate")
    p.add_argument("--layout-tsv", default="", help="Optional layout override TSV")
    p.add_argument("--wells", nargs="+", required=True, help="Well IDs (W01..W96)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# The CellenONE "xls" files (tab-separated text)
# ---------------------------------------------------------------------------

def read_cellenone_table(path: Path) -> list:
    """Read a CellenONE *.xls (really TSV) into a list of dicts.

    Handles the three things that make these files awkward: CRLF endings, a long run
    of empty trailing columns, and `=HYPERLINK("name.png")`-wrapped file references.
    """
    rows = []
    with open(path, newline="", encoding="latin-1") as fh:
        reader = csv.reader(fh, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return rows
        # Trim the empty trailing fields, and strip the padding CellenONE writes into
        # some header names ("Date        ").
        names = [h.strip() for h in header]
        while names and not names[-1]:
            names.pop()
        n = len(names)
        for parts in reader:
            if not parts or all(not c.strip() for c in parts):
                continue
            row = {}
            for i, key in enumerate(names):
                val = parts[i].strip() if i < len(parts) else ""
                m = HYPERLINK_RE.match(val)
                if m:
                    val = m.group(1)
                row[key] = val
            rows.append(row)
    return rows


def as_float(v):
    try:
        return float(str(v).strip())
    except (TypeError, ValueError):
        return None


def as_int(v):
    f = as_float(v)
    return int(f) if f is not None else None


# ---------------------------------------------------------------------------
# Line resolution
# ---------------------------------------------------------------------------

def detect_lines_from_jpeg(run_dir: Path, logger):
    """Find the green and magenta guide columns in an annotated Trans_*.jpg.

    The annotated report JPEGs are grayscale frames with coloured vertical lines drawn
    on, so a saturation threshold isolates the overlay cleanly and the hue then says
    which line is which. Returns (green_x, purple_x) with either possibly None.
    """
    iso_dir = run_dir / "cellenREPORT" / "Images_iso_det"
    if not iso_dir.is_dir():
        return None, None
    try:
        import numpy as np
        from PIL import Image
    except ImportError:
        logger.warning("numpy/pillow unavailable; skipping line auto-detect")
        return None, None

    candidates = sorted(iso_dir.glob("Trans_*.jpg"))[:8]
    greens, purples = [], []
    for path in candidates:
        try:
            arr = np.asarray(Image.open(path).convert("RGB"), dtype=np.int16)
        except OSError:
            continue
        sat = arr.max(2) - arr.min(2)
        strong = sat > 60
        cols = np.where(strong.sum(0) > arr.shape[0] * 0.5)[0]
        g_cols, p_cols = [], []
        for c in cols:
            px = arr[strong[:, c], c, :]
            if px.size == 0:
                continue
            r, g, b = px.mean(0)
            if g > r and g > b:
                g_cols.append(c)
            elif r > g and b > g:
                p_cols.append(c)
        if g_cols:
            greens.append(int(min(g_cols)))
        if p_cols:
            purples.append(int(min(p_cols)))

    def mode(vals):
        if not vals:
            return None
        counts = defaultdict(int)
        for v in vals:
            counts[v] += 1
        return max(counts.items(), key=lambda kv: (kv[1], -kv[0]))[0]

    return mode(greens), mode(purples)


def ejzone_width(run_dir: Path):
    """The BackgroundEjZone image's width IS the purple line (the zone's right edge)."""
    for path in sorted(run_dir.glob("*BackgroundEjZone*.png")):
        try:
            from PIL import Image
            with Image.open(path) as im:
                return int(im.size[0])
        except (ImportError, OSError):
            return None
    return None


def resolve_lines(run_dir: Path, cfg: dict, logger) -> dict:
    lines_cfg = (cfg.get("lines") or {})
    g_cfg, p_cfg = lines_cfg.get("green_x", "auto"), lines_cfg.get("purple_x", "auto")
    g_def = int(lines_cfg.get("green_default", 396))
    p_def = int(lines_cfg.get("purple_default", 647))

    green = purple = None
    g_src = p_src = "config"
    if g_cfg != "auto":
        green, g_src = int(g_cfg), "config_explicit"
    if p_cfg != "auto":
        purple, p_src = int(p_cfg), "config_explicit"

    if green is None or purple is None:
        gj, pj = detect_lines_from_jpeg(run_dir, logger)
        if green is None and gj is not None:
            green, g_src = gj, "jpeg_hue"
        if purple is None and pj is not None:
            purple, p_src = pj, "jpeg_hue"

    if purple is None:
        w = ejzone_width(run_dir)
        if w:
            purple, p_src = w, "backgroundejzone_width"

    if green is None:
        green, g_src = g_def, "config_default"
    if purple is None:
        purple, p_src = p_def, "config_default"

    logger.info("Ejection lines: green_x=%d (%s), purple_x=%d (%s)",
                green, g_src, purple, p_src)
    return {"green_x": green, "purple_x": purple,
            "green_source": g_src, "purple_source": p_src}


# ---------------------------------------------------------------------------
# Images, indexed by the (POS) token
# ---------------------------------------------------------------------------

def index_images(run_dir: Path, logger) -> dict:
    """{(row, col): {channel: filename}} built from the `(POS)` token in each name."""
    by_pos = defaultdict(dict)
    n_files = 0
    dupes = 0
    for path in run_dir.glob("*.png"):
        m = IMG_RE.search(path.name)
        if not m:
            continue
        try:
            pos = L384.parse_pos(m.group(1))
        except ValueError:
            continue
        channel = m.group(2).lower()
        n_files += 1
        if channel in by_pos[pos]:
            dupes += 1
            continue          # keep the first; duplicates are a known instrument quirk
        by_pos[pos][channel] = path.name
    logger.info("Indexed %d channel images across %d positions%s",
                n_files, len(by_pos),
                f" ({dupes} duplicate position/channel file(s) ignored)" if dupes else "")
    return by_pos


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    logger = setup_logging()
    args = parse_args()

    run_dir = Path(args.run_dir)
    if not run_dir.is_dir():
        sys.exit(f"ERROR: CellenONE run folder not found: {run_dir}\n"
                 f"Check `cellenone.runs.{args.plate}` in config.yaml — the run folder "
                 f"cannot be guessed from the plate name.")
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.config) as fh:
        cfg_all = yaml.safe_load(fh) or {}
    cfg = cfg_all.get("cellenone", {}) or {}

    # ---- tables ----------------------------------------------------------
    iso_files = sorted(run_dir.glob("Reordered_*_isolated.xls")) or \
        sorted(run_dir.glob("*_isolated.xls"))
    if not iso_files:
        sys.exit(f"ERROR: no *_isolated.xls in {run_dir}")
    iso_rows = read_cellenone_table(iso_files[0])
    logger.info("Read %d rows from %s", len(iso_rows), iso_files[0].name)

    geo_files = sorted(run_dir.glob("Reordered_*_geoprops.xls")) or \
        sorted(run_dir.glob("*_geoprops.xls"))
    geo_rows = read_cellenone_table(geo_files[0]) if geo_files else []
    if geo_files:
        logger.info("Read %d rows from %s", len(geo_rows), geo_files[0].name)

    # ---- per-run criteria: ALWAYS from this run's own file ----------------
    crit_keys = ["EjBound", "SedBound", "DetDiaMinTrans", "DetDiaMaxTrans",
                 "IsoDiaMinTrans", "IsoDiaMaxTrans", "IsoIntMinTrans", "IsoIntMaxTrans",
                 "IsoDiaMinFlu", "IsoDiaMaxFlu", "IsoIntMinFlu", "IsoIntMaxFlu"]
    criteria = {}
    criteria_varied = {}
    for k in crit_keys:
        vals = [as_float(r.get(k)) for r in iso_rows if r.get(k) not in (None, "")]
        vals = [v for v in vals if v is not None]
        if not vals:
            continue
        counts = Counter(vals)
        # The MODAL value, compared numerically. `sorted(["12","9"])[0]` is "12", so a
        # string sort here silently picked a threshold by alphabet — and an operator
        # really can change the criteria part-way through a run (plate22 does: 92 wells
        # at IsoDiaMin 9 um, then 292 at 12, with IsoDiaMax moving 25 -> 30 too).
        criteria[k] = counts.most_common(1)[0][0]
        if len(counts) > 1:
            criteria_varied[k] = {str(v): n for v, n in sorted(counts.items())}
            logger.warning("Criterion %s changed during the run: %s — the per-well value "
                           "is carried in wells_raw.tsv; run_meta records the modal %s",
                           k, criteria_varied[k], criteria[k])
    logger.info("Run criteria: detection >= %s um, isolation window %s-%s um",
                criteria.get("DetDiaMinTrans"),
                criteria.get("IsoDiaMinTrans"), criteria.get("IsoDiaMaxTrans"))

    # `.par` restates these as JSON; keep it as a cross-check, never as a requirement.
    par_params = {}
    par_files = sorted(run_dir.glob("*.par"))
    if par_files:
        try:
            text = par_files[0].read_text(encoding="latin-1", errors="replace")
            m = re.search(r'TParameters\s*=\s*"(\{.*?\})"', text, re.S)
            if m:
                par_params = json.loads(m.group(1).replace('""', '"'))
        except (OSError, json.JSONDecodeError) as exc:
            logger.warning("Could not parse %s: %s", par_files[0].name, exc)

    # ---- rows by position and channel -------------------------------------
    by_pos_channel = defaultdict(dict)
    unparsed = 0
    for r in iso_rows:
        teg = (r.get("Teg") or "").strip().lower()
        channel = CHANNEL_OF_TEG.get(teg)
        if channel is None:
            continue
        try:
            pos = L384.parse_pos(r.get("Well", ""))
        except ValueError:
            unparsed += 1
            continue
        by_pos_channel[pos][channel] = r
    if unparsed:
        logger.warning("%d isolated.xls row(s) had an unparseable Well token", unparsed)

    images = index_images(run_dir, logger)

    # ---- channel alignment offsets ----------------------------------------
    # The three channels share one camera and frame, but there is a small systematic
    # optical offset. It is measurable directly from the table as the median
    # (fluorescence position - transmission position), so no image registration is
    # needed: a per-run integer-pixel shift per channel is enough. A cell is ~24 px
    # across, so the ~8 px Orange offset is a third of a cell — worth correcting.
    offsets = {}
    for channel in ("blue", "green", "orange", "red"):
        dxs, dys = [], []
        for pos, chans in by_pos_channel.items():
            t, f = chans.get("trans"), chans.get(channel)
            if not t or not f:
                continue
            tx, ty = as_float(t.get("X")), as_float(t.get("Y"))
            fx, fy = as_float(f.get("X")), as_float(f.get("Y"))
            if None in (tx, ty, fx, fy):
                continue
            if abs(fx - tx) > 25 or abs(fy - ty) > 25:
                continue      # not the same object; do not let it drag the median
            dxs.append(fx - tx)
            dys.append(fy - ty)
        if dxs:
            dxs.sort(); dys.sort()
            mid = len(dxs) // 2
            offsets[channel] = {"dx": round(dxs[mid], 2), "dy": round(dys[mid], 2),
                                "n": len(dxs)}
            logger.info("Channel %s offset: dx=%+.1f dy=%+.1f px (from %d wells)",
                        channel, dxs[mid], dys[mid], len(dxs))
        else:
            offsets[channel] = {"dx": 0.0, "dy": 0.0, "n": 0}

    # ---- run-level statistics from geoprops -------------------------------
    drops_attempted = len({r.get("DropNo") for r in geo_rows if r.get("DropNo")})
    n_printed = len(by_pos_channel)
    run_stats = {
        "drops_attempted": drops_attempted,
        "wells_isolated": n_printed,
        "isolation_frequency_pct": (round(n_printed / drops_attempted * 100, 2)
                                    if drops_attempted else None),
        "geoprops_rows": len(geo_rows),
    }

    # ---- map CellenONE positions onto sequencing wells --------------------
    subplates = list(args.subplates)
    is_384 = bool(subplates)
    if is_384:
        layout = L384.build_layout(subplates, wells=args.wells,
                                   layout_tsv=args.layout_tsv or None)
    else:
        # A 96-well plate maps straight onto the A1..H12 grid.
        layout = []
        for i, w in enumerate(args.wells):
            layout.append({"id": w, "subplate": args.plate, "sub_index": 0,
                           "well": w, "row": i % 8, "col": i // 8,
                           "pos384": L384.format_pos(i % 8, i // 8),
                           "sample_id": f"{args.plate}_{w}"})

    fieldnames = [
        "id", "subplate", "well", "pos384", "drop_no",
        "x", "y", "diameter_um", "elongation", "circularity", "intensity",
        "blue_intensity", "blue_diameter_um", "green_intensity", "green_diameter_um",
        "orange_intensity", "orange_diameter_um", "red_intensity", "red_diameter_um",
        "img_trans", "img_blue", "img_green", "img_orange", "img_red",
        # This well's OWN gate, not the run summary — see criteria_varied above.
        "det_dia_min", "iso_dia_min", "iso_dia_max",
    ]
    out_rows = []
    matched = 0
    no_cellenone = []
    for e in layout:
        pos = (e["row"], e["col"])
        chans = by_pos_channel.get(pos, {})
        imgs = images.get(pos, {})
        t = chans.get("trans", {})
        blue, orange = chans.get("blue", {}), chans.get("orange", {})
        green, red = chans.get("green", {}), chans.get("red", {})
        if chans or imgs:
            matched += 1
        else:
            no_cellenone.append(e["pos384"])
        out_rows.append({
            "id": e["id"], "subplate": e["subplate"], "well": e["well"],
            "pos384": e["pos384"], "drop_no": t.get("DropNo", ""),
            "x": t.get("X", ""), "y": t.get("Y", ""),
            "diameter_um": t.get("Diameter", ""), "elongation": t.get("Elongation", ""),
            "circularity": t.get("Circularity", ""), "intensity": t.get("Intensity", ""),
            "blue_intensity": blue.get("Intensity", ""),
            "blue_diameter_um": blue.get("Diameter", ""),
            "green_intensity": green.get("Intensity", ""),
            "green_diameter_um": green.get("Diameter", ""),
            "orange_intensity": orange.get("Intensity", ""),
            "orange_diameter_um": orange.get("Diameter", ""),
            "red_intensity": red.get("Intensity", ""),
            "red_diameter_um": red.get("Diameter", ""),
            "img_trans": imgs.get("trans", ""), "img_blue": imgs.get("blue", ""),
            "img_green": imgs.get("green", ""), "img_orange": imgs.get("orange", ""),
            "img_red": imgs.get("red", ""),
            "det_dia_min": t.get("DetDiaMinTrans", ""),
            "iso_dia_min": t.get("IsoDiaMinTrans", ""),
            "iso_dia_max": t.get("IsoDiaMaxTrans", ""),
        })

    # Report both directions of mismatch rather than silently dropping either side.
    seq_positions = {(e["row"], e["col"]) for e in layout}
    orphan_cellenone = sorted(L384.format_pos(*p) for p in by_pos_channel
                              if p not in seq_positions)
    if no_cellenone:
        logger.warning("%d sequencing well(s) have no CellenONE record: %s",
                       len(no_cellenone), ", ".join(no_cellenone[:20])
                       + (" …" if len(no_cellenone) > 20 else ""))
    if orphan_cellenone:
        logger.warning("%d CellenONE position(s) have no sequencing well: %s",
                       len(orphan_cellenone), ", ".join(orphan_cellenone[:20])
                       + (" …" if len(orphan_cellenone) > 20 else ""))

    wells_tsv = outdir / "wells_raw.tsv"
    with open(wells_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(out_rows)
    logger.info("Wrote %s (%d rows, %d with CellenONE data)",
                wells_tsv, len(out_rows), matched)

    meta = {
        "plate": args.plate,
        "run_dir": str(run_dir),
        "run_name": run_dir.name,
        "lines": resolve_lines(run_dir, cfg, logger),
        "criteria": criteria,
        "criteria_varied": criteria_varied,
        "par_tparameters": par_params,
        "channel_offsets": offsets,
        "run_stats": run_stats,
        "counts": {
            "layout_wells": len(layout),
            "matched": matched,
            "missing_cellenone": len(no_cellenone),
            "orphan_cellenone": len(orphan_cellenone),
            "images_indexed": sum(len(v) for v in images.values()),
        },
        "image_dir": str(run_dir),
    }
    meta_path = outdir / "run_meta.json"
    meta_path.write_text(json.dumps(meta, indent=2) + "\n")
    logger.info("Wrote %s", meta_path)
    logger.info("Isolation frequency: %s%% (%d isolated / %d drops attempted)",
                run_stats["isolation_frequency_pct"], n_printed, drops_attempted)


if __name__ == "__main__":
    main()
