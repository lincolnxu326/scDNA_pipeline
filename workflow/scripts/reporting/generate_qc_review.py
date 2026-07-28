#!/usr/bin/env python3
"""
Build the static HTML QC review report for one plate.

This reads the per-well QC metrics already produced by the pipeline
(samtools stats, dedup / filter summaries, demux stats) plus the per-well
AneuFinder profile PNGs rendered by render_well_profiles.R, computes an
advisory read-count gate status from config.yaml `qc_review.usable_reads_*_cutoff`,
and writes a single self-contained `review.html`.

Layout, styling and behaviour live in `assets/qc_review.{css,js}` next to this
script and are pasted into the page verbatim; the markup is `HTML_TEMPLATE` at the
bottom. Python's whole job here is to build the `window.QC` payload (the data
contract in `design/handoff/README.md`) and substitute the `__UPPER_SNAKE__` tokens. The
JS owns the grid, selection, decision state, localStorage autosave, keyboard,
filters, CSV and the chromosome axis — do not reimplement any of that here.

Design constraints (deliberate):
  * Fully static. No server, no database, no external requests, no web fonts.
  * The browser only ever *generates* copyable text, a downloaded CSV, or a
    copyable shell command; it never writes into the project directory. The
    reproducible human artefact is `<PLATE_DIR>/qc_decisions.csv` (per plate),
    which the user saves by hand or via the in-page "COPY SAVE COMMAND".
  * Asset handling depends on plate format:
      96  — every plot is base64-embedded, so `review.html` is a single
            self-contained file that opens correctly over Samba or as file://.
      384 — assets are referenced relatively (`plots/…`, `assets/…`) instead.
            384 wells of embedded plots plus cell images would be a ~107 MB page
            that no browser opens comfortably; referencing them keeps the HTML
            small and makes the browser fetch only the clicked well.

384 mode
--------
A 384-well plate is four interleaved 96-well subplates. This report shows all four
in ONE 16x24 plate map, each well still labelled with the subplate it came from, and
a tab strip (`All | SL1 | SL2 | SL3 | SL4`) that dims wells outside the selection.
Well identity is `<subplate>_<well>` (e.g. `plate21_1_W01`) — identical to what a
standalone run of that subplate produces, so a cell keeps the same `sample_id`
whether it is reviewed alone or as part of the plate.
"""

import argparse
import base64
import csv
import html
import json
import logging
import os
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent))
import plate384_layout as L384

# The stylesheet and behaviour, pasted into every report verbatim. Both are copies of
# `design/handoff/` — diff against that folder before editing either, and change the prototype
# (`QC Review Redesign.dc.html`) first if the design itself needs to move.
ASSETS_DIR = Path(__file__).resolve().parent / "assets"

# Controlled vocabularies (kept in sync with validate_qc_decisions.py and config/README.md).
# The page has its own copy in assets/qc_review.js (DECISIONS / REASON_GROUPS, which also
# fixes how the reasons are grouped in the drawer) — change both together or the CSV a
# reviewer produces will not validate.
DECISIONS = ["PASS", "EXCLUDE", "REVIEW", "REPEAT"]
REASONS = [
    "low_read_count",
    "noisy_profile",
    "poor_bin_distribution",
    "low_complexity",
    "suspected_doublet_or_mixed_well",
    "sample_swap_suspected",
    "manual_exception",
    "other",
    "missing_qc_metric",
]

# Default read-count gating cutoffs (overridden by config qc_review.*)
DEFAULT_PASS_CUTOFF = 100000
DEFAULT_WARN_CUTOFF = 50000


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )
    return logging.getLogger("qc_review")


def parse_args():
    p = argparse.ArgumentParser(description="Generate the static QC review HTML report.")
    p.add_argument("--plate", required=True, help="Plate name")
    p.add_argument("--plate-dir", required=True, help="Per-plate output directory")
    p.add_argument("--multiqc-data", required=False, default="",
                   help=argparse.SUPPRESS)   # DEPRECATED: accepted but never read.
    # Kept only so an existing caller does not break. The pipeline no longer passes
    # it: the value was never used, and the plate-level `<plate>/multiqc/multiqc_data`
    # it named does not exist for a 384 plate. The header links MultiQC relatively
    # instead — see the `multiqc_href` block in main().
    p.add_argument("--plots-dir", required=True, help="Directory of per-well {well}.png")
    p.add_argument("--config", required=True, help="Pipeline config.yaml (gate cutoffs, binsize)")
    p.add_argument("--outdir", required=True, help="Output directory (writes the HTML report)")
    p.add_argument("--wells", required=True, nargs="+", help="Ordered list of well IDs")
    p.add_argument("--decisions-path", default="",
                   help="Destination path for the saved qc_decisions.csv "
                        "(default: <plate-dir>/qc_decisions.csv); used for the in-page save command")
    p.add_argument("--report-kind", choices=["review", "cn"], default="review",
                   help="'review' = first-pass decisions UI; 'cn' = read-only final CN viewer")
    p.add_argument("--included-wells", default="",
                   help="TSV of included wells (cn mode restricts to these)")
    p.add_argument("--heatmap", default="", help="Genome-wide CN heatmap PNG to embed (cn mode)")
    p.add_argument("--title", default="", help="Override report title")
    p.add_argument("--out-name", default="", help="Output HTML filename (default depends on kind)")
    # --- 384 mode (all additive; with --subplates empty every path below is unchanged)
    p.add_argument("--subplates", nargs="*", default=[],
                   help="Subplate names in SL order. Given => 384 mode.")
    p.add_argument("--layout-tsv", default="",
                   help="Optional 384 layout override TSV (subplate, well, pos384)")
    p.add_argument("--cellenone-dir", default="",
                   help="<plate>/cellenone directory (enables the cell-image panel)")
    p.add_argument("--assets-mode", choices=["embed", "sidecar"], default="embed",
                   help="embed = base64 data URIs (96); sidecar = relative refs (384)")
    p.add_argument("--assets-dir", default="",
                   help="Sidecar directory for assets that live outside --outdir "
                        "(default: <outdir>/assets)")
    return p.parse_args()


def load_included_wells(path: str) -> list:
    """Read an included_wells.tsv and return the well IDENTITIES it lists.

    Returns `<subplate>_<well>` when the file carries a subplate column (384 mode)
    and plain `<well>` otherwise, so the result always matches `wells_order`.
    """
    if not path:
        return []
    p = Path(path)
    if not p.exists():
        return []
    out = []
    with open(p, newline="") as fh:
        header = [h.strip() for h in fh.readline().rstrip("\r\n").split("\t")]
        try:
            widx = header.index("well")
        except ValueError:
            widx = 1 if len(header) > 1 else 0
        sidx = header.index("subplate") if "subplate" in header else None
        for line in fh:
            parts = [c.strip() for c in line.rstrip("\r\n").split("\t")]
            if len(parts) <= widx or not parts[widx]:
                continue
            if sidx is not None and len(parts) > sidx and parts[sidx]:
                out.append(f"{parts[sidx]}_{parts[widx]}")
            else:
                out.append(parts[widx])
    return out


def load_cellenone(cellenone_dir: str):
    """Load the CellenONE ingest/render outputs, if present.

    Returns (run_meta, {id: well_row}, {id: [object_row, …]}). Everything here is
    DISPLAY-ONLY — see the note on compute_status().
    """
    if not cellenone_dir:
        return {}, {}, {}
    cdir = Path(cellenone_dir)
    meta = {}
    meta_path = cdir / "run_meta.json"
    if meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text())
        except (json.JSONDecodeError, OSError):
            meta = {}

    wells = {}
    wells_path = cdir / "cellenone_wells.tsv"
    if wells_path.exists():
        with open(wells_path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                wells[row["id"]] = row

    objects = defaultdict(list)
    obj_path = cdir / "objects.tsv"
    if obj_path.exists():
        with open(obj_path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                objects[row["id"]].append(row)
    return meta, wells, objects


# ---------------------------------------------------------------------------
# Metric parsing
# ---------------------------------------------------------------------------

def parse_samtools_stats(stats_file: Path) -> dict:
    """Parse the SN block of a samtools stats file (post-dedup BAM), plus mean GC%.

    GC comes free: samtools already writes GCF/GCL — the GC-content *distribution*
    of first/last fragments as (gc_percent, read_count) pairs — so the per-well mean
    is a weighted average over those rows. No extra tool, no re-run.
    """
    out = {"reads": None, "mapped_reads": None, "average_quality": None,
           "mapping_rate": None, "gc_content": None}
    if not stats_file.exists():
        return out
    total = mapped = None
    gc_num = gc_den = 0.0
    with open(stats_file) as fh:
        for line in fh:
            if line.startswith("SN\traw total sequences:"):
                total = int(line.split("\t")[2])
            elif line.startswith("SN\treads mapped:"):
                mapped = int(line.split("\t")[2])
            elif line.startswith("SN\taverage quality:"):
                out["average_quality"] = float(line.split("\t")[2])
            elif line.startswith("GCF\t") or line.startswith("GCL\t"):
                parts = line.split("\t")
                try:
                    gc_pct, n = float(parts[1]), float(parts[2])
                except (IndexError, ValueError):
                    continue
                gc_num += gc_pct * n
                gc_den += n
    out["reads"] = total
    out["mapped_reads"] = mapped
    if total and total > 0 and mapped is not None:
        out["mapping_rate"] = round(mapped / total * 100, 2)
    if gc_den > 0:
        out["gc_content"] = round(gc_num / gc_den, 2)
    return out


def load_tsv_by_well(path: Path, value_cols: dict) -> dict:
    """Load a well-keyed TSV into {well: {alias: value}} for the given columns."""
    result = {}
    if not path.exists():
        return result
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            idx = {c: header.index(c) for c in ["well_id", *value_cols.keys()]}
        except ValueError:
            return result
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx.values()):
                continue
            well = parts[idx["well_id"]]
            row = {}
            for col, alias in value_cols.items():
                raw = parts[idx[col]]
                try:
                    row[alias] = float(raw)
                except (ValueError, TypeError):
                    row[alias] = None
            result[well] = row
    return result


def load_demux_totals(path: Path) -> dict:
    """Return {well: demuxed_total_reads} from demux_stats.json."""
    if not path.exists():
        return {}
    try:
        data = json.load(open(path))
    except (json.JSONDecodeError, OSError):
        return {}
    wells = data.get("wells", {})
    return {w: v.get("total") for w, v in wells.items() if w != "unassigned"}


# ---------------------------------------------------------------------------
# Automated status — read-count-only gating
# ---------------------------------------------------------------------------

# auto_status -> the gate status the report renders. The UI has three gate colours,
# not four (spec §1), so a well with no usable_reads at all is shown as FAIL: it has
# no measurement to pass on, and its `flags` still say why. The reviewer then sees the
# FAIL default (EXCLUDE), so the well is visibly written off rather than quietly.
GATE_STATUS = {"PASS": "PASS", "WARN": "WARN", "FAIL": "FAIL", "UNKNOWN": "FAIL"}


def compute_status(usable_reads, pass_cutoff: int, warn_cutoff: int) -> tuple:
    """
    Read-count-only gating. Returns (auto_status, default_decision, default_reason, flags).

    usable_reads = mapped reads in the dedup BAM. Duplication is NOT used for gating.
    Advisory only — never overrides the human call.

    Only auto_status (through GATE_STATUS) and flags reach the report: the page derives
    the pre-selected decision from the gate status itself, so default_decision and
    default_reason are here for callers and tests, not for the UI.

    DO NOT reintroduce duplication gating here, and DO NOT add the CellenONE image
    call. The image call is display-only: it gets its own labelled panel and at most
    offers a *suggested* reason chip the reviewer clicks. Letting it move auto_status
    or the default decision would silently re-gate wells on a heuristic the operator
    never opted into.
    """
    if usable_reads is None:
        return ("UNKNOWN", "REVIEW", "missing_qc_metric", ["missing usable_reads"])
    if usable_reads >= pass_cutoff:
        return ("PASS", "PASS", "", [])
    if usable_reads >= warn_cutoff:
        return ("WARN", "REVIEW", "", [f"usable_reads<{pass_cutoff}"])
    return ("FAIL", "EXCLUDE", "low_read_count", [f"usable_reads<{warn_cutoff}"])


# ---------------------------------------------------------------------------
# Plate-map geometry
# ---------------------------------------------------------------------------

def well_positions(wells: list) -> dict:
    """
    Map each well to its PHYSICAL (row, col) on the 8x12 grid.

    Three shapes, tried in order:
      A1..H12   the id *is* the position.
      W01..W96  the dispense is column-major down the 8 rows, so
                row = (n-1) % 8, col = (n-1) // 8 — the same geometry
                plate384_layout applies inside one subplate.
      anything else, or a list that will not fit 8x12: a sequential
                row-major fill, 12 columns, as many rows as needed.

    The old code used the sequential fill for W-numbered wells too, which
    transposed the plate: a bad physical column read as a bad row (spec §0).
    """
    import re
    pos = {}
    positional = True
    for w in wells:
        m = re.match(r"^([A-Ha-h])0*([0-9]+)$", w)
        if not m:
            positional = False
            break
        row = ord(m.group(1).upper()) - ord("A")
        col = int(m.group(2)) - 1
        if not (0 <= row < 8 and 0 <= col < 12):
            positional = False
            break
        pos[w] = (row, col)

    if positional and pos:
        return {"mode": "positional", "rows": 8, "cols": 12, "pos": pos}

    try:
        pos = {}
        for w in wells:
            n0 = L384.parse_well(w) - 1
            pos[w] = (n0 % L384.SUB_ROWS, n0 // L384.SUB_ROWS)
        if pos:
            return {"mode": "physical", "rows": 8, "cols": 12, "pos": pos}
    except ValueError:
        pass

    # Sequential fill: 12 columns, as many rows as needed.
    ncols = 12
    nrows = max(1, -(-len(wells) // ncols))  # ceil
    pos = {w: (i // ncols, i % ncols) for i, w in enumerate(wells)}
    return {"mode": "sequential", "rows": nrows, "cols": ncols, "pos": pos}


def well_positions_384(layout: list) -> dict:
    """The 16x24 plate map, keyed by well identity rather than well id.

    buildMap() in the template needs no structural change for this:
    `String.fromCharCode(65 + r)` already yields A..P and `c + 1` yields 1..24.
    """
    return {"mode": "384", "rows": L384.N_ROWS, "cols": L384.N_COLS,
            "pos": {e["id"]: (e["row"], e["col"]) for e in layout}}


# ---------------------------------------------------------------------------
# Assets — embedded data URIs (96) or relative references (384)
# ---------------------------------------------------------------------------

MIME_BY_SUFFIX = {".png": "image/png", ".svg": "image/svg+xml",
                  ".jpg": "image/jpeg", ".jpeg": "image/jpeg", ".gif": "image/gif"}


def data_uri(path: Path) -> str:
    mime = MIME_BY_SUFFIX.get(path.suffix.lower(), "application/octet-stream")
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"


class AssetRefs:
    """Turn a file on disk into something the HTML can point at.

    embed   -> a base64 data URI (one self-contained file; the 96-well behaviour).
    sidecar -> a relative path that NEVER points above the report directory.

               Assets already inside `<outdir>` (the report's own `plots/`) are
               referenced where they lie. Everything else — chiefly the plate-level
               `cellenone/images/` next door — is brought into `<outdir>/assets/`.

               `../cellenone/images/…` was tried and is wrong: the report is read over
               a Samba share and is routinely copied around, and the moment
               `qc_review/` travels without its parent every cell image 404s. The rule
               is simply that a report must be openable from its own directory.

               To keep that cheap, assets are HARD LINKED rather than copied — at 384
               with every channel on that is ~1,540 files, and copying meant a second
               13 MB and minutes of `copy2` on a busy shared filesystem. A hard link is
               instant, adds no bytes, and (unlike a symlink) is an ordinary file to
               Samba and to anything that later copies the tree. Falls back to a real
               copy when the link cannot be made (different filesystem, or an fs
               without hard links).
    """

    def __init__(self, mode: str, outdir: Path, assets_dir: Path = None,
                 base_dir: Path = None):
        self.mode = mode
        self.outdir = Path(outdir).resolve()
        # Retained for callers; only `outdir` decides what is referenced in place.
        self.base_dir = Path(base_dir).resolve() if base_dir else self.outdir
        self.assets_dir = Path(assets_dir) if assets_dir else self.outdir / "assets"
        self.n_embedded = 0
        self.n_linked = 0
        self.n_copied = 0
        self.n_hardlinked = 0
        self.bytes_embedded = 0

    def ref(self, path) -> str:
        if not path:
            return ""
        path = Path(path)
        if not path.exists():
            return ""
        if self.mode == "embed":
            self.n_embedded += 1
            self.bytes_embedded += path.stat().st_size
            return data_uri(path)
        resolved = path.resolve()
        try:
            rel = resolved.relative_to(self.outdir)
        except ValueError:
            pass                     # outside the report dir -> bring it inside
        else:
            self.n_linked += 1
            return str(rel).replace(os.sep, "/")

        self.assets_dir.mkdir(parents=True, exist_ok=True)
        dest = self.assets_dir / path.name
        stale = (not dest.exists()) or dest.stat().st_mtime < resolved.stat().st_mtime
        if stale:
            if dest.exists():
                dest.unlink()
            try:
                os.link(resolved, dest)          # same fs: instant, no extra bytes
                self.n_hardlinked += 1
            except OSError:
                shutil.copy2(resolved, dest)     # cross-device or no hardlink support
                self.n_copied += 1
        self.n_linked += 1
        rel_dir = os.path.relpath(self.assets_dir, self.outdir).replace(os.sep, "/")
        return f"{rel_dir}/{path.name}"

    def summary(self) -> str:
        if self.mode == "embed":
            return f"embedded {self.n_embedded} asset(s), {self.bytes_embedded / 1e6:.1f} MB"
        return (f"referenced {self.n_linked} asset(s) "
                f"({self.n_hardlinked} hard-linked, {self.n_copied} copied "
                f"into {self.assets_dir.name}/)")


def well_plot_path(plots_dir: Path, well: str, suffix: str):
    """Locate {well}_{suffix}, preferring SVG (vector) then PNG.

    Falls back to the legacy {well}.png (profile only). Returns None if absent.
    """
    svg = plots_dir / f"{well}_{suffix}.svg"
    if svg.exists():
        return svg
    png = plots_dir / f"{well}_{suffix}.png"
    if png.exists():
        return png
    if suffix == "profile":
        legacy = plots_dir / f"{well}.png"
        if legacy.exists():
            return legacy
    return None


# ---------------------------------------------------------------------------
# Rendering — build the payload, substitute the tokens, paste in the CSS/JS
# ---------------------------------------------------------------------------

# The usable-reads bar in the inspector is a FIXED scale (spec §4), so every well is
# read against the same ruler; only the two gate markers move with the config.
USABLE_SCALE = 500000


def read_label(n) -> str:
    """50000 -> '50k', 100000 -> '100k', 1500000 -> '1.5M'. Used for the gate labels."""
    if n is None:
        return "?"
    if n >= 1e6:
        return f"{n / 1e6:g}M"
    if n >= 1e3:
        return f"{n / 1e3:g}k"
    return str(int(n))


def read_label_bp(n) -> str:
    """1000000 -> '1 Mb', 500000 -> '500 kb'. The profile caption's bin size."""
    if not n:
        return ""
    if n >= 1e6:
        return f"{n / 1e6:g} Mb"
    if n >= 1e3:
        return f"{n / 1e3:g} kb"
    return f"{int(n)} bp"


def bar_pct(n) -> str:
    """Where a cutoff sits on the fixed 0-500k usable-reads bar, as a CSS length."""
    return f"{min(100.0, max(0.0, 100.0 * (n or 0) / USABLE_SCALE)):g}%"


def as_number(raw, cast):
    """TSV cells are strings and may be empty. Absent stays None, never 0."""
    if raw is None or raw == "":
        return None
    try:
        return cast(float(raw))
    except (TypeError, ValueError):
        return None


def read_asset(name: str) -> str:
    path = ASSETS_DIR / name
    if not path.exists():
        raise FileNotFoundError(
            f"missing report asset {path}. The CSS and JS are pasted into the report "
            f"verbatim; without them there is no report."
        )
    return path.read_text(encoding="utf-8")


def render(template: str, payload: dict, tokens: dict) -> str:
    """Substitute the `__UPPER_SNAKE__` tokens.

    str.replace, never %-formatting or .format(): the CSS is full of braces. The CSS,
    JS and JSON go in LAST so that nothing already substituted can be re-scanned, and
    `</` is escaped inside the JSON so a stray string can never close the <script>.
    """
    import re
    blob = json.dumps(payload, separators=(",", ":")).replace("</", "<\\/")
    late = {"__QC_CSS__": read_asset("qc_review.css"),
            "__QC_JS__": read_asset("qc_review.js"),
            "__QC_JSON__": blob}
    missing = sorted(set(re.findall(r"__[A-Z0-9_]+__", template)) - set(tokens) - set(late))
    if missing:
        raise AssertionError(f"template token(s) with nothing to substitute: {missing}")

    out = template
    for token, value in sorted(tokens.items(), key=lambda kv: -len(kv[0])):
        out = out.replace(token, value)
    for token, value in late.items():
        out = out.replace(token, value)
    return out


def main():
    logger = setup_logging()
    args = parse_args()

    plate_dir = Path(args.plate_dir)
    plots_dir = Path(args.plots_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)
    thr = cfg.get("qc_thresholds", {}) or {}
    review_cfg = cfg.get("qc_review", {}) or {}
    pass_cutoff = review_cfg.get("usable_reads_pass_cutoff", DEFAULT_PASS_CUTOFF)
    warn_cutoff = review_cfg.get("usable_reads_warn_cutoff", DEFAULT_WARN_CUTOFF)

    kind = args.report_kind

    # ---- plate layout -------------------------------------------------------
    subplates = list(args.subplates or [])
    is_384 = bool(subplates)
    if is_384:
        layout = L384.build_layout(subplates, wells=args.wells,
                                   layout_tsv=args.layout_tsv or None)
        geom = well_positions_384(layout)
    else:
        layout = [{"id": w, "subplate": args.plate, "sub_index": 0,
                   "sub_label": args.plate, "well": w, "pos384": "",
                   "sample_id": f"{args.plate}_{w}"} for w in args.wells]
        geom = well_positions(args.wells)

    ids = [e["id"] for e in layout]   # always the FULL plate (cn greys out excluded)
    included_set = set()
    if kind == "cn":
        included = load_included_wells(args.included_wells)
        if args.included_wells and not included:
            logger.warning("No included wells found in %s; CN report will mark all as excluded",
                           args.included_wells)
        included_set = set(included)
        # The cn viewer shows the WHOLE plate: included (PASS) wells keep their original
        # auto-status colour and carry second-pass plots; excluded wells are greyed out.
    logger.info("Building %s report for plate %s (%d wells%s, %d included)",
                kind, args.plate, len(ids),
                f", {len(subplates)} subplates" if is_384 else "",
                len(included_set) if kind == "cn" else len(ids))

    assets = AssetRefs(args.assets_mode, outdir,
                       Path(args.assets_dir) if args.assets_dir else None,
                       base_dir=plate_dir)

    # ---- per-subplate metric sources ---------------------------------------
    # The preprocessing outputs stay under each subplate directory in 384 mode; the
    # existing loaders are reused verbatim, just pointed at the right directory.
    sources = {}
    for sub in (subplates or [args.plate]):
        sdir = (plate_dir / sub) if is_384 else plate_dir
        sources[sub] = {
            "dir": sdir,
            "dedup": load_tsv_by_well(
                sdir / "dedup" / "dedup_summary.tsv",
                {"unique_reads": "unique_reads", "duplicate_reads": "duplicate_reads",
                 "dedup_rate": "dedup_rate"}),
            "dimers": load_tsv_by_well(
                sdir / "filtered" / "adapter_filter_summary.tsv",
                {"dimer_rate": "dimer_rate", "kept_pairs": "kept_pairs"}),
            "demux": load_demux_totals(sdir / "demux" / "demux_stats.json"),
        }

    cell_meta, cell_wells, cell_objects = load_cellenone(args.cellenone_dir)
    if args.cellenone_dir:
        logger.info("CellenONE: %d well record(s), %d with objects",
                    len(cell_wells), len(cell_objects))

    wells_payload = []
    n_profiles = 0
    n_hist = 0
    n_cell_images = 0
    status_counts = {"PASS": 0, "WARN": 0, "FAIL": 0, "UNKNOWN": 0}
    for e in layout:
        wid, sub, w = e["id"], e["subplate"], e["well"]
        src = sources[sub]
        sam = parse_samtools_stats(src["dir"] / "bam" / f"{w}.stats.txt")
        metrics = dict(sam)
        # The alignment rate has to come from the RAW (pre-dedup) BAM. umi_tools drops
        # unmapped reads, so in the dedup BAM `reads mapped` == `raw total sequences`
        # and mapped/total is 100.00% for every well by construction — which is what
        # this field used to report. The real figure is ~98.5%.
        raw = parse_samtools_stats(src["dir"] / "raw_bam" / f"{w}.stats.txt")
        metrics["mapping_rate"] = raw.get("mapping_rate")
        metrics["aligned_reads"] = raw.get("mapped_reads")
        if w in src["dedup"]:
            metrics["unique_reads"] = src["dedup"][w].get("unique_reads")
            metrics["duplicate_reads"] = src["dedup"][w].get("duplicate_reads")
            metrics["dedup_rate"] = src["dedup"][w].get("dedup_rate")
        if w in src["dimers"]:
            metrics["dimer_rate"] = src["dimers"][w].get("dimer_rate")
            metrics["kept_pairs"] = src["dimers"][w].get("kept_pairs")
        if w in src["demux"]:
            metrics["demux_reads"] = src["demux"][w]

        # UMI dedup retention (diagnostic only): unique / (unique + duplicate)
        uniq = metrics.get("unique_reads")
        dup = metrics.get("duplicate_reads")
        if uniq is not None and dup is not None and (uniq + dup) > 0:
            metrics["umi_retention"] = round(uniq / (uniq + dup) * 100, 2)

        # Read-count-only gating on usable_reads = mapped reads in the dedup BAM.
        # NOTE: the CellenONE image call is DELIBERATELY not an input here — see the
        # comment on compute_status().
        usable_reads = metrics.get("mapped_reads")
        metrics["usable_reads"] = usable_reads
        auto_status, default_decision, default_reason, flags = compute_status(
            usable_reads, pass_cutoff, warn_cutoff)
        status_counts[auto_status] = status_counts.get(auto_status, 0) + 1

        # plots (prefer SVG → PNG). In cn-mode only included wells have second-pass plots.
        is_included = (kind != "cn") or (wid in included_set)
        profile_uri = assets.ref(well_plot_path(plots_dir, wid, "profile")) if is_included else ""
        hist_uri = assets.ref(well_plot_path(plots_dir, wid, "histogram")) if is_included else ""
        if profile_uri:
            n_profiles += 1
        if hist_uri:
            n_hist += 1

        # CellenONE evidence. DISPLAY ONLY — see the note on compute_status().
        # The counts come straight from render_cellenone_images.py, which derives them
        # and the call from the same object list, so they agree by construction (spec §4).
        cell_images = {}
        call, n_obj, n_iso, right_x, right_dia = "NO_IMAGE", None, None, None, None
        cell_dia = cell_elong = cell_circ = cell_int = None
        flu_int = {}
        if wid in cell_wells:
            cw = cell_wells[wid]
            for name in (cw.get("image_paths") or "").split(";"):
                if not name:
                    continue
                channel = name.rsplit("_", 1)[-1].rsplit(".", 1)[0]
                if channel not in ("merge", "trans", "blue", "orange", "red"):
                    continue
                uri = assets.ref(Path(args.cellenone_dir) / "images" / name)
                if uri:
                    cell_images[channel] = uri
            if cell_images:
                n_cell_images += 1
            call = cw.get("image_call") or "NO_IMAGE"
            n_obj = as_number(cw.get("n_objects"), int)
            n_iso = as_number(cw.get("n_in_iso_window"), int)
            right_x = as_number(cw.get("rightmost_x"), float)
            right_dia = as_number(cw.get("rightmost_diameter_um"), float)
            # CellenONE's own measurements of the isolated cell (secondary column).
            cell_dia = as_number(cw.get("diameter_um"), float)
            cell_elong = as_number(cw.get("elongation"), float)
            cell_circ = as_number(cw.get("circularity"), float)
            cell_int = as_number(cw.get("intensity"), float)
            for _ch in ("blue", "orange", "red"):
                _v = as_number(cw.get(f"{_ch}_intensity"), float)
                # 0 means the channel recorded nothing for this cell — report it as
                # absent rather than as a measured zero. (Red was not a configured
                # channel on this run at all, so it is 0 for every well that has a
                # stray Red frame.)
                if _v:
                    flu_int[_ch] = _v

        r, c = geom["pos"][wid]
        wells_payload.append({
            "n": len(wells_payload),
            "id": e.get("sample_id") or f"{args.plate}_{w}",   # the CSV key
            "well": w,
            "subplate": sub if is_384 else "",
            "pos": e.get("pos384", "") if is_384 else "",
            "row": r,
            "col": c,
            "status": GATE_STATUS[auto_status],
            "reads": usable_reads,
            "included": is_included,
            "call": call,
            "nObj": n_obj,
            "nIso": n_iso,
            "rightmost_x": right_x,
            "rightmost_dia": right_dia,
            # Shape/intensity of the cell CellenONE isolated. Secondary detail — the
            # panel keeps these behind a disclosure, so they never compete with the
            # four primary image facts above.
            "cell_dia": cell_dia,
            "cell_elong": cell_elong,
            "cell_circ": cell_circ,
            "cell_int": cell_int,
            # Per-channel fluorescence intensity of the isolated cell, keyed by the
            # instrument's channel name. Shown against whichever channel is on screen.
            "flu": flu_int or None,
            "flags": flags,
            # An absent metric is null, never 0: the UI renders "–" for null and would
            # otherwise claim a real measurement of zero.
            "metrics": {k: metrics.get(k) for k in (
                "reads", "demux_reads", "unique_reads", "mapping_rate",
                "dedup_rate", "umi_retention", "gc_content", "dimer_rate",
                "average_quality")},
            "images": cell_images or None,
            "profile_src": profile_uri or None,
            "hist_src": hist_uri or None,
        })

    logger.info("Resolved %d/%d profile plots, %d/%d histogram plots",
                n_profiles, len(ids), n_hist, len(ids))
    # No-cell-image mode. A plate is only "imaged" if an image actually resolved for at
    # least one well — not merely because --cellenone-dir was passed. That way a run
    # folder that is present but yields nothing degrades the same way as a plate that
    # never had a CellenONE at all, instead of shipping 384 empty image frames.
    has_cell_images = bool(args.cellenone_dir) and n_cell_images > 0
    if args.cellenone_dir:
        logger.info("Cell images for %d/%d wells", n_cell_images, len(ids))
        if not has_cell_images:
            logger.warning("--cellenone-dir %s resolved no images; building the report "
                           "without the cell-image panel", args.cellenone_dir)
    else:
        logger.info("No cell images for this plate: Evidence is the copy-number "
                    "profile and the bin-count histogram only")
    logger.info("Auto-status: PASS=%d WARN=%d FAIL=%d UNKNOWN=%d",
                status_counts["PASS"], status_counts["WARN"],
                status_counts["FAIL"], status_counts["UNKNOWN"])

    heatmap_uri = ""
    if kind == "cn" and args.heatmap:
        heatmap_uri = assets.ref(Path(args.heatmap))
        if heatmap_uri:
            logger.info("Genome heatmap: %s", args.heatmap)

    if args.title:
        title = args.title
    else:
        title = (f"{args.plate} final copy number" if kind == "cn"
                 else f"{args.plate} QC review")

    # MultiQC lives per subplate, but the header has ONE link slot, so 384 points at SL1
    # and says how many there are; the other three sit beside it in the same directory.
    if is_384:
        multiqc_href = f"../{subplates[0]}/multiqc/multiqc_report.html"
        multiqc_label = f"{len(subplates)} subplates"
    else:
        multiqc_href = "../multiqc/multiqc_report.html"
        multiqc_label = args.plate

    # The gate labels and the two markers on the usable-reads bar come from the
    # configured cutoffs, not from literals, so the report cannot claim a gate it did
    # not apply. At the defaults (50k/100k) this reproduces the prototype exactly.
    warn_label, pass_label = read_label(warn_cutoff), read_label(pass_cutoff)
    if pass_cutoff > USABLE_SCALE:
        logger.warning("usable_reads_pass_cutoff (%d) is off the fixed 0-%s bar scale; "
                       "its marker is pinned at 100%%", pass_cutoff, read_label(USABLE_SCALE))

    subplate_tabs = "".join(
        f'<button type="button" class="tab" data-sub="{html.escape(s, quote=True)}">SL{i + 1}</button>'
        for i, s in enumerate(subplates)
    )

    decisions_path = args.decisions_path or str(Path(args.plate_dir) / "qc_decisions.csv")
    bins_label = read_label_bp((cfg.get("aneufinder", {}) or {}).get("binsize"))

    # The data contract — design/handoff/README.md § Data contract. Everything the page does
    # with it (grid, decisions, autosave, CSV, keyboard) belongs to assets/qc_review.js.
    payload = {
        "plate": args.plate,
        "size": 384 if is_384 else 96,
        "mode": kind,
        "gate": {"warn": warn_cutoff, "pass": pass_cutoff},
        "bins_label": bins_label,
        # False => qc_review.js drops the cell-image row and the droplet call.
        "cell_images": has_cell_images,
        "decisions_path": decisions_path,
        "decisions_dir": str(Path(decisions_path).parent),
        "heatmap": heatmap_uri,
        "wells": wells_payload,
    }
    tokens = {
        "__TITLE__": html.escape(title),
        "__WELL_COUNT__": str(len(wells_payload)),
        "__GATE_SUMMARY__": f"&#8805;{warn_label} warn &#183; &#8805;{pass_label} pass",
        "__PASS_NOTE__": f"&#8805; {pass_label} usable reads",
        "__REVIEW_NOTE__": f"{warn_label} &#8211; {pass_label} reads",
        "__FAIL_NOTE__": f"&lt; {warn_label} reads",
        "__WARN_PCT__": bar_pct(warn_cutoff),
        "__PASS_PCT__": bar_pct(pass_cutoff),
        "__WARN_TICK__": warn_label,
        "__PASS_TICK__": pass_label,
        "__MULTIQC_HREF__": html.escape(multiqc_href, quote=True),
        "__MULTIQC_LABEL__": html.escape(multiqc_label),
        "__DECISIONS_PATH__": html.escape(decisions_path),
        "__SUBPLATE_TABS__": subplate_tabs,
    }

    page = render(HTML_TEMPLATE, payload, tokens)
    out_name = args.out_name or ("cn_review.html" if kind == "cn" else "review.html")
    out_html = outdir / out_name
    out_html.write_text(page, encoding="utf-8")
    logger.info("Wrote %s (%.1f MB HTML; %s)",
                out_html, out_html.stat().st_size / 1e6, assets.summary())


# ---------------------------------------------------------------------------
# The page itself. Markup from design/handoff/qc_review_template.html; the CSS and JS are
# read verbatim from assets/ and pasted in, so the report is one self-contained file.
# Substitution is plain str.replace — %-formatting and .format() both choke on the CSS.
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<!-- QC review report — generated by workflow/scripts/reporting/generate_qc_review.py.
     Markup: design/handoff/qc_review_template.html. Placeholders are upper-snake tokens between
     double underscores, substituted with str.replace (NOT %-formatting or .format():
     the CSS is full of braces). render() fails the build on an unknown one.
     __QC_CSS__ / __QC_JS__ are the verbatim contents of assets/qc_review.{css,js} next to
     this script, so the report stays one self-contained file that opens from file://. -->
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>__TITLE__</title>
<style>
__QC_CSS__
</style>
</head>
<body>
<div class="page">

  <header class="hdr">
    <h1>__TITLE__</h1>
    <div class="hdr-facts">
      <div><div class="micro">Wells</div><div class="v">__WELL_COUNT__</div></div>
      <div><div class="micro">Read gate</div><div class="v">__GATE_SUMMARY__</div></div>
      <div><div class="micro">MultiQC</div><div class="v"><a href="__MULTIQC_HREF__">__MULTIQC_LABEL__ &#8599;</a></div></div>
      <div><div class="micro">Needing you</div><div class="v mono" id="progress">&#8211;</div></div>
    </div>
    <input class="search" id="search" placeholder="well, position, sample id&#8230;">
    <!-- cn only; qc_review.js unhides it and hides the search box -->
    <div class="hdr-cn" id="hdr-cn" hidden>
      <div class="micro">Cells in second pass</div>
      <div class="n" id="included-count">&#8211;</div>
    </div>
  </header>

  <!-- ============ cn only: genome-wide copy number ============ -->
  <section class="bp cnmap" id="cnmap" hidden>
    <i class="corner tl"></i><i class="corner tr"></i><i class="corner bl"></i><i class="corner br"></i>
    <div class="cnmap-head">
      <h2>Genome-wide copy number &#183; <span id="cn-cells">&#8211;</span> cells</h2>
      <span class="note">natural size &#183; scroll horizontally</span>
    </div>
    <div class="cnmap-scroll"><img id="cn-heatmap" alt="genome-wide copy-number heatmap"></div>
    <div class="cnmap-legend">
      <span>copy number</span>
      <span><i style="background:#8c2b22"></i>0</span>
      <span><i style="background:#d5928a"></i>1</span>
      <span><i style="background:#e8e6e2"></i>2</span>
      <span><i style="background:#a3c1dc"></i>3</span>
      <span><i style="background:#4a7fae"></i>4</span>
      <span><i style="background:#22496f"></i>5+</span>
    </div>
  </section>

  <div class="cols">

    <!-- ============ left: plate map ============ -->
    <div class="col-map">
      <div class="bp map-panel">
        <i class="corner tl"></i><i class="corner tr"></i><i class="corner bl"></i><i class="corner br"></i>

        <div class="map-head">
          <h2>Plate map</h2>
          <div id="subtabs" style="display:flex;gap:8px">
            <button type="button" class="tab on">ALL</button>
            __SUBPLATE_TABS__
          </div>
        </div>

        <div class="stats">
          <div class="stat stat-pass" id="stat-pass">
            <div class="top"></div><div class="micro">Pass rate</div>
            <div class="row"><span class="big">&#8211;</span><span class="sub"></span></div>
            <div class="track"><i></i></div>
            <div class="note">__PASS_NOTE__</div>
          </div>
          <div class="stat stat-review" id="stat-review">
            <div class="top"></div><div class="micro">Review rate</div>
            <div class="row"><span class="big">&#8211;</span><span class="sub"></span></div>
            <div class="track"><i></i></div>
            <div class="note">__REVIEW_NOTE__</div>
          </div>
          <div class="stat stat-fail" id="stat-fail">
            <div class="top"></div><div class="micro">Fail rate</div>
            <div class="row"><span class="big">&#8211;</span><span class="sub"></span></div>
            <div class="track"><i></i></div>
            <div class="note">__FAIL_NOTE__</div>
          </div>
        </div>

        <div class="grid-cols"><div style="flex:0 0 18px"></div><div class="colhdr" id="colhdr"></div></div>
        <div class="grid-wrap">
          <div class="rowhdr" id="rowhdr"></div>
          <div class="wells" id="wells"></div>
        </div>

        <div class="hoverline mono" id="hoverline"></div>

        <div class="legend">
          <span><i style="background:#2e7d4f"></i>PASS</span>
          <span><i style="background:#d9e7f7;border:1px solid #2f6cad"></i>REVIEW</span>
          <span><i style="background:#a8711c"></i>REPEAT</span>
          <span><i style="background:#f7dfdc;border:1px solid #b23a2f"></i>EXCLUDE</span>
          <span><i style="border:1px dotted rgba(29,31,32,.45)"></i>decided &#9632; / seen &#9633;</span>
        </div>
      </div>

      <div class="export" id="export">
        <span class="micro" style="color:rgba(231,231,234,.55)">Export</span>
        <span class="sum" id="export-summary"></span>
        <button type="button" class="primary" id="copy-cmd">COPY SAVE COMMAND</button>
        <button type="button" class="ghost" id="copy-csv">COPY CSV</button>
        <button type="button" class="ghost" id="download-csv">DOWNLOAD CSV</button>
        <span class="flash" id="flash"></span>
        <span class="path">__DECISIONS_PATH__</span>
      </div>
    </div>

    <!-- ============ right: inspector ============ -->
    <div class="col-insp">
      <div class="bp insp">
        <i class="corner tl"></i><i class="corner tr"></i><i class="corner bl"></i><i class="corner br"></i>

        <div class="insp-head">
          <div>
            <h2 id="well-id">&#8211;</h2>
            <div class="chips">
              <span class="chip" id="chip-decision"></span>
              <span class="chip" id="chip-source"></span>
              <span class="chip c-pos" id="chip-pos"></span>
            </div>
          </div>
          <div class="dec" id="dec">
            <div class="dec-grid" id="dec-grid">
              <button type="button" class="dec-btn" data-dec="PASS"><i class="s-PASS"></i><span>PASS</span></button>
              <button type="button" class="dec-btn" data-dec="EXCLUDE"><i class="s-EXCLUDE"></i><span>EXCLUDE</span></button>
              <button type="button" class="dec-btn" data-dec="REVIEW"><i class="s-REVIEW"></i><span>REVIEW</span></button>
              <button type="button" class="dec-btn" data-dec="REPEAT"><i class="s-REPEAT"></i><span>REPEAT</span></button>
            </div>
            <div class="dec-save">
              <button type="button" class="ok" id="confirm" title="Save this decision and open the next well">&#10003;</button>
              <button type="button" class="caret" id="caret" title="Reasons and notes">&#9660;</button>
            </div>
          </div>
        </div>

        <div class="drawer" id="drawer" hidden>
          <div class="micro">Reasons</div>
          <div class="rgroups" id="rgroups"></div>
          <textarea class="notes" id="notes" placeholder="Notes (optional) &#8212; saved as you type"></textarea>
        </div>

        <div class="sect"><h3>Diagnostics</h3><span class="rule"></span></div>
        <div class="usable">
          <div class="row">
            <span class="micro">Usable reads</span>
            <span class="v" id="usable-val">&#8211;</span>
          </div>
          <div class="track">
            <i class="fill" id="usable-fill"></i>
            <i class="tri tri-warn" style="left:__WARN_PCT__"></i>
            <i class="tri tri-pass" style="left:__PASS_PCT__"></i>
          </div>
          <div class="ticks">
            <span class="t-warn" style="left:__WARN_PCT__">__WARN_TICK__</span>
            <span class="t-pass" style="left:__PASS_PCT__">__PASS_TICK__</span>
          </div>
        </div>
        <div class="diag" id="diag"></div>

        <div class="sect"><h3>Evidence</h3><span class="rule"></span></div>

        <!-- The whole cell-image layer. qc_review.js hides this row when the payload
             says the plate has no images, leaving Evidence as profile + histogram. -->
        <div class="img-row" id="img-row">
          <div class="img-panel">
            <div class="img-head">
              <span class="micro">CellenONE image</span>
              <div class="seg" id="seg">
                <button type="button" data-ch="merge" class="on">MERGE</button>
                <button type="button" data-ch="trans">TRANSMISSION</button>
              </div>
            </div>
            <!-- No .vline/.vlbl overlay: render_cellenone_images.py already draws the
                 isolation and ejection lines into the image at the run's real x, and a
                 second pair at fixed 38%/66% would contradict them. -->
            <div class="plate" id="plate">
              <img id="plate-img" alt="CellenONE image" hidden>
              <span class="nozzle">nozzle &#8592;</span>
              <span class="cap" id="plate-cap"></span>
            </div>
            <!-- Outline colours follow cellenONE's own particle taxonomy
                 (cellenREPORT S4.3), so a reviewer reading both sees the same words. -->
            <div class="objkey micro" id="objkey">
              <span><i style="background:#3cdc5a"></i>isolated</span>
              <span><i style="background:#ffd228"></i>fits criteria</span>
              <span><i style="background:#ff69b4"></i>detected</span>
            </div>
          </div>
          <div class="img-metrics" id="img-metrics">
            <div class="hd micro">
              <span>Metrics</span>
              <!-- Reveals the isolated-cell shape column. Auto-open on a wide window
                   (see the media query), so the control only matters when space is
                   tight. -->
              <button type="button" class="more" id="img-more"
                      title="Isolated-cell shape (→)" aria-expanded="false">&#9654;</button>
            </div>
            <div class="img-metrics-cols">
              <div id="img-metrics-body"></div>
              <div id="img-metrics-more" class="more-col" hidden></div>
            </div>
          </div>
        </div>

        <figure>
          <div class="frame">
            <div style="display:flex;align-items:baseline;gap:9px;margin-bottom:7px">
              <span class="micro">Copy-number profile</span>
              <span class="mono" style="font-size:12px;color:rgba(29,31,32,.4)" id="profile-note"></span>
            </div>
            <img id="profile-img" alt="copy-number profile" hidden>
            <div class="mono" id="profile-placeholder" style="font-size:12px;color:rgba(29,31,32,.45)" hidden></div>
            <div class="axis" id="profile-axis"></div>
          </div>
        </figure>

        <figure>
          <div class="frame">
            <div class="micro" style="margin-bottom:6px">Bin read counts</div>
            <img id="hist-img" alt="bin read counts" hidden>
            <div class="mono" id="hist-placeholder" style="font-size:12px;color:rgba(29,31,32,.45)" hidden></div>
          </div>
        </figure>

      </div>
    </div>
  </div>
</div>

<script>window.QC = __QC_JSON__;</script>
<script>
__QC_JS__
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
