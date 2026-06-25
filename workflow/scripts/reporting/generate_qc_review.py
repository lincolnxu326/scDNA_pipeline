#!/usr/bin/env python3
"""
Build the static HTML QC review report for one plate.

This reads the per-well QC metrics already produced by the pipeline
(samtools stats, dedup / filter summaries, demux stats) plus the per-well
AneuFinder profile PNGs rendered by render_well_profiles.R, computes an
advisory automated QC status from config.yaml `qc_thresholds`, and writes a
single self-contained `review.html`.

Design constraints (deliberate):
  * Fully static. No server, no database, no external requests.
  * Per-well plot PNGs are base64-embedded so the page opens correctly over
    Samba / as a local file:// document.
  * The browser only ever *generates* copyable text, a downloaded CSV, or a
    copyable shell command; it never writes into the project directory. The
    reproducible human artefact is `<PLATE_DIR>/qc_decisions.csv` (per plate),
    which the user saves by hand or via the in-page "Copy terminal save command".
"""

import argparse
import base64
import json
import logging
import sys
from pathlib import Path

import yaml

# Controlled vocabularies (kept in sync with validate_qc_decisions.py and config/README.md)
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
                   help="multiqc_data directory (optional enrichment)")
    p.add_argument("--plots-dir", required=True, help="Directory of per-well {well}.png")
    p.add_argument("--config", required=True, help="Pipeline config.yaml (for qc_thresholds)")
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
    return p.parse_args()


def load_included_wells(path: str) -> list:
    """Read an included_wells.tsv (sample_id\\twell) and return the well column."""
    if not path:
        return []
    p = Path(path)
    if not p.exists():
        return []
    wells = []
    with open(p) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        try:
            widx = header.index("well")
        except ValueError:
            widx = 1 if len(header) > 1 else 0
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) > widx and parts[widx].strip():
                wells.append(parts[widx].strip())
    return wells


# ---------------------------------------------------------------------------
# Metric parsing
# ---------------------------------------------------------------------------

def parse_samtools_stats(stats_file: Path) -> dict:
    """Parse the SN block of a samtools stats file (post-dedup BAM)."""
    out = {"reads": None, "mapped_reads": None, "average_quality": None,
           "mapping_rate": None}
    if not stats_file.exists():
        return out
    total = mapped = None
    with open(stats_file) as fh:
        for line in fh:
            if line.startswith("SN\traw total sequences:"):
                total = int(line.split("\t")[2])
            elif line.startswith("SN\treads mapped:"):
                mapped = int(line.split("\t")[2])
            elif line.startswith("SN\taverage quality:"):
                out["average_quality"] = float(line.split("\t")[2])
    out["reads"] = total
    out["mapped_reads"] = mapped
    if total and total > 0 and mapped is not None:
        out["mapping_rate"] = round(mapped / total * 100, 2)
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

# auto_status -> default decision mapping (read-count-only gating)
STATUS_DEFAULT_DECISION = {
    "PASS": "PASS",
    "WARN": "REVIEW",
    "FAIL": "EXCLUDE",
    "UNKNOWN": "REVIEW",
}


def compute_status(usable_reads, pass_cutoff: int, warn_cutoff: int) -> tuple:
    """
    Read-count-only gating. Returns (auto_status, default_decision, default_reason, flags).

    usable_reads = mapped reads in the dedup BAM. Duplication is NOT used for gating.
    Advisory only — never overrides the human call.
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
    Map each well to a (row, col) grid position.

    Uses true positional IDs when every well matches A-H + 01-12; otherwise
    falls back to a sequential 8x12 (then wider) fill in the given order.
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
        ncols = 12
        nrows = 8
        return {"mode": "positional", "rows": nrows, "cols": ncols, "pos": pos}

    # Sequential fill: 12 columns, as many rows as needed.
    ncols = 12
    nrows = max(1, -(-len(wells) // ncols))  # ceil
    pos = {w: (i // ncols, i % ncols) for i, w in enumerate(wells)}
    return {"mode": "sequential", "rows": nrows, "cols": ncols, "pos": pos}


# ---------------------------------------------------------------------------
# HTML assembly
# ---------------------------------------------------------------------------

def embed_png(path: Path) -> str:
    if not path.exists():
        return ""
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{data}"


def embed_well_plot(plots_dir: Path, well: str, suffix: str) -> str:
    """Embed {well}_{suffix} as a data URI, preferring SVG (vector, responsive) then PNG.

    Also falls back to the legacy {well}.png (profile only)."""
    svg = plots_dir / f"{well}_{suffix}.svg"
    if svg.exists():
        data = base64.b64encode(svg.read_bytes()).decode("ascii")
        return f"data:image/svg+xml;base64,{data}"
    png = plots_dir / f"{well}_{suffix}.png"
    if png.exists():
        return embed_png(png)
    if suffix == "profile":
        return embed_png(plots_dir / f"{well}.png")
    return ""


def build_html(app: dict) -> str:
    data_blob = json.dumps(app, separators=(",", ":"))
    return HTML_TEMPLATE.replace("__DATA_BLOB__", data_blob)


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
    wells = args.wells          # always the FULL plate (cn viewer greys out excluded wells)
    included_set = set()
    if kind == "cn":
        included = load_included_wells(args.included_wells)
        if args.included_wells and not included:
            logger.warning("No included wells found in %s; CN report will mark all as excluded",
                           args.included_wells)
        included_set = set(included)
        # The cn viewer shows the WHOLE plate: included (PASS) wells keep their original
        # auto-status colour and carry second-pass plots; excluded wells are greyed out.
    logger.info("Building %s report for plate %s (%d wells, %d included)",
                kind, args.plate, len(wells), len(included_set) if kind == "cn" else len(wells))

    # ---- aggregate per-well metrics ----------------------------------------
    dedup = load_tsv_by_well(
        plate_dir / "dedup" / "dedup_summary.tsv",
        {"unique_reads": "unique_reads", "duplicate_reads": "duplicate_reads",
         "dedup_rate": "dedup_rate"},
    )
    dimers = load_tsv_by_well(
        plate_dir / "filtered" / "adapter_filter_summary.tsv",
        {"dimer_rate": "dimer_rate", "kept_pairs": "kept_pairs"},
    )
    demux_totals = load_demux_totals(plate_dir / "demux" / "demux_stats.json")

    geom = well_positions(wells)

    well_data = {}
    n_profiles = 0
    n_hist = 0
    status_counts = {"PASS": 0, "WARN": 0, "FAIL": 0, "UNKNOWN": 0}
    for w in wells:
        sam = parse_samtools_stats(plate_dir / "bam" / f"{w}.stats.txt")
        metrics = dict(sam)
        if w in dedup:
            metrics["unique_reads"] = dedup[w].get("unique_reads")
            metrics["duplicate_reads"] = dedup[w].get("duplicate_reads")
            metrics["dedup_rate"] = dedup[w].get("dedup_rate")
        if w in dimers:
            metrics["dimer_rate"] = dimers[w].get("dimer_rate")
            metrics["kept_pairs"] = dimers[w].get("kept_pairs")
        if w in demux_totals:
            metrics["demux_reads"] = demux_totals[w]

        # UMI dedup retention (diagnostic only): unique / (unique + duplicate)
        uniq = metrics.get("unique_reads")
        dup = metrics.get("duplicate_reads")
        if uniq is not None and dup is not None and (uniq + dup) > 0:
            metrics["umi_retention"] = round(uniq / (uniq + dup) * 100, 2)

        # Read-count-only gating on usable_reads = mapped reads in the dedup BAM.
        usable_reads = metrics.get("mapped_reads")
        metrics["usable_reads"] = usable_reads
        auto_status, default_decision, default_reason, flags = compute_status(
            usable_reads, pass_cutoff, warn_cutoff)
        status_counts[auto_status] = status_counts.get(auto_status, 0) + 1

        # plots (prefer SVG → PNG). In cn-mode only included wells have second-pass plots.
        is_included = (kind != "cn") or (w in included_set)
        profile_uri = embed_well_plot(plots_dir, w, "profile") if is_included else ""
        hist_uri = embed_well_plot(plots_dir, w, "histogram") if is_included else ""
        if profile_uri:
            n_profiles += 1
        if hist_uri:
            n_hist += 1
        r, c = geom["pos"][w]
        well_data[w] = {
            "well": w,
            "sample_id": f"{args.plate}_{w}",
            "row": r,
            "col": c,
            "metrics": metrics,
            "status": auto_status,
            "default_decision": default_decision,
            "default_reason": default_reason,
            "flags": flags,
            "included": is_included,
            "plot": profile_uri,
            "histogram": hist_uri,
        }

    logger.info("Embedded %d/%d profile PNGs, %d/%d histogram PNGs",
                n_profiles, len(wells), n_hist, len(wells))
    logger.info("Auto-status: PASS=%d WARN=%d FAIL=%d UNKNOWN=%d",
                status_counts["PASS"], status_counts["WARN"],
                status_counts["FAIL"], status_counts["UNKNOWN"])

    heatmap_uri = ""
    if kind == "cn" and args.heatmap:
        heatmap_uri = embed_png(Path(args.heatmap))
        if heatmap_uri:
            logger.info("Embedded genome heatmap from %s", args.heatmap)

    if args.title:
        title = args.title
    elif kind == "cn":
        title = "Final CN Review"
    else:
        title = "QC Review"

    decisions_path = args.decisions_path or str(Path(args.plate_dir) / "qc_decisions.csv")
    app = {
        "plate": args.plate,
        "kind": kind,
        "title": title,
        "decisions": DECISIONS,
        "reasons": REASONS,
        "plate_dir": args.plate_dir,
        "decisions_path": decisions_path,
        "decisions_dir": str(Path(decisions_path).parent),
        "cutoffs": {"pass": pass_cutoff, "warn": warn_cutoff},
        "grid": {"rows": geom["rows"], "cols": geom["cols"], "mode": geom["mode"]},
        "wells_order": wells,
        "wells": well_data,
        "heatmap": heatmap_uri,
    }

    html = build_html(app)
    out_name = args.out_name or ("cn_review.html" if kind == "cn" else "review.html")
    out_html = outdir / out_name
    out_html.write_text(html, encoding="utf-8")
    logger.info("Wrote %s", out_html)


# ---------------------------------------------------------------------------
# Self-contained HTML/CSS/JS template. __DATA_BLOB__ is replaced with JSON.
# ---------------------------------------------------------------------------

HTML_TEMPLATE = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>QC Review</title>
<style>
  * { box-sizing: border-box; }
  body { margin:0; font-family: system-ui, -apple-system, Segoe UI, Roboto, Helvetica, Arial, sans-serif;
         color:#1a1a1a; background:#f4f5f7; }
  header { background:#fff; border-bottom:1px solid #ddd; padding:12px 20px;
           display:flex; align-items:center; gap:20px; flex-wrap:wrap; }
  header h1 { font-size:18px; margin:0; }
  header .meta { color:#555; font-size:13px; }
  header a { color:#1565c0; }
  .layout { display:flex; gap:12px; padding:12px; align-items:flex-start; }
  .left { flex:0 0 auto; max-width:470px; }   /* ~ plate-map width; hint/legend wrap, no dead space */
  .right { flex:1 1 auto; min-width:0; }
  .panel { background:#fff; border:1px solid #ddd; border-radius:6px; padding:12px; }
  table.grid { border-collapse:collapse; }
  table.grid th { font-size:11px; color:#666; font-weight:600; padding:2px 4px; text-align:center; }
  table.grid td { padding:2px; }
  .cell { position:relative; width:30px; height:24px; border-radius:3px; border:2px solid transparent;
          color:#fff; font-size:9px; font-weight:600; cursor:pointer; display:flex;
          flex-direction:column; align-items:center; justify-content:center; line-height:1.05;
          user-select:none; }
  .cell:hover { outline:1px solid #333; }
  .cell.sel { border-color:#111; box-shadow:0 0 0 2px #111 inset; }
  .cell.manual::after { content:""; position:absolute; top:2px; right:2px; width:6px; height:6px;
                        border-radius:50%; background:#fff; box-shadow:0 0 0 1px #111; }
  .cell .dec { font-size:7px; letter-spacing:0.2px; opacity:0.95; }
  .legend { display:flex; gap:14px; margin-top:10px; font-size:12px; flex-wrap:wrap; }
  .legend span { display:inline-flex; align-items:center; gap:5px; }
  .swatch { width:13px; height:13px; border-radius:3px; display:inline-block; }
  .dot { width:7px; height:7px; border-radius:50%; background:#fff; box-shadow:0 0 0 1px #111; display:inline-block; }
  .kv { width:100%; border-collapse:collapse; font-size:13px; margin-bottom:12px; }
  .kv td { padding:3px 8px; border-bottom:1px solid #eee; }
  .kv td.k { color:#666; width:46%; }
  .kv tr.gate td { background:#f0f7ff; font-weight:600; }
  .kv .diag { color:#888; font-size:11px; }
  .badge { display:inline-block; padding:2px 8px; border-radius:10px; color:#fff; font-size:12px; font-weight:600; }
  .flags { color:#b71c1c; font-size:12px; }
  .statebadge { display:inline-block; padding:2px 9px; border-radius:10px; font-size:12px; font-weight:600; margin-left:8px; }
  .statebadge.auto { background:#eee; color:#555; }
  .statebadge.saved { background:#1565c0; color:#fff; }
  .ctl { margin:14px 0; }
  .ctl > .lbl { font-size:12px; color:#555; margin-bottom:5px; }
  .btnrow { display:flex; gap:8px; flex-wrap:wrap; }
  .decbtn { font:inherit; font-weight:600; padding:7px 16px; border:2px solid #bbb; background:#fff;
            color:#333; border-radius:5px; cursor:pointer; }
  .decbtn.active { color:#fff; border-color:transparent; }
  .reasonbtn { font:inherit; font-size:12px; padding:5px 11px; border:1px solid #bbb; background:#fff;
               color:#444; border-radius:14px; cursor:pointer; }
  .reasonbtn.active { background:#1565c0; color:#fff; border-color:#1565c0; }
  textarea { font:inherit; padding:6px 8px; border:1px solid #bbb; border-radius:4px; width:100%; min-height:42px; }
  .savebtn { font:inherit; font-weight:600; padding:8px 18px; border:1px solid #2e7d32; background:#2e7d32;
             color:#fff; border-radius:5px; cursor:pointer; margin-top:10px; }
  .plots { display:flex; gap:12px; flex-wrap:wrap; margin-top:12px; }
  .plots.stack { flex-direction:column; flex-wrap:nowrap; }
  .plotbox { flex:1 1 460px; min-width:320px; border:1px solid #eee; border-radius:4px; background:#fafafa; }
  .plots.stack .plotbox { flex:1 1 auto; width:100%; max-width:1200px; }
  .plotbox.narrow { max-width:820px; }
  .plotbox .cap { font-size:12px; color:#666; padding:4px 8px; border-bottom:1px solid #eee; }
  .plotbox img { width:100%; display:block; }
  .cell.excl { opacity:0.5; }
  .noplot { padding:20px; text-align:center; color:#999; font-style:italic; }
  .toolbar { background:#fff; border:1px solid #ddd; border-radius:6px; padding:12px 14px; margin:0 16px 16px; }
  .toolbar button { font:inherit; padding:7px 14px; margin-right:10px; border:1px solid #1565c0;
                    background:#1565c0; color:#fff; border-radius:4px; cursor:pointer; }
  .toolbar button.secondary { background:#fff; color:#1565c0; }
  .toolbar code { background:#eef; padding:1px 5px; border-radius:3px; }
  .hint { font-size:12px; color:#555; margin-top:8px; line-height:1.5; }
  .progress { font-size:12px; color:#444; margin-left:auto; }
  .mapwrap { display:flex; flex-direction:column; gap:12px; align-items:flex-start; }
  .statcol { display:flex; flex-direction:row; gap:10px; flex-wrap:wrap; }
  .statcard { border:1px solid #ddd; border-radius:8px; padding:12px 18px; min-width:150px;
              text-align:center; background:#fafafa; }
  .statcard .lbl { font-size:12px; color:#666; text-transform:uppercase; letter-spacing:.5px; }
  .statcard .pct { font-size:34px; font-weight:800; line-height:1.05; margin-top:6px; }
  .statcard .cnt { font-size:13px; color:#555; margin-top:2px; }
  .statcard.pass .pct { color:#2e7d32; }
  .statcard.review .pct { color:#1565c0; }
</style>
</head>
<body>
<script>const APP = __DATA_BLOB__;</script>
<header>
  <h1><span id="reportTitle">QC Review</span> &mdash; <span id="plateName"></span></h1>
  <span class="meta"><span id="wellCount"></span> wells &middot; grid: <span id="gridMode"></span></span>
  <span class="meta">gating: usable_reads &ge; <span id="cutPass"></span> = PASS, &ge; <span id="cutWarn"></span> = WARN</span>
  <span class="meta"><a id="multiqcLink" href="../multiqc/multiqc_report.html" target="_blank" rel="noopener">Open MultiQC report &#8599;</a></span>
  <span class="progress" id="progress"></span>
</header>

<div class="panel" id="heatmapPanel" style="margin:0 12px 12px; display:none">
  <h2 style="margin:0 0 10px; font-size:16px">Genome-wide copy-number heatmap (second pass)</h2>
  <img id="heatmapImg" style="width:100%; border:1px solid #eee; border-radius:4px" alt="genome heatmap">
</div>

<div class="layout">
  <div class="left panel">
    <div class="mapwrap">
      <div id="plateMap"></div>
      <div class="statcol" id="statcol">
        <div class="statcard pass">
          <div class="lbl">Pass rate</div>
          <div class="pct"><span id="passPct">0</span>%</div>
          <div class="cnt"><span id="passCount">0</span> / <span id="passTotal">0</span> wells</div>
        </div>
        <div class="statcard review">
          <div class="lbl">Review rate</div>
          <div class="pct"><span id="revPct">0</span>%</div>
          <div class="cnt"><span id="revCount">0</span> / <span id="revTotal">0</span> wells</div>
        </div>
      </div>
    </div>
    <div class="legend" id="legend">
      <span><span class="swatch" style="background:#2e7d32"></span>PASS</span>
      <span><span class="swatch" style="background:#c62828"></span>EXCLUDE</span>
      <span><span class="swatch" style="background:#1565c0"></span>REVIEW</span>
      <span><span class="swatch" style="background:#ef6c00"></span>REPEAT</span>
      <span><span class="dot"></span>manually saved</span>
    </div>
    <div class="hint" id="mapHint">Cell colour = current decision (auto-default until you Save). The detail
      panel shows the automated read-count status. A dot marks wells you have saved.</div>
  </div>

  <div class="right panel" id="detail">
    <p style="color:#888">Select a well from the plate map to review it.</p>
  </div>
</div>

<div class="toolbar" id="toolbar">
  <button id="btnCopyCmd">Copy terminal save command</button>
  <button id="btnCopy" class="secondary">Copy CSV only</button>
  <button id="btnDownload" class="secondary">Download qc_decisions.csv</button>
  <span class="hint" id="copyStatus"></span>
  <div class="hint">
    <b>Recommended:</b> click <b>Copy terminal save command</b>, then paste it into a terminal
    on the cluster &mdash; it writes your decisions straight to
    <code id="decPath"></code>.
    <b>⚠ This overwrites any existing file at that path</b> (the command prints a warning if one exists).<br>
    Then run <code>MODE=post_review</code> (validate &rarr; second-pass AneuFinder &rarr; cn_review.html).
    The browser cannot write into the project; the command (or the downloaded file) is how the CSV gets saved.
  </div>
</div>

<script>
(function () {
  var DEC_COLOR = { PASS:"#2e7d32", EXCLUDE:"#c62828", REVIEW:"#1565c0", REPEAT:"#ef6c00" };
  var ST_COLOR  = { PASS:"#2e7d32", WARN:"#f9a825", FAIL:"#c62828", UNKNOWN:"#9e9e9e" };

  var decisions = {};   // well -> {decision, reasons:[], notes}
  var manual = {};      // well -> bool (user has Saved this well)
  var selected = null;

  function initDecisions() {
    APP.wells_order.forEach(function (w) {
      var d = APP.wells[w];
      var r = d.default_reason ? [d.default_reason] : [];
      decisions[w] = { decision: d.default_decision || "REVIEW", reasons: r, notes: "" };
      manual[w] = false;
    });
  }

  function fmt(v, suffix) {
    if (v === null || v === undefined || v === "") return "&ndash;";
    if (typeof v === "number") {
      var s = Number.isInteger(v) ? v.toLocaleString() : v.toFixed(2);
      return s + (suffix || "");
    }
    return v;
  }

  function buildMap() {
    var g = APP.grid;
    var html = '<table class="grid"><thead><tr><th></th>';
    for (var c = 0; c < g.cols; c++) html += '<th>' + (c + 1) + '</th>';
    html += '</tr></thead><tbody>';
    var byPos = {};
    APP.wells_order.forEach(function (w) { var d = APP.wells[w]; byPos[d.row + ":" + d.col] = w; });
    for (var r = 0; r < g.rows; r++) {
      html += '<tr><th>' + String.fromCharCode(65 + r) + '</th>';
      for (var c2 = 0; c2 < g.cols; c2++) {
        var w = byPos[r + ":" + c2];
        if (!w) { html += '<td></td>'; continue; }
        html += '<td><div class="cell" data-well="' + w + '" id="cell-' + w + '">' +
                '<span>' + w + '</span><span class="dec" id="dec-' + w + '"></span></div></td>';
      }
      html += '</tr>';
    }
    html += '</tbody></table>';
    document.getElementById("plateMap").innerHTML = html;
    document.querySelectorAll(".cell").forEach(function (el) {
      el.addEventListener("click", function () { select(el.getAttribute("data-well")); });
    });
    APP.wells_order.forEach(refreshCell);
  }

  function refreshCell(w) {
    var cell = document.getElementById("cell-" + w);
    if (!cell) return;
    if (APP.kind === 'cn') {
      var dd = APP.wells[w];
      // included (PASS) wells keep their original auto-status colour; excluded -> grey
      cell.style.background = dd.included ? (ST_COLOR[dd.status] || "#9e9e9e") : "#d6d6d6";
      if (dd.included) cell.classList.remove("excl"); else cell.classList.add("excl");
      return;
    }
    var dec = decisions[w];
    cell.style.background = DEC_COLOR[dec.decision] || "#9e9e9e";
    var lbl = document.getElementById("dec-" + w);
    if (lbl) lbl.textContent = dec.decision;
    if (manual[w]) cell.classList.add("manual"); else cell.classList.remove("manual");
  }

  function metricRow(k, v, cls) {
    return '<tr' + (cls ? ' class="' + cls + '"' : '') + '><td class="k">' + k + '</td><td>' + v + '</td></tr>';
  }

  function plotBox(cap, uri, well, what, cls) {
    var c = "plotbox" + (cls ? " " + cls : "");
    if (uri) return '<div class="' + c + '"><div class="cap">' + cap + '</div><img src="' + uri + '" alt="' + what + ' ' + well + '"></div>';
    return '<div class="' + c + '"><div class="cap">' + cap + '</div><div class="noplot">No ' + what + ' for ' + well + '</div></div>';
  }

  function select(w) {
    selected = w;
    document.querySelectorAll(".cell").forEach(function (el) { el.classList.remove("sel"); });
    var cell = document.getElementById("cell-" + w);
    if (cell) cell.classList.add("sel");

    var d = APP.wells[w], m = d.metrics || {}, dec = decisions[w];
    var stColor = ST_COLOR[d.status] || "#9e9e9e";

    var html = '';
    html += '<h2 style="margin:0 0 10px">' + d.sample_id;
    if (APP.kind === 'review') {
      html += '<span class="statebadge ' + (manual[w] ? 'saved' : 'auto') + '" id="stateBadge">' +
              (manual[w] ? 'saved (manual)' : 'auto-default') + '</span>';
    }
    html += '</h2>';

    html += '<table class="kv">';
    html += metricRow("Automated status",
      '<span class="badge" style="background:' + stColor + '">' + d.status + '</span>' +
      (d.flags && d.flags.length ? ' <span class="flags">' + d.flags.join(", ") + '</span>' : ''));
    html += metricRow("Usable reads (mapped, dedup BAM)", fmt(m.usable_reads), "gate");
    html += metricRow("Total reads (dedup BAM)", fmt(m.reads));
    html += metricRow("Mapping rate", fmt(m.mapping_rate, "%"));
    html += metricRow("UMI-tools duplicate rate <span class=\"diag\">(diagnostic, not gated)</span>", fmt(m.dedup_rate, "%"));
    html += metricRow("UMI dedup retention <span class=\"diag\">(diagnostic)</span>", fmt(m.umi_retention, "%"));
    html += metricRow("Adapter-dimer rate", fmt(m.dimer_rate, "%"));
    html += metricRow("Unique reads", fmt(m.unique_reads));
    html += metricRow("Demux reads", fmt(m.demux_reads));
    html += metricRow("Avg quality", fmt(m.average_quality));
    html += '</table>';

    if (APP.kind === 'review') {
      // decision buttons
      html += '<div class="ctl"><div class="lbl">Decision</div><div class="btnrow" id="decBtns">';
      APP.decisions.forEach(function (x) {
        var on = x === dec.decision;
        html += '<button type="button" class="decbtn' + (on ? ' active' : '') + '" data-dec="' + x +
                '" style="' + (on ? 'background:' + DEC_COLOR[x] + ';' : '') + '">' + x + '</button>';
      });
      html += '</div></div>';
      // reason multi-select buttons
      html += '<div class="ctl"><div class="lbl">Reasons (select any)</div><div class="btnrow" id="reasonBtns">';
      APP.reasons.forEach(function (x) {
        var on = dec.reasons.indexOf(x) !== -1;
        html += '<button type="button" class="reasonbtn' + (on ? ' active' : '') + '" data-reason="' + x + '">' + x + '</button>';
      });
      html += '</div></div>';
      // notes + save
      html += '<div class="ctl"><div class="lbl">Notes</div><textarea id="txtNotes">' + (dec.notes || '') + '</textarea></div>';
      html += '<button type="button" class="savebtn" id="btnSave">Save decision</button>';
    }

    // plots: review stacks profile/histogram vertically; cn keeps them side-by-side
    html += '<div class="plots' + (APP.kind === 'review' ? ' stack' : '') + '">';
    if (APP.kind === 'cn' && !d.included) {
      html += '<div class="plotbox"><div class="noplot">Excluded at review &mdash; not included in the second AneuFinder pass.</div></div>';
    } else {
      html += plotBox("Copy-number profile", d.plot, d.well, "profile");
      html += plotBox("Bin read-count histogram", d.histogram, d.well, "histogram", "narrow");
    }
    html += '</div>';

    document.getElementById("detail").innerHTML = html;
    if (APP.kind !== 'review') return;

    // wire decision buttons (live recolour; marks manual only on Save)
    document.querySelectorAll("#decBtns .decbtn").forEach(function (b) {
      b.addEventListener("click", function () {
        decisions[w].decision = b.getAttribute("data-dec");
        document.querySelectorAll("#decBtns .decbtn").forEach(function (x) {
          var on = x.getAttribute("data-dec") === decisions[w].decision;
          x.classList.toggle("active", on);
          x.style.background = on ? DEC_COLOR[decisions[w].decision] : "";
        });
        refreshCell(w); updateProgress();
      });
    });
    // wire reason buttons (toggle, multi-select)
    document.querySelectorAll("#reasonBtns .reasonbtn").forEach(function (b) {
      b.addEventListener("click", function () {
        var rs = decisions[w].reasons, val = b.getAttribute("data-reason");
        var i = rs.indexOf(val);
        if (i === -1) { rs.push(val); b.classList.add("active"); }
        else { rs.splice(i, 1); b.classList.remove("active"); }
      });
    });
    document.getElementById("txtNotes").addEventListener("input", function (e) {
      decisions[w].notes = e.target.value;
    });
    document.getElementById("btnSave").addEventListener("click", function () {
      manual[w] = true;
      refreshCell(w); updateProgress();
      var sb = document.getElementById("stateBadge");
      if (sb) { sb.className = "statebadge saved"; sb.textContent = "saved (manual)"; }
      flash("Saved " + w + ": " + decisions[w].decision +
            (decisions[w].reasons.length ? " (" + decisions[w].reasons.join(";") + ")" : ""));
    });
  }

  function csvCell(v) {
    v = (v === null || v === undefined) ? "" : String(v);
    if (/[",\n]/.test(v)) return '"' + v.replace(/"/g, '""') + '"';
    return v;
  }

  function buildCSV() {
    var lines = ["sample_id,well,decision,reason,notes"];
    APP.wells_order.forEach(function (w) {
      var d = decisions[w];
      lines.push([
        csvCell(APP.wells[w].sample_id), csvCell(w),
        csvCell(d.decision), csvCell(d.reasons.join(";")), csvCell(d.notes)
      ].join(","));
    });
    return lines.join("\n") + "\n";
  }

  function updateProgress() {
    var counts = { PASS: 0, EXCLUDE: 0, REVIEW: 0, REPEAT: 0 }, saved = 0;
    APP.wells_order.forEach(function (w) { counts[decisions[w].decision]++; if (manual[w]) saved++; });
    document.getElementById("progress").textContent =
      "PASS " + counts.PASS + " / EXCLUDE " + counts.EXCLUDE +
      " / REVIEW " + counts.REVIEW + " / REPEAT " + counts.REPEAT + "  •  saved " + saved;
    var _tot = APP.wells_order.length;
    var _set = function (id, v) { var e = document.getElementById(id); if (e) e.textContent = v; };
    var _pct = function (n) { return _tot ? (n / _tot * 100).toFixed(1) : "0"; };
    _set("passPct", _pct(counts.PASS)); _set("passCount", counts.PASS); _set("passTotal", _tot);
    _set("revPct", _pct(counts.REVIEW)); _set("revCount", counts.REVIEW); _set("revTotal", _tot);
  }

  function flash(msg) {
    var el = document.getElementById("copyStatus");
    el.textContent = msg;
    setTimeout(function () { el.textContent = ""; }, 2500);
  }

  // header + boot
  document.getElementById("reportTitle").textContent = APP.title || "QC Review";
  document.getElementById("plateName").textContent = APP.plate;
  document.getElementById("wellCount").textContent = APP.wells_order.length;
  document.getElementById("gridMode").textContent = APP.grid.mode;
  document.getElementById("cutPass").textContent = (APP.cutoffs && APP.cutoffs.pass != null) ? APP.cutoffs.pass.toLocaleString() : "?";
  document.getElementById("cutWarn").textContent = (APP.cutoffs && APP.cutoffs.warn != null) ? APP.cutoffs.warn.toLocaleString() : "?";

  if (APP.heatmap) {
    document.getElementById("heatmapImg").src = APP.heatmap;
    document.getElementById("heatmapPanel").style.display = "";
  }

  initDecisions();
  buildMap();

  if (APP.kind === 'review') {
    updateProgress();
  } else {
    // cn viewer: read-only. Hide the decisions toolbar + the pass/review cards, and
    // swap the legend to the status scheme used to colour the kept wells.
    var tb = document.getElementById("toolbar"); if (tb) tb.style.display = "none";
    var sc = document.getElementById("statcol"); if (sc) sc.style.display = "none";
    var lg = document.getElementById("legend");
    if (lg) lg.innerHTML =
      '<span><span class="swatch" style="background:#2e7d32"></span>PASS</span>' +
      '<span><span class="swatch" style="background:#f9a825"></span>WARN</span>' +
      '<span><span class="swatch" style="background:#c62828"></span>FAIL</span>' +
      '<span><span class="swatch" style="background:#9e9e9e"></span>UNKNOWN</span>' +
      '<span><span class="swatch" style="background:#d6d6d6"></span>excluded (not in 2nd pass)</span>';
    var mh = document.getElementById("mapHint");
    if (mh) mh.textContent = "Cell colour = original automated QC status of the wells kept for the " +
      "second pass; excluded wells are greyed out. Click a well to see its final CN profile + histogram.";
  }

  if (APP.wells_order.length) {
    select(APP.wells_order[0]);
  } else {
    document.getElementById("detail").innerHTML = '<p style="color:#888">No wells to display.</p>';
  }

  // Shell command that writes the current decisions straight to the per-plate path.
  // Quoted heredoc delimiter => no shell expansion of CSV/notes content.
  function buildSaveCommand() {
    var path = APP.decisions_path || "qc_decisions.csv";
    var dir = APP.decisions_dir || ".";
    // Overwrite (cat > ...) — writes the full CSV fresh each time. Warn first if the
    // file already exists so an existing decisions file is not clobbered silently.
    return "mkdir -p '" + dir + "'\n" +
           "[ -e '" + path + "' ] && echo 'WARNING: overwriting existing " + path + "'\n" +
           "cat > '" + path + "' <<'QC_DECISIONS_EOF'\n" +
           buildCSV() + "QC_DECISIONS_EOF\n" +
           "echo 'Wrote " + path + "'\n";
  }

  function copyText(text, okMsg) {
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(text).then(
        function () { flash(okMsg); },
        function () { fallbackCopy(text); });
    } else { fallbackCopy(text); }
  }

  if (APP.kind === 'review') {
    var dp = document.getElementById("decPath");
    if (dp) dp.textContent = APP.decisions_path || "<PLATE_DIR>/qc_decisions.csv";

    document.getElementById("btnCopyCmd").addEventListener("click", function () {
      copyText(buildSaveCommand(), "Copied save command — paste it into a cluster terminal.");
    });
    document.getElementById("btnCopy").addEventListener("click", function () {
      copyText(buildCSV(), "Copied " + APP.wells_order.length + " CSV rows to clipboard.");
    });
    document.getElementById("btnDownload").addEventListener("click", function () {
      var blob = new Blob([buildCSV()], { type: "text/csv" });
      var url = URL.createObjectURL(blob);
      var a = document.createElement("a");
      a.href = url; a.download = "qc_decisions.csv";
      document.body.appendChild(a); a.click(); document.body.removeChild(a);
      URL.revokeObjectURL(url);
      flash("Downloaded qc_decisions.csv");
    });
  }

  function fallbackCopy(text) {
    var ta = document.createElement("textarea");
    ta.value = text; document.body.appendChild(ta); ta.select();
    try { document.execCommand("copy"); flash("Copied to clipboard."); }
    catch (e) { flash("Copy failed &mdash; use Download instead."); }
    document.body.removeChild(ta);
  }
})();
</script>
</body>
</html>
"""


if __name__ == "__main__":
    main()
