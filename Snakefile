"""
Single Cell DNA Sequencing Analysis Pipeline
Processes one plate at a time via SLURM submission
"""

import re
import sys

import pandas as pd
from pathlib import Path

# ========================================================================
# Configuration
# ========================================================================

configfile: "config.yaml"

# Get plate directory from command line (passed by submit_pipeline.sh)
if "plate_dir" in config:
    PLATE_DIR = Path(config["plate_dir"]).resolve()
    PLATE = PLATE_DIR.name
    OUTPUT_BASE = PLATE_DIR.parent
  # print(f"Processing plate: {PLATE}")
  # print(f"Plate directory: {PLATE_DIR}")
else:
    raise ValueError("No plate directory specified. Use --config plate_dir=PATH")

# Pipeline directory
PIPELINE_DIR = Path(workflow.basedir).resolve()

# The 384 layout helper is stdlib-only precisely so it can be imported here, at
# parse time, as well as from the reporting scripts (which run in other conda envs).
sys.path.insert(0, str(PIPELINE_DIR / "workflow" / "scripts" / "reporting"))
import plate384_layout as L384

# ========================================================================
# Plate format: 96 (one plate = one subplate = one FASTQ pair) or 384
# ========================================================================
#
# A 384-well plate is dispensed as four interleaved 96-well subplates, each of
# which is sequenced as its own FASTQ pair. 384 mode keeps that per-subplate
# preprocessing exactly as-is (under `<PLATE_DIR>/<subplate>/`) and adds ONE
# plate-level AneuFinder + QC review over all 384 wells.
#
# The whole 96/384 branch reduces to two strings:
#
#   SUB_BASE  the per-subplate output root. In 384 mode it carries a `{sub}`
#             wildcard; in 96 mode it is literally `str(PLATE_DIR)`, so every
#             96-mode path renders CHARACTER-IDENTICAL to before this existed
#             (same on-disk outputs, same .snakemake/metadata keys, no re-runs).
#   LOG_TAG   the per-subplate log filename tag: `{sub}` vs the plate name.
#
# `wid()` is the matching identity function for well ids. Together they let one
# set of rule bodies serve both formats.
PLATE_FORMAT = int(config.get("plate_format", 96))
if PLATE_FORMAT not in (96, 384):
    raise ValueError(f"plate_format must be 96 or 384, got {PLATE_FORMAT}")
IS_384 = (PLATE_FORMAT == 384)

if IS_384:
    SUBPLATES = L384.discover_subplates(PLATE_DIR, expected=config.get("subplates") or None)
    SUB_BASE = str(PLATE_DIR / "{sub}")
    LOG_TAG = "{sub}"
    print(f"Plate format: 384 — {len(SUBPLATES)} subplate(s): {', '.join(SUBPLATES)}")
else:
    SUBPLATES = [PLATE]
    SUB_BASE = str(PLATE_DIR)          # NO wildcard -> zero DAG change
    LOG_TAG = PLATE
    assert len(SUBPLATES) == 1, "96 mode must have exactly one (virtual) subplate"
SUB_LOG = SUB_BASE + "/logs"

# Snakemake wildcards default to `.+`, which matches `/`. Without these the `{sub}`
# wildcard would happily swallow path separators and match the wrong files.
wildcard_constraints:
    sub = "|".join(re.escape(s) for s in SUBPLATES),
    well = r"W\d+"


def wid(sub, well):
    """The single definition of a well's plate-level identity.

    This one string is the staged BAM basename -> AneuFinder model id ->
    `.RData` basename -> plot basename -> HTML well key -> included_wells.tsv row.
    That chain is exactly why `run_aneufinder.R` needs no edits for 384 mode.
    """
    return f"{sub}_{well}" if IS_384 else well

# Per plate directories
LOG_DIR = PLATE_DIR / "logs"
MAPPABILITY_DIR = PLATE_DIR / "mappability"

RESOURCE_DIR = PIPELINE_DIR / "resources"
REFERENCE_RESOURCE_DIR = RESOURCE_DIR / "reference"
MAPPABILITY_RESOURCE_DIR = REFERENCE_RESOURCE_DIR / "mappability"

MAPPABILITY_CONFIG = config.get("mappability", {})
ANEUFINDER_CONFIG = config.get("aneufinder", {})
REFERENCE_ASSEMBLY = ANEUFINDER_CONFIG.get("assembly", "hg38")

VARIABLE_WIDTH_REFERENCE = ANEUFINDER_CONFIG.get("variable_width_reference") or (
    f"{MAPPABILITY_CONFIG['reference_bam']}.bed"
)

SHARED_BLACKLIST = MAPPABILITY_CONFIG.get("blacklist") or str(
    MAPPABILITY_RESOURCE_DIR / "blacklist.bed.gz"
)

GC_RDS = ANEUFINDER_CONFIG.get("gc_rds") or str(
    REFERENCE_RESOURCE_DIR
    / f"{REFERENCE_ASSEMBLY}_binsize{config['aneufinder']['binsize']}_variable_bins_with_GC.rds"
)

ANEUFINDER_CHROMS = ",".join(config["aneufinder"]["chromosomes"])
MAPPABILITY_DIR = MAPPABILITY_RESOURCE_DIR
PLATE_GC_RDS = GC_RDS

# Load barcodes (prefer shared resources, then fall back to the plate directory).
# In 384 mode the fallback must resolve PER SUBPLATE: `<PLATE_DIR>/barcodes.tsv`
# does not exist for a 384 plate — only `<PLATE_DIR>/<subplate>/barcodes.tsv` does.
def _barcode_candidates(d):
    return [
        RESOURCE_DIR / "barcodes.tsv",
        RESOURCE_DIR / "barcodes" / "barcodes.tsv",
        Path(d) / "barcodes.tsv",
    ]


def _find_barcodes(d):
    return next((p for p in _barcode_candidates(d) if p.exists()), None)


if IS_384:
    BARCODES_BY_SUB = {s: _find_barcodes(PLATE_DIR / s) for s in SUBPLATES}
    _missing = [s for s, p in BARCODES_BY_SUB.items() if p is None]
    if _missing:
        raise FileNotFoundError(
            "No barcodes.tsv found for subplate(s): " + ", ".join(_missing) + ".\n"
            + "\n".join(
                f"  {s}: checked " + ", ".join(str(p) for p in _barcode_candidates(PLATE_DIR / s))
                for s in _missing
            )
        )
    # The plate-level well set only makes sense if every subplate uses the same
    # well_id ordering (true today: barcodes.tsv is byte-identical across subplates).
    _orders = {s: pd.read_csv(p, sep="\t", comment="#")["well_id"].tolist()
               for s, p in BARCODES_BY_SUB.items()}
    _ref_sub = SUBPLATES[0]
    for s in SUBPLATES[1:]:
        if _orders[s] != _orders[_ref_sub]:
            raise ValueError(
                f"barcodes.tsv well_id order differs between {_ref_sub} and {s}. "
                "384 mode requires an identical well set/order across subplates "
                f"({BARCODES_BY_SUB[_ref_sub]} vs {BARCODES_BY_SUB[s]})."
            )
    barcodes_file = BARCODES_BY_SUB[_ref_sub]
    WELLS = _orders[_ref_sub]
    print(f"Using {len(WELLS)} well barcodes x {len(SUBPLATES)} subplates "
          f"= {len(WELLS) * len(SUBPLATES)} cells (barcodes from {barcodes_file})")
else:
    BARCODES_BY_SUB = {PLATE: _find_barcodes(PLATE_DIR)}
    barcodes_file = BARCODES_BY_SUB[PLATE]
    if barcodes_file is None:
        searched = ", ".join(str(path) for path in _barcode_candidates(PLATE_DIR))
        raise FileNotFoundError(f"No barcodes.tsv found. Checked: {searched}")
    barcodes_df = pd.read_csv(barcodes_file, sep="\t", comment="#")
    WELLS = barcodes_df["well_id"].tolist()
    print(f"Using {len(WELLS)} well barcodes from {barcodes_file}")


def barcodes_for(wildcards):
    """Per-subplate barcodes.tsv (a params function so `{sub}` resolves at job time)."""
    return str(BARCODES_BY_SUB[wildcards.sub] if IS_384 else barcodes_file)


# QC review sub-stage (static HTML human-review report; runs after multiqc + aneufinder)
REVIEW_CONFIG = config.get("qc_review", {})
REVIEW_DIR = PLATE_DIR / "qc_review"
REVIEW_METHOD = REVIEW_CONFIG.get("method") or config["aneufinder"]["method"][0]
# Human decisions live WITH the plate data (per-plate), not in the shared pipeline dir.
# `decisions_file` is resolved relative to PLATE_DIR unless given as an absolute path.
_decisions_cfg = REVIEW_CONFIG.get("decisions_file") or "qc_decisions.csv"
DECISIONS_FILE = _decisions_cfg if Path(_decisions_cfg).is_absolute() else str(PLATE_DIR / _decisions_cfg)

# Two-pass AneuFinder: the post-review (second) pass reruns only PASS wells into a
# separate directory so the first-pass outputs the review was based on stay intact.
REVIEWED_DIR = PLATE_DIR / (ANEUFINDER_CONFIG.get("reviewed_dir") or "aneufinder_reviewed")
# Final CN deliverable lives in its own clearly-named folder (the last step).
CN_REVIEW_DIR = PLATE_DIR / "CN_review"
REVIEWED_PLOTS_DIR = CN_REVIEW_DIR / "plots"
REVIEWED_HEATMAP = CN_REVIEW_DIR / "genome_heatmap.png"
INCLUDED_WELLS_TSV = REVIEW_DIR / "included_wells.tsv"
INCLUDE_REVIEW_FLAG = "--include-review" if REVIEW_CONFIG.get("include_review") else ""

# Per-well plot resolution. F3: render_well_profiles.R declared --res but never used
# it (PNG_DPI was hardcoded at 150); wiring it is the cheapest lever on 384 asset
# weight — 384 renders at a lower DPI, taking the plots dir from ~81 MB to ~36 MB.
# The 96 default stays 150 so existing reports render pixel-identically.
PLOT_RES = int((REVIEW_CONFIG.get("plot_res_384") or 100) if IS_384
               else (REVIEW_CONFIG.get("plot_res") or 150))
# A 384-well heatmap is far heavier than a 96-well one (F9), so plate-level profile
# rendering gets its own resource block.
QC_REVIEW_RES = (config["resources"].get("qc_review_384", config["resources"]["qc_review"])
                 if IS_384 else config["resources"]["qc_review"])
# `auto` => self-contained single file at 96, sidecar assets at 384 (a 384 report with
# everything base64-embedded would be ~107 MB and unopenable).
_embed_cfg = REVIEW_CONFIG.get("embed_assets", "auto")
ASSETS_MODE = ("sidecar" if IS_384 else "embed") if _embed_cfg == "auto" else (
    "embed" if _embed_cfg else "sidecar")

# ------------------------------------------------------------------------
# CellenONE cell images (optional; display-only QC layer)
# ------------------------------------------------------------------------
CELLENONE_CONFIG = config.get("cellenone", {}) or {}
CELLENONE_DIR = PLATE_DIR / "cellenone"
# The CellenONE run folder -> plate mapping is NOT derivable from any name; it must
# be given explicitly in config.yaml. Missing entry => the feature is simply off.
CELLENONE_RUN = (CELLENONE_CONFIG.get("runs") or {}).get(PLATE, "")
CELLENONE_ENABLED = bool(CELLENONE_CONFIG.get("enable", True)) and bool(CELLENONE_RUN)
CELL_IMAGES_ENABLED = CELLENONE_ENABLED and bool(REVIEW_CONFIG.get("cell_images", True))

# ========================================================================
# Helper Functions
# ========================================================================

def get_input_fastqs(wildcards):
    """The FASTQ pair for this (sub)plate.

    The naming convention requires the FASTQ basename to equal its directory name,
    so 96 mode reads `<PLATE_DIR>/<plate>_R{1,2}.fastq.gz` and 384 mode reads
    `<PLATE_DIR>/<sub>/<sub>_R{1,2}.fastq.gz`.
    """
    if IS_384:
        sub = wildcards.sub
        base = PLATE_DIR / sub / sub
    else:
        base = PLATE_DIR / PLATE
    return {"r1": f"{base}_R1.fastq.gz", "r2": f"{base}_R2.fastq.gz"}

# ========================================================================
# Main Rules
# ========================================================================

rule all:
    """Complete pipeline."""
    input:
        # blacklist generated for this plate
        str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        # plate specific GC template rds (must be created in relaxed env beforehand)
        # aneufinder done (ONE plate-level pass, over all subplates in 384 mode)
        str(PLATE_DIR / "aneufinder" / "complete.flag"),
        # multiqc done (per subplate)
        expand(SUB_BASE + "/multiqc/multiqc_report.html", sub=SUBPLATES),
        # static human QC review report (after multiqc + aneufinder)
        str(REVIEW_DIR / "review.html")

rule all_preprocessing:
    """Run up to alignment."""
    input:
        expand(SUB_BASE + "/bam/{well}.bam", sub=SUBPLATES, well=WELLS)

rule generate_blacklist_plots:
    """Generate blacklist plots only."""
    input:
        str(MAPPABILITY_DIR / "mappability_plot.pdf"),
        str(MAPPABILITY_DIR / "blacklist_diagnosis.txt")

rule all_qc:
    """Generate QC reports only."""
    input:
        expand(SUB_BASE + "/multiqc/multiqc_report.html", sub=SUBPLATES)

rule cellenone_report:
    """CellenONE cell-image ingest + rendering only (iterate on images cheaply)."""
    input:
        str(CELLENONE_DIR / "complete.flag")

rule qc_review_report:
    """Build the static human QC review report (after multiqc + aneufinder)."""
    input:
        str(REVIEW_DIR / "review.html")

rule all_aneufinder_first:
    """First-pass AneuFinder (alias of the existing run; symmetric MODE naming)."""
    input:
        str(PLATE_DIR / "aneufinder" / "complete.flag")

rule all_aneufinder_reviewed:
    """Second (post-review) AneuFinder pass on PASS wells only."""
    input:
        str(REVIEWED_DIR / "complete.flag")

rule cn_review_report:
    """Final copy-number review viewer over the second AneuFinder pass."""
    input:
        str(CN_REVIEW_DIR / "cn_review.html")

rule post_review:
    """All post-review steps in one target (the normal MODE=post_review entry point).

    Reuses the existing chain: validate_qc_decisions -> derive_included_wells ->
    run_aneufinder_reviewed -> render_reviewed_profiles -> cn_review. Requires the
    human-saved <PLATE_DIR>/qc_decisions.csv; no pre-review target depends on it.
    """
    input:
        str(CN_REVIEW_DIR / "cn_review.html")

# ------------------------------------------------------------------------
# Preprocessing Rules
# ------------------------------------------------------------------------

rule demultiplex:
    input:
        unpack(get_input_fastqs)
    output:
        fastqs = expand(SUB_BASE + "/demux/{well}_R{read}.fastq.gz",
                        well=WELLS, read=[1, 2], allow_missing=True),
        stats = SUB_BASE + "/demux/demux_stats.json"
    params:
        outdir = SUB_BASE + "/demux",
        logdir = SUB_LOG,
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "preprocessing" / "extract_umi_barcode.py"),
        barcodes = barcodes_for,
        batch_size = config["batch_size"],
        umi_length = config["preprocessing"]["umi_length"],
        bc1_length = config["preprocessing"]["barcode1_length"],
        bc2_length = config["preprocessing"]["barcode2_length"],
        bc1_offset = config["preprocessing"]["barcode1_offset"],
        bc2_offset = config["preprocessing"]["barcode2_offset"],
        r1_trim = config["preprocessing"]["r1_trim_5prime"],
        r2_trim = config["preprocessing"]["r2_trim_5prime"]
    log:
        SUB_LOG + f"/demux_{LOG_TAG}.log"
    threads: config["resources"]["demux"]["threads"]
    resources:
        mem_mb = config["resources"]["demux"]["mem_mb"],
        runtime = config["resources"]["demux"]["time"],
        partition = config["resources"]["demux"]["partition"]
    conda:
        "workflow/envs/preprocessing.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        python {params.script} \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --barcodes {params.barcodes} \
            --outdir {params.outdir} \
            --batch-size {params.batch_size} \
            --umi-length {params.umi_length} \
            --bc1-length {params.bc1_length} \
            --bc2-length {params.bc2_length} \
            --bc1-offset {params.bc1_offset} \
            --bc2-offset {params.bc2_offset} \
            --r1-trim {params.r1_trim} \
            --r2-trim {params.r2_trim} \
            --stats {output.stats} \
            --log {log} 2>&1
        """

rule filter_dimers:
    input:
        SUB_BASE + "/demux/demux_stats.json"
    output:
        summary = SUB_BASE + "/filtered/adapter_filter_summary.tsv",
        stats = SUB_BASE + "/filtered/filter_stats.json",
        fastqs = expand(SUB_BASE + "/filtered/{well}_R{read}.filtered.fastq.gz",
                        well=WELLS, read=[1, 2], allow_missing=True)
    params:
        indir = SUB_BASE + "/demux",
        outdir = SUB_BASE + "/filtered",
        logdir = SUB_LOG,
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "preprocessing" / "filter_adapter_dimers.py"),
        adapter = config["preprocessing"]["adapter_sequence"],
        case_insensitive = "--case-insensitive" if config["preprocessing"]["case_insensitive"] else "",
        both_reads = "--both-reads" if config["preprocessing"]["filter_both_reads"] else ""
    log:
        SUB_LOG + f"/filter_{LOG_TAG}.log"
    threads: config["resources"]["filter"]["threads"]
    resources:
        mem_mb = config["resources"]["filter"]["mem_mb"],
        runtime = config["resources"]["filter"]["time"],
        partition = config["resources"]["filter"]["partition"]
    conda:
        "workflow/envs/preprocessing.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        python {params.script} \
            --indir {params.indir} \
            --outdir {params.outdir} \
            --adapter {params.adapter} \
            {params.case_insensitive} \
            {params.both_reads} \
            --suffix .fastq.gz \
            --out-suffix .filtered.fastq.gz \
            --stats {output.stats} \
            --log {log} 2>&1
        """

# ------------------------------------------------------------------------
# QC Rules
# ------------------------------------------------------------------------

rule fastqc:
    input:
        r1 = SUB_BASE + "/filtered/{well}_R1.filtered.fastq.gz",
        r2 = SUB_BASE + "/filtered/{well}_R2.filtered.fastq.gz"
    output:
        html1 = SUB_BASE + "/fastqc/{well}_R1_fastqc.html",
        html2 = SUB_BASE + "/fastqc/{well}_R2_fastqc.html",
        zip1 = SUB_BASE + "/fastqc/{well}_R1_fastqc.zip",
        zip2 = SUB_BASE + "/fastqc/{well}_R2_fastqc.zip"
    params:
        outdir = SUB_BASE + "/fastqc",
        logdir = SUB_LOG
    log:
        SUB_LOG + "/fastqc/" + LOG_TAG + "_{well}.log"
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 30,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}/fastqc
        fastqc -o {params.outdir} -t {threads} {input.r1} {input.r2} 2>&1 | tee {log}
        for f in {params.outdir}/*_R[12].filtered_fastqc.*; do
            if [ -f "$f" ]; then
                newname=$(echo "$f" | sed 's/.filtered_fastqc/_fastqc/')
                mv "$f" "$newname" 2>/dev/null || true
            fi
        done
        """

# ------------------------------------------------------------------------
# Alignment Rules
# ------------------------------------------------------------------------

rule build_index:
    input:
        fasta = config["genome"]["fasta"]
    output:
        expand(config["genome"]["index_prefix"] + ".{ext}.bt2",
               ext=["1", "2", "3", "4", "rev.1", "rev.2"])
    params:
        prefix = config["genome"]["index_prefix"]
    log:
        str(LOG_DIR / "bowtie2_index.log")
    threads: 8
    resources:
        mem_mb = 32000,
        runtime = 120,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        """
        mkdir -p {LOG_DIR}
        bowtie2-build --threads {threads} {input.fasta} {params.prefix} 2>&1 | tee {log}
        """

rule align:
    input:
        r1 = SUB_BASE + "/filtered/{well}_R1.filtered.fastq.gz",
        r2 = SUB_BASE + "/filtered/{well}_R2.filtered.fastq.gz"
    output:
        bam = SUB_BASE + "/raw_bam/{well}.bam",
        bai = SUB_BASE + "/raw_bam/{well}.bam.bai",
        stats = SUB_BASE + "/raw_bam/{well}.stats.txt",
        flagstat = SUB_BASE + "/raw_bam/{well}.flagstat.txt"
    params:
        index_prefix = config["genome"]["index_prefix"],
        bowtie2_params = config["alignment"]["bowtie2_params"],
        workdir = SUB_BASE,
        logdir = SUB_LOG,
        # RG-ID gains the subplate in 384 mode (it already carried the plate name).
        rgid_ = lambda wildcards: f"{getattr(wildcards, 'sub', PLATE)}_{wildcards.well}",
        # DO NOT touch rgsm_: changing SM would invalidate every existing BAM header
        # and force a realignment of all 768 BAMs.
        rgsm_ = lambda wildcards: wildcards.well
    log:
        SUB_LOG + "/align/" + LOG_TAG + "_{well}.log"
    threads: config["resources"]["alignment"]["threads"]
    resources:
        mem_mb = config["resources"]["alignment"]["mem_mb"],
        runtime = config["resources"]["alignment"]["time"],
        partition = config["resources"]["alignment"]["partition"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        """
        mkdir -p {params.workdir}/raw_bam {params.logdir}/align
        (bowtie2 -x {params.index_prefix} \
            -1 {input.r1} -2 {input.r2} \
            --threads {threads} {params.bowtie2_params} \
            --rg-id {params.rgid_} \
            --rg "SM:{params.rgsm_}" \
            --rg "PL:ILLUMINA" \
            2> {log} \
        | samtools view -bS - \
        | samtools sort -@ {threads} -o {output.bam} -)
        samtools index {output.bam}
        samtools stats {output.bam} > {output.stats}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule deduplicate_bam:
    input:
        bam = SUB_BASE + "/raw_bam/{well}.bam",
        bai = SUB_BASE + "/raw_bam/{well}.bam.bai"
    output:
        bam = SUB_BASE + "/bam/{well}.bam",
        bai = SUB_BASE + "/bam/{well}.bam.bai",
        stats = SUB_BASE + "/bam/{well}.stats.txt",
        flagstat = SUB_BASE + "/bam/{well}.flagstat.txt"
    params:
        method = config["preprocessing"]["umi_tools_method"],
        workdir = SUB_BASE,
        logdir = SUB_LOG,
        tmpdir = SUB_BASE + "/dedup/tmp/{well}"
    log:
        SUB_LOG + "/dedup/" + LOG_TAG + "_{well}.log"
    threads: config["resources"]["dedup"]["threads"]
    resources:
        mem_mb = config["resources"]["dedup"]["mem_mb"],
        runtime = config["resources"]["dedup"]["time"],
        partition = config["resources"]["dedup"]["partition"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        """
        mkdir -p {params.workdir}/bam {params.tmpdir} {params.logdir}/dedup
        umi_tools dedup \
            --stdin {input.bam} \
            --stdout {output.bam} \
            --paired \
            --extract-umi-method=read_id \
            --umi-separator "_" \
            --method {params.method} \
            --temp-dir {params.tmpdir} \
            --log={log}
        samtools index {output.bam}
        samtools stats {output.bam} > {output.stats}
        samtools flagstat {output.bam} > {output.flagstat}
        """

rule summarize_dedup:
    input:
        raw_bams = expand(SUB_BASE + "/raw_bam/{well}.bam", well=WELLS, allow_missing=True),
        dedup_bams = expand(SUB_BASE + "/bam/{well}.bam", well=WELLS, allow_missing=True)
    output:
        summary = SUB_BASE + "/dedup/dedup_summary.tsv",
        stats = SUB_BASE + "/dedup/dedup_stats.json"
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "preprocessing" / "summarize_umi_tools_dedup.py"),
        raw_dir = SUB_BASE + "/raw_bam",
        dedup_dir = SUB_BASE + "/bam",
        outdir = SUB_BASE + "/dedup",
        logdir = SUB_LOG,
        method = config["preprocessing"]["umi_tools_method"],
        wells = " ".join(WELLS)
    log:
        SUB_LOG + f"/dedup_summary_{LOG_TAG}.log"
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        python {params.script} \
            --raw-dir {params.raw_dir} \
            --dedup-dir {params.dedup_dir} \
            --outdir {params.outdir} \
            --method {params.method} \
            --stats {output.stats} \
            --log {log} \
            --wells {params.wells}
        """


# ------------------------------------------------------------------------
# Blacklist Generation
# ------------------------------------------------------------------------

rule get_mappability_bam:
    output:
        bam = str(MAPPABILITY_DIR / "mappability.bam"),
        bai = str(MAPPABILITY_DIR / "mappability.bam.bai")
    params:
        reference_bam = config.get("mappability", {}).get("reference_bam", "")
    log:
        str(LOG_DIR / "mappability.log")
    threads: 1
    resources:
        mem_mb = 1000,
        runtime = 10,
        partition = config["resources"]["default"]["partition"]
    shell:
        r"""
        mkdir -p {MAPPABILITY_DIR} {LOG_DIR}

        REF="{params.reference_bam}"

        if [ -z "$REF" ]; then
            echo "ERROR: Mappability reference BAM is required!" | tee {log}
            exit 1
        fi

        if [ ! -f "$REF" ]; then
            echo "ERROR: Reference BAM not found: $REF" | tee {log}
            exit 1
        fi

        echo "Using existing mappability BAM: $REF" | tee {log}
        cp "$REF" {output.bam}

        if [ -f "${{REF}}.bai" ]; then
            cp "${{REF}}.bai" {output.bai}
        else
            REF_NOBAM="${{REF%.bam}}"
            if [ -f "${{REF_NOBAM}}.bai" ]; then
                cp "${{REF_NOBAM}}.bai" {output.bai}
            else
                echo "Creating BAM index..." | tee -a {log}
                samtools index {output.bam}
            fi
        fi
        """

rule generate_blacklist:
    input:
        bam = str(MAPPABILITY_DIR / "mappability.bam")
    output:
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        plot = str(MAPPABILITY_DIR / "mappability_plot.pdf"),
        stats = str(MAPPABILITY_DIR / "blacklist_diagnosis.txt")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "generate_blacklist.R"),
        binsize = config["blacklist"]["binsize"],
        lower = config["blacklist"]["lower_quantile"],
        upper = config["blacklist"]["upper_quantile"],
        chromosomes = ",".join(config["aneufinder"]["chromosomes"])
    log:
        str(LOG_DIR / "blacklist.log")
    threads: 4
    resources:
        mem_mb = 16000,
        runtime = 120,
        partition = config["resources"]["default"]["partition"]
    shell:
        """
        mkdir -p {MAPPABILITY_DIR} {LOG_DIR}
        Rscript {params.script} \
            --input {input.bam} \
            --binsize {params.binsize} \
            --lower-quantile {params.lower} \
            --upper-quantile {params.upper} \
            --chromosomes {params.chromosomes} \
            --output-bed {output.blacklist} \
            --output-plot {output.plot} \
            --output-stats {output.stats} \
            2>&1 | tee {log}
        """

# ------------------------------------------------------------------------
# AneuFinder
# ------------------------------------------------------------------------

rule all_aneufinder:
    input:
        str(PLATE_DIR / "aneufinder" / "complete.flag")

rule check_gc_rds:
    input:
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz")
    output:
        flag = str(PLATE_DIR / "logs" / "gc_rds_ready.flag")
    params:
        gc_ready = str(PLATE_DIR / "logs" / "gc_rds_ready.flag")
    log:
        str(PLATE_DIR / "logs" / "check_gc_rds.log")
    shell:
        r"""
        mkdir -p $(dirname {output.flag})
        if [ ! -f "{params.gc_rds}" ]; then
            echo "Missing GC RDS: {params.gc_rds}" > {log}
            echo "Generate it manually with ad_hoc_checks/generate_gc_corr.R after blacklist is created." >> {log}
            exit 1
        fi
        echo "Found GC RDS: {params.gc_rds}" > {log}
        touch {output.flag}
        """

# ------------------------------------------------------------------------
# 384 model reuse: the fast path
# ------------------------------------------------------------------------
#
# AneuFinder is strictly per-cell — each BAM is independently binned, GC-corrected
# and segmented — so a model computed in a 96-well subplate run is IDENTICAL to the
# one a plate-level run would produce. Only the cluster PDF / heatmap are batch-wide,
# and render_well_profiles.R regenerates those at plate level from the .RData anyway.
#
# So in `reuse` mode we symlink the four subplates' existing models into the
# plate-level MODELS dir under their namespaced ids and skip AneuFinder entirely.
# Because render_well_profiles.R derives the plot name from the .RData basename,
# that alone yields correctly-named plate-level plots at zero compute cost.
ANEUFINDER_METHODS = list(config["aneufinder"]["method"])
ANEUFINDER_384_MODE = str(ANEUFINDER_CONFIG.get("plate384_models", "auto")).lower()
if ANEUFINDER_384_MODE not in ("auto", "reuse", "rerun"):
    raise ValueError("aneufinder.plate384_models must be auto|reuse|rerun, "
                     f"got {ANEUFINDER_384_MODE!r}")


def _subplate_models_present():
    """True when every subplate already has models for every configured method."""
    for s in SUBPLATES:
        for m in ANEUFINDER_METHODS:
            d = PLATE_DIR / s / "aneufinder" / "MODELS" / f"method-{m}"
            if not d.is_dir() or not any(d.glob("*.RData")):
                return False
    return True


# Decided at PARSE time so `threads`/`resources` match what the job will actually do.
if IS_384 and ANEUFINDER_384_MODE == "reuse":
    REUSE_384_MODELS = True
elif IS_384 and ANEUFINDER_384_MODE == "auto":
    REUSE_384_MODELS = _subplate_models_present()
else:
    REUSE_384_MODELS = False

if REUSE_384_MODELS:
    print("AneuFinder: reusing existing per-subplate models (no AneuFinder compute)")
    ANEUFINDER_RES = config["resources"].get(
        "aneufinder_link",
        {"threads": 1, "mem_mb": 4000, "time": 30,
         "partition": config["resources"]["default"]["partition"]})
elif IS_384:
    ANEUFINDER_RES = config["resources"].get("aneufinder_384", config["resources"]["aneufinder"])
else:
    ANEUFINDER_RES = config["resources"]["aneufinder"]


def _bam_dir(sub):
    return (PLATE_DIR / sub / "bam") if IS_384 else (PLATE_DIR / "bam")


# `--bam-dir SUB=PATH` once per subplate; in 96 mode that is a single entry.
BAM_DIR_ARGS = " ".join(f"--bam-dir {s}={_bam_dir(s)}" for s in SUBPLATES)
# The stager's basename template, derived from wid() itself so there is exactly ONE
# definition of well identity: wid("{sub}", "{well}") -> "{sub}_{well}" @384, "{well}" @96.
# That basename becomes the AneuFinder model id, the .RData name, the plot name and
# the HTML well key.
ID_TEMPLATE = wid("{sub}", "{well}")


rule run_aneufinder:
    """Run AneuFinder for copy number analysis (requires plate specific GC template).

    ONE plate-level pass. In 384 mode it covers all subplates at once, with wells
    namespaced as `<subplate>_<well>`; with `plate384_models: reuse` it instead
    symlinks the subplates' already-computed models and runs no AneuFinder at all.
    """
    input:
        bams = expand(SUB_BASE + "/bam/{well}.bam", sub=SUBPLATES, well=WELLS),
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        gc_rds = str(PLATE_GC_RDS)
    output:
        flag = str(PLATE_DIR / "aneufinder" / "complete.flag")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "run_aneufinder.R"),
        stager = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "stage_aneufinder_input.py"),
        linker = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "link_subplate_models.py"),
        bam_dir_args = BAM_DIR_ARGS,
        # lambda-wrapped: a params STRING containing `{well}` would otherwise be
        # treated as a wildcard reference and fail to expand in a wildcard-free rule.
        id_template = lambda wildcards: ID_TEMPLATE,
        subplate_dirs = " ".join(str(PLATE_DIR / s / "aneufinder") for s in SUBPLATES),
        subplate_names = " ".join(SUBPLATES),
        reuse_models = "TRUE" if REUSE_384_MODELS else "FALSE",
        outdir = str(PLATE_DIR / "aneufinder"),
        logdir = str(LOG_DIR),
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        gc_rds = str(PLATE_GC_RDS),
        method = ",".join(ANEUFINDER_METHODS),
        binsize = config["aneufinder"]["binsize"],
        chromosomes = ",".join(config["aneufinder"]["chromosomes"]),
        num_cpu = config["aneufinder"]["num_cpu"],
        cluster_plots = "TRUE" if config["aneufinder"]["cluster_plots"] else "FALSE",
        reuse = "TRUE" if config["aneufinder"]["reuse_existing"] else "FALSE",
        min_reads = config["aneufinder"].get("min_reads_for_model", 100)
    log:
        str(LOG_DIR / f"aneufinder_{PLATE}.log")
    threads: ANEUFINDER_RES["threads"]
    resources:
        mem_mb = ANEUFINDER_RES["mem_mb"],
        runtime = ANEUFINDER_RES["time"],
        partition = ANEUFINDER_RES["partition"]
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        if [ "{params.reuse_models}" = "TRUE" ]; then
            # Fast path: the per-subplate 96-well runs already produced every model.
            python {params.linker} \
                --plate-dir {PLATE_DIR} \
                --subplates {params.subplate_names} \
                --methods {params.method} \
                --outdir {params.outdir} \
                2>&1 | tee {log}
        else
            # AneuFinder aborts the ENTIRE batch if any single well has ~no reads
            # (GC-loess yields "NAs found in reads"). The stager keeps only wells above
            # --min-reads; dead wells sit ~1000x below the read-count PASS gate
            # (auto-FAIL) and carry no usable CN signal. run_aneufinder.R is unchanged -
            # it just globs this staged dir, and the staged basenames become the model ids.
            STAGE={params.outdir}/input_bams
            python {params.stager} \
                {params.bam_dir_args} \
                --id-template '{params.id_template}' \
                --min-reads {params.min_reads} \
                --out-dir "$STAGE" \
                --report {params.outdir}/input_bams.tsv \
                2>&1 | tee {log}
            Rscript {params.script} \
                --input "$STAGE" \
                --output {params.outdir} \
                --blacklist {params.blacklist} \
                --gc-rds {params.gc_rds} \
                --method {params.method} \
                --binsize {params.binsize} \
                --chromosomes {params.chromosomes} \
                --numcpu {params.num_cpu} \
                --cluster-plots {params.cluster_plots} \
                --reuse-existing {params.reuse} \
                2>&1 | tee -a {log}
        fi
        touch {output.flag}
        """

# ------------------------------------------------------------------------
# MultiQC
# ------------------------------------------------------------------------

rule multiqc:
    input:
        expand(SUB_BASE + "/fastqc/{well}_R{read}_fastqc.html",
               well=WELLS, read=[1, 2], allow_missing=True),
        SUB_BASE + "/demux/demux_stats.json",
        SUB_BASE + "/dedup/dedup_summary.tsv",
        SUB_BASE + "/filtered/adapter_filter_summary.tsv",
        expand(SUB_BASE + "/bam/{well}.stats.txt", well=WELLS, allow_missing=True),
        expand(SUB_BASE + "/bam/{well}.flagstat.txt", well=WELLS, allow_missing=True)
    output:
        report = SUB_BASE + "/multiqc/multiqc_report.html",
        data = directory(SUB_BASE + "/multiqc/multiqc_data")
    params:
        outdir = SUB_BASE + "/multiqc",
        # In 384 mode MultiQC must scan only this subplate, not the whole plate.
        workdir = SUB_BASE,
        logdir = SUB_LOG,
        config_file = str(PIPELINE_DIR / "multiqc_config.yaml"),
        extra = config.get("multiqc", {}).get("extra_params", "")
    log:
        SUB_LOG + f"/multiqc_{LOG_TAG}.log"
    threads: config["resources"]["multiqc"]["threads"]
    resources:
        mem_mb = config["resources"]["multiqc"]["mem_mb"],
        runtime = config["resources"]["multiqc"]["time"],
        partition = config["resources"]["multiqc"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        multiqc \
            --force \
            --outdir {params.outdir} \
            --config {params.config_file} \
            {params.extra} \
            {params.workdir} \
            2>&1 | tee {log}
        """

# ------------------------------------------------------------------------
# CellenONE cell images (optional, display-only)
# ------------------------------------------------------------------------
#
# The CellenONE dispenser photographs every printed drop in three channels
# (Transmission / Blue / Orange) before ejecting it into the well. Those images say
# whether one cell, several, or none actually made it into the well — information the
# read counts cannot provide.
#
# DISPLAY ONLY. The image call is its own labelled panel and, at most, a *suggested*
# reason chip the human may click. Read-count gating remains the sole driver of
# auto_status and the default decision (see compute_status() in generate_qc_review.py).

rule ingest_cellenone:
    """Parse the CellenONE run folder into per-well records (parse-only, no images)."""
    output:
        wells = str(CELLENONE_DIR / "wells_raw.tsv"),
        meta = str(CELLENONE_DIR / "run_meta.json")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "ingest_cellenone.py"),
        run_dir = CELLENONE_RUN,
        outdir = str(CELLENONE_DIR),
        logdir = str(LOG_DIR),
        plate = PLATE,
        config_file = str(PIPELINE_DIR / "config.yaml"),
        subplates_arg = ("--subplates " + " ".join(SUBPLATES)) if IS_384 else "",
        layout_arg = (f"--layout-tsv {config['plate_layout_tsv']}"
                      if IS_384 and config.get("plate_layout_tsv") else ""),
        wells = " ".join(WELLS)
    log:
        str(LOG_DIR / f"cellenone_ingest_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/cellenone.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        python {params.script} \
            --run-dir '{params.run_dir}' \
            --plate {params.plate} \
            --config {params.config_file} \
            --outdir {params.outdir} \
            {params.subplates_arg} {params.layout_arg} \
            --wells {params.wells} \
            2>&1 | tee {log}
        """

rule render_cellenone_images:
    """Detect objects, apply the ejection-line QC rule, and render the composite JPEGs."""
    input:
        wells = str(CELLENONE_DIR / "wells_raw.tsv"),
        meta = str(CELLENONE_DIR / "run_meta.json")
    output:
        flag = str(CELLENONE_DIR / "complete.flag"),
        wells = str(CELLENONE_DIR / "cellenone_wells.tsv"),
        objects = str(CELLENONE_DIR / "objects.tsv")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "render_cellenone_images.py"),
        outdir = str(CELLENONE_DIR),
        logdir = str(LOG_DIR),
        config_file = str(PIPELINE_DIR / "config.yaml")
    log:
        str(LOG_DIR / f"cellenone_images_{PLATE}.log")
    threads: config["resources"].get("cellenone", {}).get("threads", 8)
    resources:
        mem_mb = config["resources"].get("cellenone", {}).get("mem_mb", 16000),
        runtime = config["resources"].get("cellenone", {}).get("time", 120),
        partition = config["resources"].get("cellenone", {}).get(
            "partition", config["resources"]["default"]["partition"])
    conda:
        "workflow/envs/cellenone.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        python {params.script} \
            --cellenone-dir {params.outdir} \
            --config {params.config_file} \
            --threads {threads} \
            2>&1 | tee {log}
        touch {output.flag}
        """

# ------------------------------------------------------------------------
# QC Review (static HTML, after MultiQC + AneuFinder)
# ------------------------------------------------------------------------

rule render_well_profiles:
    """Render one AneuFinder copy-number profile PNG per well from existing .RData models.

    Read-only consumer of run_aneufinder.R output (does not modify it). Runs via
    renv (no conda:), exactly like run_aneufinder / generate_blacklist.
    """
    input:
        flag = str(PLATE_DIR / "aneufinder" / "complete.flag")
    output:
        manifest = str(REVIEW_DIR / "plots" / "manifest.json")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "render_well_profiles.R"),
        models_dir = str(PLATE_DIR / "aneufinder" / "MODELS"),
        outdir = str(REVIEW_DIR / "plots"),
        logdir = str(LOG_DIR),
        res = PLOT_RES,
        method = REVIEW_METHOD
    log:
        str(LOG_DIR / f"qc_review_plots_{PLATE}.log")
    threads: 4
    resources:
        mem_mb = QC_REVIEW_RES["mem_mb"],
        runtime = QC_REVIEW_RES["time"],
        partition = config["resources"]["default"]["partition"]
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        Rscript {params.script} \
            --input {params.models_dir} \
            --outdir {params.outdir} \
            --method {params.method} \
            --format png \
            --res {params.res} \
            2>&1 | tee {log}
        """

def qc_review_inputs(wildcards=None):
    """Inputs for the plate-level review: per-subplate MultiQC + BAM stats, plus plots."""
    ins = {
        "multiqc": expand(SUB_BASE + "/multiqc/multiqc_report.html", sub=SUBPLATES),
        "multiqc_data": expand(SUB_BASE + "/multiqc/multiqc_data", sub=SUBPLATES),
        "manifest": str(REVIEW_DIR / "plots" / "manifest.json"),
        "bam_stats": expand(SUB_BASE + "/bam/{well}.stats.txt", sub=SUBPLATES, well=WELLS),
    }
    if CELL_IMAGES_ENABLED:
        ins["cellenone"] = str(CELLENONE_DIR / "complete.flag")
    return ins


rule qc_review:
    """Build the static HTML QC review report (metadata + per-well plots).

    96 mode embeds every asset as a data URI (one self-contained file). 384 mode
    writes assets alongside as `qc_review/assets/` and references them relatively,
    so the page stays small and the browser fetches only the clicked well.
    """
    input:
        unpack(qc_review_inputs)
    output:
        html = str(REVIEW_DIR / "review.html")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "generate_qc_review.py"),
        plate = PLATE,
        plate_dir = str(PLATE_DIR),
        plots_dir = str(REVIEW_DIR / "plots"),
        config_file = str(PIPELINE_DIR / "config.yaml"),
        outdir = str(REVIEW_DIR),
        logdir = str(LOG_DIR),
        decisions_path = DECISIONS_FILE,
        subplates_arg = ("--subplates " + " ".join(SUBPLATES)) if IS_384 else "",
        layout_arg = (f"--layout-tsv {config['plate_layout_tsv']}"
                      if IS_384 and config.get("plate_layout_tsv") else ""),
        cellenone_arg = (f"--cellenone-dir {CELLENONE_DIR}") if CELL_IMAGES_ENABLED else "",
        assets_arg = f"--assets-mode {ASSETS_MODE}",
        wells = " ".join(WELLS)
    log:
        str(LOG_DIR / f"qc_review_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = QC_REVIEW_RES.get("mem_mb", 8000) if IS_384 else 8000,
        runtime = 30 if not IS_384 else QC_REVIEW_RES.get("time", 180),
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        python {params.script} \
            --plate {params.plate} \
            --plate-dir {params.plate_dir} \
            --plots-dir {params.plots_dir} \
            --config {params.config_file} \
            --outdir {params.outdir} \
            --decisions-path {params.decisions_path} \
            {params.subplates_arg} {params.layout_arg} \
            {params.cellenone_arg} {params.assets_arg} \
            --wells {params.wells} \
            2>&1 | tee {log}
        """

rule validate_qc_decisions:
    """Validate the human-edited <PLATE_DIR>/qc_decisions.csv (run on demand, not part of `all`)."""
    input:
        decisions = DECISIONS_FILE
    output:
        flag = str(REVIEW_DIR / "qc_decisions.validated.flag")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "validate_qc_decisions.py"),
        subplates_arg = ("--subplates " + " ".join(SUBPLATES)) if IS_384 else "",
        wells = " ".join(WELLS)
    log:
        str(LOG_DIR / f"qc_decisions_validate_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 10,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {REVIEW_DIR} {LOG_DIR}
        python {params.script} \
            --decisions {input.decisions} \
            --wells {params.wells} \
            {params.subplates_arg} \
            --out-flag {output.flag} \
            2>&1 | tee {log}
        """

# ------------------------------------------------------------------------
# Two-pass AneuFinder: post-review (second) pass on PASS wells only
# ------------------------------------------------------------------------

rule derive_included_wells:
    """Derive the PASS (+ optional REVIEW) wells from a validated decisions CSV."""
    input:
        flag = str(REVIEW_DIR / "qc_decisions.validated.flag"),   # ensures validation ran
        decisions = DECISIONS_FILE
    output:
        included = str(INCLUDED_WELLS_TSV)
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "derive_included_wells.py"),
        include_review = INCLUDE_REVIEW_FLAG
    log:
        str(LOG_DIR / f"derive_included_wells_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = 2000,
        runtime = 10,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {REVIEW_DIR} {LOG_DIR}
        python {params.script} \
            --decisions {input.decisions} {params.include_review} \
            --out {output.included} \
            2>&1 | tee {log}
        """

rule run_aneufinder_reviewed:
    """Second AneuFinder pass on PASS wells only.

    Reuses run_aneufinder.R UNCHANGED by pointing --input at a directory of symlinks
    to the included wells' BAMs and --output at a separate aneufinder_reviewed/ dir.
    Original BAMs and the first-pass aneufinder/ outputs are never touched.
    """
    input:
        included = str(INCLUDED_WELLS_TSV),
        bams = expand(SUB_BASE + "/bam/{well}.bam", sub=SUBPLATES, well=WELLS),
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        gc_rds = str(PLATE_GC_RDS)
    output:
        flag = str(REVIEWED_DIR / "complete.flag")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "run_aneufinder.R"),
        stager = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "stage_aneufinder_input.py"),
        bam_dir_args = BAM_DIR_ARGS,
        # lambda-wrapped: a params STRING containing `{well}` would otherwise be
        # treated as a wildcard reference and fail to expand in a wildcard-free rule.
        id_template = lambda wildcards: ID_TEMPLATE,
        logdir = str(LOG_DIR),
        link_dir = str(REVIEWED_DIR / "input_bams"),
        outdir = str(REVIEWED_DIR),
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        gc_rds = str(PLATE_GC_RDS),
        method = ",".join(config["aneufinder"]["method"]),
        binsize = config["aneufinder"]["binsize"],
        chromosomes = ",".join(config["aneufinder"]["chromosomes"]),
        num_cpu = config["aneufinder"]["num_cpu"],
        cluster_plots = "TRUE" if config["aneufinder"]["cluster_plots"] else "FALSE"
    log:
        str(LOG_DIR / f"aneufinder_reviewed_{PLATE}.log")
    threads: config["resources"].get("aneufinder_reviewed", config["resources"]["aneufinder"])["threads"]
    resources:
        mem_mb = config["resources"].get("aneufinder_reviewed", config["resources"]["aneufinder"])["mem_mb"],
        runtime = config["resources"].get("aneufinder_reviewed", config["resources"]["aneufinder"])["time"],
        partition = config["resources"].get("aneufinder_reviewed", config["resources"]["aneufinder"])["partition"]
    shell:
        r"""
        mkdir -p {params.link_dir} {params.outdir} {params.logdir}
        # Symlink-only input set for the included wells; never copy or mutate the
        # originals. --min-reads 0 because the read-count gate already happened at
        # review time and the human call wins. The stager understands both the
        # 2-column and the 384-mode 3-column included_wells.tsv (the shell loop this
        # replaced silently read the SUBPLATE column as the well name).
        python {params.stager} \
            {params.bam_dir_args} \
            --id-template '{params.id_template}' \
            --wells-file {input.included} \
            --min-reads 0 \
            --out-dir {params.link_dir} \
            --report {params.outdir}/input_bams.tsv \
            --label "Second-pass AneuFinder" \
            2>&1 | tee {log}
        Rscript {params.script} \
            --input {params.link_dir} \
            --output {params.outdir} \
            --blacklist {params.blacklist} \
            --gc-rds {params.gc_rds} \
            --method {params.method} \
            --binsize {params.binsize} \
            --chromosomes {params.chromosomes} \
            --numcpu {params.num_cpu} \
            --cluster-plots {params.cluster_plots} \
            --reuse-existing FALSE \
            2>&1 | tee -a {log}
        touch {output.flag}
        """

rule render_reviewed_profiles:
    """Render second-pass per-well PNGs + a genome-wide CN heatmap for the CN viewer."""
    input:
        flag = str(REVIEWED_DIR / "complete.flag")
    output:
        manifest = str(REVIEWED_PLOTS_DIR / "manifest.json"),
        heatmap = str(REVIEWED_HEATMAP)
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "render_well_profiles.R"),
        models_dir = str(REVIEWED_DIR / "MODELS"),
        outdir = str(REVIEWED_PLOTS_DIR),
        logdir = str(LOG_DIR),
        method = REVIEW_METHOD
    log:
        str(LOG_DIR / f"cn_review_plots_{PLATE}.log")
    threads: 4
    resources:
        mem_mb = QC_REVIEW_RES["mem_mb"],
        runtime = QC_REVIEW_RES["time"],
        partition = config["resources"]["default"]["partition"]
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        Rscript {params.script} \
            --input {params.models_dir} \
            --outdir {params.outdir} \
            --method {params.method} \
            --format svg \
            --heatmap {output.heatmap} \
            2>&1 | tee {log}
        """

rule cn_review:
    """Build the final copy-number review viewer (read-only) over the second pass."""
    input:
        unpack(lambda wc: dict(
            {"manifest": str(REVIEWED_PLOTS_DIR / "manifest.json"),
             "heatmap": str(REVIEWED_HEATMAP),
             "included": str(INCLUDED_WELLS_TSV),
             "multiqc": expand(SUB_BASE + "/multiqc/multiqc_report.html", sub=SUBPLATES)},
            **({"cellenone": str(CELLENONE_DIR / "complete.flag")} if CELL_IMAGES_ENABLED else {})
        ))
    output:
        html = str(CN_REVIEW_DIR / "cn_review.html")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "generate_qc_review.py"),
        plate = PLATE,
        plate_dir = str(PLATE_DIR),
        plots_dir = str(REVIEWED_PLOTS_DIR),
        config_file = str(PIPELINE_DIR / "config.yaml"),
        outdir = str(CN_REVIEW_DIR),
        logdir = str(LOG_DIR),
        subplates_arg = ("--subplates " + " ".join(SUBPLATES)) if IS_384 else "",
        layout_arg = (f"--layout-tsv {config['plate_layout_tsv']}"
                      if IS_384 and config.get("plate_layout_tsv") else ""),
        cellenone_arg = (f"--cellenone-dir {CELLENONE_DIR}") if CELL_IMAGES_ENABLED else "",
        assets_arg = f"--assets-mode {ASSETS_MODE}",
        wells = " ".join(WELLS)
    log:
        str(LOG_DIR / f"cn_review_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = QC_REVIEW_RES.get("mem_mb", 8000) if IS_384 else 8000,
        runtime = 30 if not IS_384 else QC_REVIEW_RES.get("time", 180),
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {params.logdir}
        python {params.script} \
            --report-kind cn \
            --plate {params.plate} \
            --plate-dir {params.plate_dir} \
            --plots-dir {params.plots_dir} \
            --included-wells {input.included} \
            --heatmap {input.heatmap} \
            --config {params.config_file} \
            --outdir {params.outdir} \
            {params.subplates_arg} {params.layout_arg} \
            {params.cellenone_arg} {params.assets_arg} \
            --wells {params.wells} \
            2>&1 | tee {log}
        """

# ========================================================================
# Utilities
# ========================================================================

rule clean:
    """Remove this plate's generated files.

    DELIBERATELY NON-RECURSIVE. In 384 mode only the PLATE-LEVEL outputs are removed;
    the `<plate>/<subplate>/` directories are left completely alone. Recursing would
    delete all 768 per-subplate BAMs — days of alignment — for what reads like a
    routine cleanup. Clean a subplate by running `clean` against it directly.
    """
    shell:
        """
        echo "Cleaning plate-level outputs in {PLATE_DIR} (subplate dirs are NOT touched)..."
        rm -rf {PLATE_DIR}/demux
        rm -rf {PLATE_DIR}/dedup
        rm -rf {PLATE_DIR}/filtered
        rm -rf {PLATE_DIR}/fastqc
        rm -rf {PLATE_DIR}/raw_bam
        rm -rf {PLATE_DIR}/bam
        rm -rf {PLATE_DIR}/aneufinder
        rm -rf {REVIEWED_DIR}
        rm -rf {PLATE_DIR}/multiqc
        rm -rf {PLATE_DIR}/qc_review
        rm -rf {PLATE_DIR}/CN_review
        rm -rf {PLATE_DIR}/cellenone
        rm -rf {PLATE_DIR}/mappability
        rm -rf {PLATE_DIR}/logs
        """

onsuccess:
    print("\n" + "="*60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*60)
    print(f"Plate: {PLATE} ({PLATE_FORMAT}-well)")
    print(f"Results in: {PLATE_DIR}")
    # MultiQC is per subplate, so 384 has one report per subplate rather than one
    # at the plate level.
    for _s in SUBPLATES:
        print(f"MultiQC report: {SUB_BASE.format(sub=_s)}/multiqc/multiqc_report.html")
    print(f"Logs in: {PLATE_DIR}/logs")

onerror:
    print("\n" + "="*60)
    print("PIPELINE ERROR")
    print("="*60)
    print(f"Plate: {PLATE}")
    print(f"Check logs in: {PLATE_DIR}/logs/")
