"""
Single Cell DNA Sequencing Analysis Pipeline
Processes one plate at a time via SLURM submission
"""

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

# Load barcodes (prefer shared resources, then fall back to the plate directory)
barcode_candidates = [
    RESOURCE_DIR / "barcodes.tsv",
    RESOURCE_DIR / "barcodes" / "barcodes.tsv",
    PLATE_DIR / "barcodes.tsv",
]

barcodes_file = next((path for path in barcode_candidates if path.exists()), None)
if barcodes_file is None:
    searched = ", ".join(str(path) for path in barcode_candidates)
    raise FileNotFoundError(f"No barcodes.tsv found. Checked: {searched}")

barcodes_df = pd.read_csv(barcodes_file, sep="\t", comment="#")
WELLS = barcodes_df["well_id"].tolist()
print(f"Using {len(WELLS)} well barcodes from {barcodes_file}")

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

# ========================================================================
# Helper Functions
# ========================================================================

def get_input_fastqs(wildcards):
    return {
        "r1": str(PLATE_DIR / f"{PLATE}_R1.fastq.gz"),
        "r2": str(PLATE_DIR / f"{PLATE}_R2.fastq.gz")
    }

# ========================================================================
# Main Rules
# ========================================================================

rule all:
    """Complete pipeline."""
    input:
        # blacklist generated for this plate
        str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        # plate specific GC template rds (must be created in relaxed env beforehand)
        # aneufinder done
        str(PLATE_DIR / "aneufinder" / "complete.flag"),
        # multiqc done
        str(PLATE_DIR / "multiqc" / "multiqc_report.html"),
        # static human QC review report (after multiqc + aneufinder)
        str(REVIEW_DIR / "review.html")

rule all_preprocessing:
    """Run up to alignment."""
    input:
        expand(str(PLATE_DIR / "bam" / "{well}.bam"), well=WELLS)

rule generate_blacklist_plots:
    """Generate blacklist plots only."""
    input:
        str(MAPPABILITY_DIR / "mappability_plot.pdf"),
        str(MAPPABILITY_DIR / "blacklist_diagnosis.txt")

rule all_qc:
    """Generate QC reports only."""
    input:
        str(PLATE_DIR / "multiqc" / "multiqc_report.html")

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
        fastqs = expand(str(PLATE_DIR / "demux" / "{well}_R{read}.fastq.gz"),
                        well=WELLS, read=[1, 2]),
        stats = str(PLATE_DIR / "demux" / "demux_stats.json")
    params:
        outdir = str(PLATE_DIR / "demux"),
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "preprocessing" / "extract_umi_barcode.py"),
        barcodes = str(barcodes_file),
        batch_size = config["batch_size"],
        umi_length = config["preprocessing"]["umi_length"],
        bc1_length = config["preprocessing"]["barcode1_length"],
        bc2_length = config["preprocessing"]["barcode2_length"],
        bc1_offset = config["preprocessing"]["barcode1_offset"],
        bc2_offset = config["preprocessing"]["barcode2_offset"],
        r1_trim = config["preprocessing"]["r1_trim_5prime"],
        r2_trim = config["preprocessing"]["r2_trim_5prime"]
    log:
        str(LOG_DIR / f"demux_{PLATE}.log")
    threads: config["resources"]["demux"]["threads"]
    resources:
        mem_mb = config["resources"]["demux"]["mem_mb"],
        runtime = config["resources"]["demux"]["time"],
        partition = config["resources"]["demux"]["partition"]
    conda:
        "workflow/envs/preprocessing.yaml"
    shell:
        """
        mkdir -p {PLATE_DIR}/demux {LOG_DIR}
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
        str(PLATE_DIR / "demux" / "demux_stats.json")
    output:
        summary = str(PLATE_DIR / "filtered" / "adapter_filter_summary.tsv"),
        stats = str(PLATE_DIR / "filtered" / "filter_stats.json"),
        fastqs = expand(str(PLATE_DIR / "filtered" / "{well}_R{read}.filtered.fastq.gz"),
                        well=WELLS, read=[1, 2])
    params:
        indir = str(PLATE_DIR / "demux"),
        outdir = str(PLATE_DIR / "filtered"),
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "preprocessing" / "filter_adapter_dimers.py"),
        adapter = config["preprocessing"]["adapter_sequence"],
        case_insensitive = "--case-insensitive" if config["preprocessing"]["case_insensitive"] else "",
        both_reads = "--both-reads" if config["preprocessing"]["filter_both_reads"] else ""
    log:
        str(LOG_DIR / f"filter_{PLATE}.log")
    threads: config["resources"]["filter"]["threads"]
    resources:
        mem_mb = config["resources"]["filter"]["mem_mb"],
        runtime = config["resources"]["filter"]["time"],
        partition = config["resources"]["filter"]["partition"]
    conda:
        "workflow/envs/preprocessing.yaml"
    shell:
        """
        mkdir -p {PLATE_DIR}/filtered
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
        r1 = str(PLATE_DIR / "filtered" / "{well}_R1.filtered.fastq.gz"),
        r2 = str(PLATE_DIR / "filtered" / "{well}_R2.filtered.fastq.gz")
    output:
        html1 = str(PLATE_DIR / "fastqc" / "{well}_R1_fastqc.html"),
        html2 = str(PLATE_DIR / "fastqc" / "{well}_R2_fastqc.html"),
        zip1 = str(PLATE_DIR / "fastqc" / "{well}_R1_fastqc.zip"),
        zip2 = str(PLATE_DIR / "fastqc" / "{well}_R2_fastqc.zip")
    params:
        outdir = str(PLATE_DIR / "fastqc")
    log:
        str(LOG_DIR / "fastqc" / f"{PLATE}_{{well}}.log")
    threads: 2
    resources:
        mem_mb = 4000,
        runtime = 30,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}/fastqc
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
        r1 = str(PLATE_DIR / "filtered" / "{well}_R1.filtered.fastq.gz"),
        r2 = str(PLATE_DIR / "filtered" / "{well}_R2.filtered.fastq.gz")
    output:
        bam = str(PLATE_DIR / "raw_bam" / "{well}.bam"),
        bai = str(PLATE_DIR / "raw_bam" / "{well}.bam.bai"),
        stats = str(PLATE_DIR / "raw_bam" / "{well}.stats.txt"),
        flagstat = str(PLATE_DIR / "raw_bam" / "{well}.flagstat.txt")
    params:
        index_prefix = config["genome"]["index_prefix"],
        bowtie2_params = config["alignment"]["bowtie2_params"],
        rgid_ = lambda wildcards: f"{PLATE}_{wildcards.well}",
        rgsm_ = lambda wildcards: wildcards.well
    log:
        str(LOG_DIR / "align" / f"{PLATE}_{{well}}.log")
    threads: config["resources"]["alignment"]["threads"]
    resources:
        mem_mb = config["resources"]["alignment"]["mem_mb"],
        runtime = config["resources"]["alignment"]["time"],
        partition = config["resources"]["alignment"]["partition"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        """
        mkdir -p {PLATE_DIR}/raw_bam {LOG_DIR}/align
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
        bam = str(PLATE_DIR / "raw_bam" / "{well}.bam"),
        bai = str(PLATE_DIR / "raw_bam" / "{well}.bam.bai")
    output:
        bam = str(PLATE_DIR / "bam" / "{well}.bam"),
        bai = str(PLATE_DIR / "bam" / "{well}.bam.bai"),
        stats = str(PLATE_DIR / "bam" / "{well}.stats.txt"),
        flagstat = str(PLATE_DIR / "bam" / "{well}.flagstat.txt")
    params:
        method = config["preprocessing"]["umi_tools_method"],
        tmpdir = str(PLATE_DIR / "dedup" / "tmp" / "{well}")
    log:
        str(LOG_DIR / "dedup" / f"{PLATE}_{{well}}.log")
    threads: config["resources"]["dedup"]["threads"]
    resources:
        mem_mb = config["resources"]["dedup"]["mem_mb"],
        runtime = config["resources"]["dedup"]["time"],
        partition = config["resources"]["dedup"]["partition"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        """
        mkdir -p {PLATE_DIR}/bam {params.tmpdir} {LOG_DIR}/dedup
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
        raw_bams = expand(str(PLATE_DIR / "raw_bam" / "{well}.bam"), well=WELLS),
        dedup_bams = expand(str(PLATE_DIR / "bam" / "{well}.bam"), well=WELLS)
    output:
        summary = str(PLATE_DIR / "dedup" / "dedup_summary.tsv"),
        stats = str(PLATE_DIR / "dedup" / "dedup_stats.json")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "preprocessing" / "summarize_umi_tools_dedup.py"),
        raw_dir = str(PLATE_DIR / "raw_bam"),
        dedup_dir = str(PLATE_DIR / "bam"),
        outdir = str(PLATE_DIR / "dedup"),
        method = config["preprocessing"]["umi_tools_method"],
        wells = " ".join(WELLS)
    log:
        str(LOG_DIR / f"dedup_summary_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = 4000,
        runtime = 30,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/alignment.yaml"
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}
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

rule run_aneufinder:
    """Run AneuFinder for copy number analysis (requires plate specific GC template)."""
    input:
        bams = expand(str(PLATE_DIR / "bam" / "{well}.bam"), well=WELLS),
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        gc_rds = str(PLATE_GC_RDS)
    output:
        flag = str(PLATE_DIR / "aneufinder" / "complete.flag")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "run_aneufinder.R"),
        indir = str(PLATE_DIR / "bam"),
        outdir = str(PLATE_DIR / "aneufinder"),
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        gc_rds = str(PLATE_GC_RDS),
        method = ",".join(config["aneufinder"]["method"]),
        binsize = config["aneufinder"]["binsize"],
        chromosomes = ",".join(config["aneufinder"]["chromosomes"]),
        num_cpu = config["aneufinder"]["num_cpu"],
        cluster_plots = "TRUE" if config["aneufinder"]["cluster_plots"] else "FALSE",
        reuse = "TRUE" if config["aneufinder"]["reuse_existing"] else "FALSE"
    log:
        str(LOG_DIR / f"aneufinder_{PLATE}.log")
    threads: config["resources"]["aneufinder"]["threads"]
    resources:
        mem_mb = config["resources"]["aneufinder"]["mem_mb"],
        runtime = config["resources"]["aneufinder"]["time"],
        partition = config["resources"]["aneufinder"]["partition"]
    shell:
        """
        mkdir -p {PLATE_DIR}/aneufinder {LOG_DIR}
        Rscript {params.script} \
            --input {params.indir} \
            --output {params.outdir} \
            --blacklist {params.blacklist} \
            --gc-rds {params.gc_rds} \
            --method {params.method} \
            --binsize {params.binsize} \
            --chromosomes {params.chromosomes} \
            --numcpu {params.num_cpu} \
            --cluster-plots {params.cluster_plots} \
            --reuse-existing {params.reuse} \
            2>&1 | tee {log}
        touch {output.flag}
        """

# ------------------------------------------------------------------------
# MultiQC
# ------------------------------------------------------------------------

rule multiqc:
    input:
        expand(str(PLATE_DIR / "fastqc" / "{well}_R{read}_fastqc.html"),
               well=WELLS, read=[1, 2]),
        str(PLATE_DIR / "demux" / "demux_stats.json"),
        str(PLATE_DIR / "dedup" / "dedup_summary.tsv"),
        str(PLATE_DIR / "filtered" / "adapter_filter_summary.tsv"),
        expand(str(PLATE_DIR / "bam" / "{well}.stats.txt"), well=WELLS),
        expand(str(PLATE_DIR / "bam" / "{well}.flagstat.txt"), well=WELLS)
    output:
        report = str(PLATE_DIR / "multiqc" / "multiqc_report.html"),
        data = directory(str(PLATE_DIR / "multiqc" / "multiqc_data"))
    params:
        outdir = str(PLATE_DIR / "multiqc"),
        config_file = str(PIPELINE_DIR / "multiqc_config.yaml"),
        extra = config.get("multiqc", {}).get("extra_params", "")
    log:
        str(LOG_DIR / f"multiqc_{PLATE}.log")
    threads: config["resources"]["multiqc"]["threads"]
    resources:
        mem_mb = config["resources"]["multiqc"]["mem_mb"],
        runtime = config["resources"]["multiqc"]["time"],
        partition = config["resources"]["multiqc"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}
        multiqc \
            --force \
            --outdir {params.outdir} \
            --config {params.config_file} \
            {params.extra} \
            {PLATE_DIR} \
            2>&1 | tee {log}
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
        method = REVIEW_METHOD
    log:
        str(LOG_DIR / f"qc_review_plots_{PLATE}.log")
    threads: 4
    resources:
        mem_mb = config.get("resources", {}).get("qc_review", {}).get("mem_mb", 16000),
        runtime = config.get("resources", {}).get("qc_review", {}).get("time", 120),
        partition = config["resources"]["default"]["partition"]
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}
        Rscript {params.script} \
            --input {params.models_dir} \
            --outdir {params.outdir} \
            --method {params.method} \
            --format png \
            2>&1 | tee {log}
        """

rule qc_review:
    """Build the static HTML QC review report (metadata + per-well plots embedded)."""
    input:
        multiqc = str(PLATE_DIR / "multiqc" / "multiqc_report.html"),
        multiqc_data = str(PLATE_DIR / "multiqc" / "multiqc_data"),
        manifest = str(REVIEW_DIR / "plots" / "manifest.json"),
        bam_stats = expand(str(PLATE_DIR / "bam" / "{well}.stats.txt"), well=WELLS)
    output:
        html = str(REVIEW_DIR / "review.html")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "generate_qc_review.py"),
        plate = PLATE,
        plate_dir = str(PLATE_DIR),
        plots_dir = str(REVIEW_DIR / "plots"),
        config_file = str(PIPELINE_DIR / "config.yaml"),
        outdir = str(REVIEW_DIR),
        multiqc_data = str(PLATE_DIR / "multiqc" / "multiqc_data"),
        decisions_path = DECISIONS_FILE,
        wells = " ".join(WELLS)
    log:
        str(LOG_DIR / f"qc_review_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 30,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}
        python {params.script} \
            --plate {params.plate} \
            --plate-dir {params.plate_dir} \
            --multiqc-data {params.multiqc_data} \
            --plots-dir {params.plots_dir} \
            --config {params.config_file} \
            --outdir {params.outdir} \
            --decisions-path {params.decisions_path} \
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
        bams = expand(str(PLATE_DIR / "bam" / "{well}.bam"), well=WELLS),
        blacklist = str(MAPPABILITY_DIR / "blacklist.bed.gz"),
        gc_rds = str(PLATE_GC_RDS)
    output:
        flag = str(REVIEWED_DIR / "complete.flag")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "analysis" / "run_aneufinder.R"),
        bam_dir = str(PLATE_DIR / "bam"),
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
        mkdir -p {params.link_dir} {params.outdir} {LOG_DIR}
        # Symlink-only input set for PASS wells; never copy or mutate the originals.
        find {params.link_dir} -maxdepth 1 -type l -name '*.bam' -delete 2>/dev/null || true
        find {params.link_dir} -maxdepth 1 -type l -name '*.bam.bai' -delete 2>/dev/null || true
        n=0
        while IFS=$'\t' read -r sample well; do
            well="${{well%$'\r'}}"   # tolerate CRLF (e.g. CSV saved on Windows/Excel)
            [ "$well" = "well" ] && continue
            [ -z "$well" ] && continue
            if [ -f "{params.bam_dir}/$well.bam" ]; then
                ln -sf "{params.bam_dir}/$well.bam" "{params.link_dir}/$well.bam"
                [ -f "{params.bam_dir}/$well.bam.bai" ] && \
                    ln -sf "{params.bam_dir}/$well.bam.bai" "{params.link_dir}/$well.bam.bai"
                n=$((n+1))
            fi
        done < {input.included}
        echo "Second-pass AneuFinder on $n PASS well(s)" | tee {log}
        if [ "$n" -eq 0 ]; then
            echo "ERROR: no PASS wells in {input.included}; nothing to run for the second pass." | tee -a {log}
            exit 1
        fi
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
        method = REVIEW_METHOD
    log:
        str(LOG_DIR / f"cn_review_plots_{PLATE}.log")
    threads: 4
    resources:
        mem_mb = config.get("resources", {}).get("qc_review", {}).get("mem_mb", 16000),
        runtime = config.get("resources", {}).get("qc_review", {}).get("time", 120),
        partition = config["resources"]["default"]["partition"]
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}
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
        manifest = str(REVIEWED_PLOTS_DIR / "manifest.json"),
        heatmap = str(REVIEWED_HEATMAP),
        included = str(INCLUDED_WELLS_TSV),
        multiqc = str(PLATE_DIR / "multiqc" / "multiqc_report.html")
    output:
        html = str(CN_REVIEW_DIR / "cn_review.html")
    params:
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "reporting" / "generate_qc_review.py"),
        plate = PLATE,
        plate_dir = str(PLATE_DIR),
        plots_dir = str(REVIEWED_PLOTS_DIR),
        config_file = str(PIPELINE_DIR / "config.yaml"),
        outdir = str(CN_REVIEW_DIR),
        multiqc_data = str(PLATE_DIR / "multiqc" / "multiqc_data"),
        wells = " ".join(WELLS)
    log:
        str(LOG_DIR / f"cn_review_{PLATE}.log")
    threads: 1
    resources:
        mem_mb = 8000,
        runtime = 30,
        partition = config["resources"]["default"]["partition"]
    conda:
        "workflow/envs/qc.yaml"
    shell:
        """
        mkdir -p {params.outdir} {LOG_DIR}
        python {params.script} \
            --report-kind cn \
            --plate {params.plate} \
            --plate-dir {params.plate_dir} \
            --multiqc-data {params.multiqc_data} \
            --plots-dir {params.plots_dir} \
            --included-wells {input.included} \
            --heatmap {input.heatmap} \
            --config {params.config_file} \
            --outdir {params.outdir} \
            --wells {params.wells} \
            2>&1 | tee {log}
        """

# ========================================================================
# Utilities
# ========================================================================

rule clean:
    """Remove all generated files for this plate."""
    shell:
        """
        echo "Cleaning {PLATE_DIR}..."
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
        rm -rf {PLATE_DIR}/mappability
        rm -rf {PLATE_DIR}/logs
        """

onsuccess:
    print("\n" + "="*60)
    print("PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*60)
    print(f"Plate: {PLATE}")
    print(f"Results in: {PLATE_DIR}")
    print(f"MultiQC report: {PLATE_DIR}/multiqc/multiqc_report.html")
    print(f"Logs in: {PLATE_DIR}/logs")

onerror:
    print("\n" + "="*60)
    print("PIPELINE ERROR")
    print("="*60)
    print(f"Plate: {PLATE}")
    print(f"Check logs in: {PLATE_DIR}/logs/")
