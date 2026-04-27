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
        str(PLATE_DIR / "multiqc" / "multiqc_report.html")

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

rule deduplicate:
    input:
        str(PLATE_DIR / "demux" / "demux_stats.json")
    output:
        summary = str(PLATE_DIR / "dedup" / "dedup_summary.tsv"),
        stats = str(PLATE_DIR / "dedup" / "dedup_stats.json"),
        fastqs = expand(str(PLATE_DIR / "dedup" / "{well}_R{read}.dedup.fastq.gz"),
                        well=WELLS, read=[1, 2])
    params:
        indir = str(PLATE_DIR / "demux"),
        outdir = str(PLATE_DIR / "dedup"),
        script = str(PIPELINE_DIR / "workflow" / "scripts" / "preprocessing" / "dedup.py"),
        method = config["preprocessing"]["dedup_method"],
        min_length = config["preprocessing"]["min_read_length"]
    log:
        str(LOG_DIR / f"dedup_{PLATE}.log")
    threads: config["resources"]["dedup"]["threads"]
    resources:
        mem_mb = config["resources"]["dedup"]["mem_mb"],
        runtime = config["resources"]["dedup"]["time"],
        partition = config["resources"]["dedup"]["partition"]
    conda:
        "workflow/envs/preprocessing.yaml"
    shell:
        """
        mkdir -p {PLATE_DIR}/dedup
        python {params.script} \
            --indir {params.indir} \
            --outdir {params.outdir} \
            --method {params.method} \
            --min-length {params.min_length} \
            --stats {output.stats} \
            --log {log} 2>&1
        """

rule filter_dimers:
    input:
        str(PLATE_DIR / "dedup" / "dedup_summary.tsv")
    output:
        summary = str(PLATE_DIR / "filtered" / "adapter_filter_summary.tsv"),
        stats = str(PLATE_DIR / "filtered" / "filter_stats.json"),
        fastqs = expand(str(PLATE_DIR / "filtered" / "{well}_R{read}.filtered.fastq.gz"),
                        well=WELLS, read=[1, 2])
    params:
        indir = str(PLATE_DIR / "dedup"),
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
            --suffix .dedup.fastq.gz \
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
        bam = str(PLATE_DIR / "bam" / "{well}.bam"),
        bai = str(PLATE_DIR / "bam" / "{well}.bam.bai"),
        stats = str(PLATE_DIR / "bam" / "{well}.stats.txt"),
        flagstat = str(PLATE_DIR / "bam" / "{well}.flagstat.txt")
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
        mkdir -p {PLATE_DIR}/bam {LOG_DIR}/align
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
        rm -rf {PLATE_DIR}/bam
        rm -rf {PLATE_DIR}/aneufinder
        rm -rf {PLATE_DIR}/multiqc
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
