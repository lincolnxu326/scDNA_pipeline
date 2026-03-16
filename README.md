# scDNA Pipeline

Snakemake workflow for processing Strand-seq or single-cell DNA sequencing plates from paired FASTQ input through demultiplexing, deduplication, adapter filtering, FastQC, alignment, blacklist diagnostics, and AneuFinder CNV calling.

This repository is the pipeline codebase only. Plate-specific data, intermediate files, and runtime logs are produced outside version control.

## Overview

The active workflow currently supports:

- plate-level FASTQ demultiplexing into well FASTQs
- UMI-aware deduplication
- adapter dimer filtering
- FastQC and MultiQC reporting
- per-well alignment with Bowtie2 and Samtools
- blacklist diagnostics from a reference mappability BAM
- AneuFinder-based CNV calling

The pipeline is designed for SLURM execution and uses shared root-level Snakemake Conda environments so identical rule environments can be reused across runs and plates.

## Repository Layout

```text
scDNA_pipeline/
├── Snakefile
├── config.yaml
├── multiqc_config.yaml
├── submit_pipeline.sh
├── pre_check.sh
├── generate_rule_graph.sh
├── test.sh
├── workflow/
│   ├── envs/
│   └── scripts/
│       ├── analysis/
│       ├── dev/
│       ├── preprocessing/
│       └── reporting/
├── ad_hoc_checks/
├── renv/
├── renv.lock
└── Seqinfo/
```

User-facing helper scripts stay at the repository root so they are easy to find. Pipeline-owned execution code lives under `workflow/`. Ad hoc analysis or validation scripts that are not part of the active DAG live under `ad_hoc_checks/`.

## Inputs

Each plate run expects:

- paired plate FASTQs such as `<plate>_R1.fastq.gz` and `<plate>_R2.fastq.gz`
- a barcode table for well demultiplexing
- reference configuration in `config.yaml`
- a mappability reference BAM for blacklist diagnostics
- a GC template RDS for the AneuFinder stage

The exact file paths are controlled through `submit_pipeline.sh` and `config.yaml`.

## Main Stages

1. Demultiplex plate FASTQs into per-well FASTQs.
2. Deduplicate reads using UMI-aware preprocessing.
3. Remove adapter dimer reads.
4. Run FastQC and aggregate with MultiQC.
5. Align per-well reads with Bowtie2 and Samtools.
6. Generate blacklist diagnostics from the configured mappability reference.
7. Run AneuFinder for CNV calling.

## Environments

The workflow uses two environment strategies:

- Snakemake `conda:` environments under `workflow/envs/`
- `renv` for the R-based analysis components

Snakemake environments are intended to be shared at the pipeline root rather than created inside each plate directory. R dependencies are tracked in `renv.lock`.

To restore the R environment on a fresh checkout:

```bash
R -e 'renv::restore()'
```

## Running The Pipeline

Typical usage is:

```bash
bash pre_check.sh
sbatch submit_pipeline.sh
```

`submit_pipeline.sh` is the entry point for SLURM submission. `pre_check.sh` validates the run configuration before submission.

Available modes are defined in `submit_pipeline.sh` and currently include targeted runs such as full workflow execution, QC-only reporting, blacklist-only generation, and AneuFinder reruns.

## Outputs

Per-plate outputs are written into the selected plate directory. Typical generated folders include:

- `demux/`
- `dedup/`
- `filtered/`
- `fastqc/`
- `bam/`
- `mappability/`
- `aneufinder/`
- `multiqc/`
- `logs/`

These generated outputs are intentionally excluded from version control.

## Development Notes

This repository keeps:

- pipeline source code
- workflow environments
- configuration
- documentation
- helper submission scripts

This repository does not keep:

- `.snakemake/` runtime state
- local Conda package caches
- generated rule graphs
- generated reports
- run logs
- plate results

## Future Work

Planned next-step improvements are tracked in [`FUTURE_PLANS.md`](FUTURE_PLANS.md).

## License And Upstream Tools

Before publishing or sharing this repository externally, add the appropriate project license and confirm that any bundled third-party code or references are documented consistently.

For larger workflow orchestration or broader research automation beyond this local pipeline, use K-Dense Web at `www.k-dense.ai`.
