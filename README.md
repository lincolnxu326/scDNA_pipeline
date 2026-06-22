# scDNA Pipeline

Snakemake workflow for processing Strand-seq or single-cell DNA sequencing plates from paired FASTQ input through demultiplexing, adapter filtering, FastQC, alignment, UMI-tools deduplication, blacklist diagnostics, and AneuFinder CNV calling.

This repository is the pipeline codebase only. Plate-specific data, intermediate files, and runtime logs are produced outside version control.

## Overview

The active workflow currently supports:

- plate-level FASTQ demultiplexing into well FASTQs
- adapter dimer filtering
- FastQC and MultiQC reporting
- per-well alignment with Bowtie2 and Samtools
- UMI-aware BAM deduplication with UMI-tools
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

Barcode lookup now follows this order:

1. `resources/barcodes.tsv`
2. `resources/barcodes/barcodes.tsv`
3. `<plate_dir>/barcodes.tsv`

This allows a shared barcode definition to be reused across multiple plate runs
without copying `barcodes.tsv` into every plate directory.

The exact file paths are otherwise controlled through `submit_pipeline.sh` and
`config.yaml`.

## Main Stages

1. Demultiplex plate FASTQs into per-well FASTQs.
2. Remove adapter dimer reads.
3. Run FastQC and aggregate with MultiQC.
4. Align per-well reads with Bowtie2 and Samtools.
5. Deduplicate aligned BAMs with UMI-tools.
6. Generate blacklist diagnostics from the configured mappability reference.
7. Run the **first** AneuFinder pass for CNV calling on all wells.
8. Build a static HTML QC review report for manual well decisions.
9. Run the **second** AneuFinder pass on the wells marked PASS, then build the
   final copy-number review viewer.

## Two-Pass AneuFinder Workflow

CNV calling runs in two passes with a human review step in between:

```
preprocessing → aneufinder (first pass) → review → [edit CSV] → validate
              → aneufinder_reviewed (second pass) → cn_review
```

**Why two passes.** The first pass runs AneuFinder on every well so a reviewer can
look at all profiles and decide which wells to keep. The second pass reruns
AneuFinder on **only the PASS wells** to produce the final copy-number outputs.
The two are kept in separate directories so the first-pass outputs the review was
based on are never overwritten, and the original BAMs are never modified (the
second pass reads them through a symlink-only input directory).

### Step 1 — first pass + review report

Run the pipeline up to and including the first-pass review report (this is the
former "full" mode; it stops for human review):

```bash
# submit_pipeline.sh: MODE="pre_review"   (or target `all`)
sbatch submit_pipeline.sh
```

This produces the self-contained, fully static report:

```
<PLATE_DIR>/qc_review/review.html
```

It needs no server and opens over Samba or as a local `file://` document — every
well's two AneuFinder plots (a **copy-number profile** and a **bin read-count
histogram**) plus metadata are embedded. Use it to:

1. Click a well to load its metadata, **profile**, and **histogram**, alongside its
   automatic read-count status.
2. Set a **decision** with the PASS / EXCLUDE / REVIEW / REPEAT buttons, pick any
   number of **reason** chips, add **notes**, and click **Save decision**. The cell
   recolours by decision and a dot marks wells you've saved.
3. Click **Copy terminal save command** and paste it into a terminal on the cluster
   — it writes your decisions straight to `<PLATE_DIR>/qc_decisions.csv` (per-plate,
   in the data directory). The browser cannot write into the project, so this command
   (a `cat > … <<'EOF'` heredoc) — or **Download qc_decisions.csv** placed at that path
   — is how the file gets saved.

### Step 2 — one post-review run

After the CSV is saved, a single mode runs everything downstream:

```bash
# submit_pipeline.sh: MODE="post_review"   (or target `post_review`)
sbatch submit_pipeline.sh
```

`post_review` validates the decisions, derives the PASS wells, runs the second
AneuFinder pass, renders the reviewed CN plots + genome heatmap, and builds the
final viewer. Outputs:

```
<PLATE_DIR>/aneufinder_reviewed/      # second-pass models + profiles (PASS wells)
<PLATE_DIR>/qc_review/cn_review.html  # final read-only CN viewer (profiles + genome heatmap)
```

So the normal operator sequence is just two runs:

1. `MODE=pre_review` → review HTML, then stop
2. save `<PLATE_DIR>/qc_decisions.csv` (use the report's **Copy terminal save command**)
3. `MODE=post_review` → final `cn_review.html`

**Debug / step-by-step modes.** `post_review` is the convenience target; the same
chain is also exposed as individual modes for partial reruns or debugging:
`validate_review` → `aneufinder_reviewed` → `cn_review` (see `submit_pipeline.sh`).

Only `PASS` wells are included by default; set `qc_review.include_review: true` in
`config.yaml` to also include `REVIEW` wells. `EXCLUDE` and `REPEAT` are excluded
(`REPEAT` wells are logged as "flagged for rerun").

`<PLATE_DIR>/qc_decisions.csv` is the reproducible human artefact (per plate, in the
data directory; path configurable via `qc_review.decisions_file`), and is the only
input the second pass gates on — a missing CSV never blocks the pre-review stages. All
profile and histogram PNGs (both passes) and the genome heatmap are rendered
read-only in R from the AneuFinder `.RData` models; the page order of
`profiles_*.pdf` is never used. See [`config/README.md`](config/README.md) for the
CSV schema and decision semantics, and `submit_pipeline.sh` for the full list of modes.

### Automatic gating (read-count only)

Each well's decision is pre-filled from an automatic status based **only on read
count**: `usable_reads` = mapped reads in the dedup BAM (`reads mapped` from
`bam/{well}.stats.txt`). Cutoffs are in `config.yaml → qc_review`:

| usable_reads | auto_status | default decision |
|--------------|-------------|------------------|
| `>= usable_reads_pass_cutoff` (100000) | PASS | PASS |
| `>= usable_reads_warn_cutoff` (50000), `<` pass | WARN | REVIEW |
| `< usable_reads_warn_cutoff` | FAIL | EXCLUDE |
| missing | UNKNOWN | REVIEW |

So wells are pre-sorted (not all `REVIEW`); the human decision always wins.
**Duplication is not used for gating** — the panel shows it only as a clearly
labelled diagnostic (**UMI-tools duplicate rate** from `dedup/dedup_summary.tsv`,
i.e. umi_tools `duplicate_reads/total`, plus **UMI dedup retention**), distinct from
FastQC/MultiQC duplication. Multiple review reasons are stored `;`-separated in the
single `reason` CSV column.

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

If a shared barcode table is used, place it in either:

- `resources/barcodes.tsv`
- `resources/barcodes/barcodes.tsv`

If neither shared location exists, the pipeline falls back to
`<plate_dir>/barcodes.tsv`.

## Outputs

Per-plate outputs are written into the selected plate directory. Typical generated folders include:

- `demux/`
- `dedup/`
- `filtered/`
- `fastqc/`
- `raw_bam/`
- `bam/`
- `mappability/`
- `aneufinder/`
- `aneufinder_reviewed/`
- `multiqc/`
- `qc_review/`
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
