# scDNA Pipeline

Snakemake workflow for single-cell DNA sequencing plates. It takes paired plate FASTQs
through well demultiplexing, adapter-dimer filtering, alignment, UMI-aware
deduplication and AneuFinder copy-number calling. A human review step sits between a
first and a second AneuFinder pass. CellenONE dispenser images can be shown alongside
each well during review.

The repository holds code only. All plate data, intermediate files and results live in
the plate directory you point the pipeline at, outside version control.

## Contents

1. [Quick start](#quick-start)
2. [What the pipeline does](#what-the-pipeline-does)
3. [Inputs](#inputs)
4. [CellenONE run folder](#cellenone-run-folder)
5. [Outputs](#outputs)
6. [Running modes](#running-modes)
7. [Review workflow](#review-workflow)
8. [384-well plates](#384-well-plates)
9. [Configuration reference](#configuration-reference)
10. [Environments](#environments)
11. [Repository layout](#repository-layout)
12. [Further documentation](#further-documentation)

## Quick start

```bash
# 1. Check the setup
bash pre_check.sh

# 2. First pass: preprocessing, AneuFinder on all wells, review report. Stops here.
MODE=pre_review PLATE_DIR=/path/to/plate17_4 sbatch submit_pipeline.sh

# 3. Open <PLATE_DIR>/qc_review/review.html, make decisions, save qc_decisions.csv

# 4. Second pass: AneuFinder on PASS wells, final copy-number viewer
MODE=post_review PLATE_DIR=/path/to/plate17_4 sbatch submit_pipeline.sh
```

For a 384-well plate add `PLATE_FORMAT=384` and point `PLATE_DIR` at the plate
directory that holds the four subplates:

```bash
PLATE_FORMAT=384 MODE=pre_review PLATE_DIR=/path/to/384_well/plate21 sbatch submit_pipeline.sh
```

To run several plates at the same time from this directory, add
`SNAKEMAKE_EXTRA=--nolock`. Snakemake locks the pipeline directory, but each plate
writes to its own `PLATE_DIR`, so the runs do not collide.

## What the pipeline does

| # | Stage | Tool | Output folder |
|---|-------|------|---------------|
| 1 | Demultiplex plate FASTQs into well FASTQs, extract UMI | Python | `demux/` |
| 2 | Remove adapter-dimer reads | Python | `filtered/` |
| 3 | Read QC | FastQC, MultiQC | `fastqc/`, `multiqc/` |
| 4 | Align each well | Bowtie2, Samtools | `raw_bam/` |
| 5 | UMI-aware deduplication | UMI-tools (`directional`) | `bam/`, `dedup/` |
| 6 | Blacklist diagnostics from a mappability reference | R | `mappability/` |
| 7 | First AneuFinder pass, all wells | AneuFinder (`edivisive`, 1 Mb bins) | `aneufinder/` |
| 8 | Static HTML review report | Python, R | `qc_review/` |
| 9 | CellenONE image ingest and rendering (optional) | Python | `cellenone/` |
| 10 | Human review, saved as a CSV | browser | `qc_decisions.csv` |
| 11 | Second AneuFinder pass, PASS wells only | AneuFinder | `aneufinder_reviewed/` |
| 12 | Final copy-number viewer | Python, R | `CN_review/` |

Stages 1 to 9 run under `MODE=pre_review`. Stages 11 and 12 run under
`MODE=post_review`.

## Inputs

### Plate directory (96-well)

One directory per 96-well plate. The FASTQ names must match the directory name,
because the pipeline reads `<PLATE_DIR>/<dir name>_R1.fastq.gz`.

```text
plate17_4/
├── plate17_4_R1.fastq.gz      required
├── plate17_4_R2.fastq.gz      required
└── barcodes.tsv               required unless a shared table exists (see below)
```

Symlinks to FASTQs stored elsewhere are fine.

### Plate directory (384-well)

A 384-well plate is four 96-well subplates, each sequenced as its own FASTQ pair.
Subplate directories are named `<plate>_1` to `<plate>_4` and each follows the
96-well layout above.

```text
plate21/
├── plate21_1/
│   ├── plate21_1_R1.fastq.gz
│   ├── plate21_1_R2.fastq.gz
│   └── barcodes.tsv
├── plate21_2/  (same)
├── plate21_3/  (same)
├── plate21_4/  (same)
└── <CellenONE .Run folder>    optional, may also live elsewhere
```

Naming convention: lowercase `plateNN` for the plate and `plateNN_k` for the
subplates. Subplates are found automatically as `<plate>_<n>` directories that hold
a FASTQ pair. If the names do not follow that pattern, name the subplates
explicitly. For example, `384_well/plate17_umi/` holds `plate17_1` to `plate17_4`:

```bash
PLATE_FORMAT=384 PLATE_DIR=/path/to/384_well/plate17_umi \
SUBPLATES="plate17_1 plate17_2 plate17_3 plate17_4" sbatch submit_pipeline.sh
```

`SUBPLATES` sets the `subplates:` list in `config.yaml`, and the submit script checks
that each directory exists.

### Read structure

Set under `preprocessing:` in `config.yaml`. These values come from the original
lab scripts and should not change unless the library chemistry changes.

| Read | Layout |
|------|--------|
| R1 | 8 bp well barcode, 14 bp linker, then genomic sequence (22 bp trimmed) |
| R2 | 12 bp UMI, 8 bp well barcode, 14 bp linker, then genomic sequence (34 bp trimmed) |

Wells are assigned by the R1 barcode.

### barcodes.tsv

Tab-separated, with a header, one row per well:

```text
well_id	barcode
W01	TCTCATCG
W02	CCAACAGT
...
W96	...
```

The pipeline looks for the table in this order and uses the first one found:

1. `resources/barcodes.tsv` (shared, in this repo)
2. `resources/barcodes/barcodes.tsv` (shared, in this repo)
3. `<PLATE_DIR>/barcodes.tsv` (per plate; per subplate in 384 mode)

In 384 mode all four subplates must list `well_id` in the same order.

### Reference files

Set in `config.yaml`. All must exist before a run.

| Key | What it is |
|-----|------------|
| `genome.fasta` | Reference FASTA. Chromosomes must be named `chr1` to `chr22`, `chrX`, `chrY`. |
| `genome.index_prefix` | Bowtie2 index prefix for that FASTA |
| `mappability.reference_bam` | Euploid reference BAM used for blacklist diagnostics |
| GC template RDS | `resources/reference/hg38_binsize1000000_variable_bins_with_GC.rds`, or set `aneufinder.gc_rds` |

The GC template must match `aneufinder.binsize`. `run_aneufinder.R` is set up for
1 Mb bins.

## CellenONE run folder

This step is optional and display-only. It never changes a well's automatic status or
its default decision. It adds a cell-image panel to the review report so a reviewer
can see what the dispenser put in each well.

### Enabling it

Map the plate name to its CellenONE `.Run` folder in `config.yaml`:

```yaml
cellenone:
  enable: true
  runs:
    plate21: "/path/to/P21_22/plate_2/K1563_plate_2_20260707_135300_812.Run"
    plate24: "/path/to/384_well/plate24/P1_K1570_V_20260902_142546_324.Run"
```

The key is the `PLATE_DIR` name (for 384 mode, the plate, not a subplate). The
pipeline never guesses this mapping, because CellenONE names do not match ours
(CellenONE `plate_2` is our `plate21`). A plate with no entry simply runs without
cell images.

### What the folder must contain

Point the config at the `.Run` folder itself, exactly as the CellenONE software wrote
it. Do not rename files inside it.

| File | Required | Used for |
|------|----------|----------|
| `Reordered_*_isolated.xls` | **yes** | Per-well table: position, diameter, shape, intensity, and the run's detection and isolation criteria. Falls back to `*_isolated.xls` if no `Reordered_` copy exists. |
| `*_Printed_*_(<POS>)_Trans_*.png` | **yes** | Transmission image of each printed drop |
| `*_Printed_*_(<POS>)Blue_*.png`, `Green`, `Orange`, `Red` | optional | Fluorescence images. Only channels the run recorded are used. |
| `Reordered_*_geoprops.xls` | optional | Run-level statistics (drops attempted, isolation rate) |
| `*BackgroundEjZone*.png` | optional | Its width gives the purple ejection line |
| `cellenREPORT/Images_iso_det/Trans_*.jpg` | optional | Annotated frames; the green and purple lines are detected from them |
| `*.par`, `Tscatter.xls`, `Fscatter.xls`, `Clonality/`, `*.log` | not used | Left in place, ignored |

Notes on these files:

- The `.xls` files are tab-separated text, not Excel. Open them in a text editor if
  you need to check them.
- Images are matched to wells by the `(<POS>)` token in the filename, for example
  `(K-22)`. They are never matched by count or order. A run often has a few extra
  or missing images, and that is fine.
- `<POS>` uses the 384-plate grid (`A-1` to `P-24`) for 384 runs and the 96-plate
  grid (`A-1` to `H-12`) for 96 runs. In 96 mode `W01` is `A1`, `W02` is `B1`, and
  so on down each column.
- Detection and isolation thresholds are read from the run's own table, per well.
  They differ between runs and can change partway through a run (plate22 used three
  settings).
- If neither the annotated JPEGs nor the `BackgroundEjZone` image is present, the
  ejection lines fall back to `cellenone.lines.green_default` (396 px) and
  `purple_default` (647 px).

### Image call

The transmission image is scored against the two ejection lines. The nozzle is on
the left; objects to the right of the purple line are still up the capillary and are
not ejected.

| Condition | Call |
|-----------|------|
| exactly 1 object | `SINGLE` |
| 2 or more objects, rightmost right of the purple line | `PASS` |
| rightmost object between the green and purple lines | `CONTAMINATION` |
| all objects left of the green line | `FAIL` |
| 0 objects | `NO_OBJECT` |
| no image for this well | `NO_IMAGE` |

This image `PASS`/`FAIL` is not the sequencing `PASS`/`FAIL`. See
[`docs/CELLENONE_AND_QC.md`](docs/CELLENONE_AND_QC.md) for how the two differ and
where each number comes from.

### CellenONE outputs

Written to `<PLATE_DIR>/cellenone/`:

| File | Contents |
|------|----------|
| `wells_raw.tsv` | One row per well: CellenONE measurements and that well's own criteria |
| `objects.tsv` | One row per detected object: position, diameter, category |
| `cellenone_wells.tsv` | Per-well summary the report reads: object counts, image call, image paths |
| `run_meta.json` | Ejection lines and their source, criteria, channel offsets, pixel scale, run stats |
| `images/` | Rendered JPEGs, one per well per channel |

Rebuild only this layer with `MODE=cellenone`.

### Plates with no cell images

Nothing needs switching off. A plate with no entry in `cellenone.runs` (for example
plate17, plate19, and every plate dispensed before CellenONE) never builds a
`cellenone/` directory. The review report then leaves out the image panel and the
image call, and the per-well evidence is the copy-number profile and the bin
read-count histogram. The read gate, decisions and CSV are unchanged.

The same layout is used if a run folder is configured but gives no usable image; the
report logs a warning instead of failing. Set `qc_review.cell_images: false` to force
this layout for a plate that does have images.

## Outputs

Everything is written inside `PLATE_DIR`. Nothing is written into the repo except
logs under `logs/`.

### 96-well plate

```text
plate17_4/
├── demux/                 per-well FASTQs and demultiplexing stats
├── filtered/              per-well FASTQs after adapter-dimer removal
├── fastqc/                FastQC per well
├── multiqc/               MultiQC report for the plate
├── raw_bam/               aligned BAMs before deduplication
├── bam/                   deduplicated BAMs: W01.bam, .bai, .stats.txt, .flagstat.txt
├── dedup/                 UMI-tools logs and dedup_summary.tsv
├── mappability/           blacklist diagnostic plots
├── aneufinder/            first pass, all wells
│   ├── MODELS/            per-well .RData models
│   ├── profiles_edivisive.pdf
│   └── Genome_heatmap_cluster_1Mb_bins_edivisive.pdf
├── cellenone/             only if CellenONE is configured
├── qc_review/
│   ├── review.html        review report (open this)
│   └── plots/             per-well profile and histogram images
├── qc_decisions.csv       written by you after review
├── aneufinder_reviewed/   second pass, PASS wells only
├── CN_review/
│   ├── cn_review.html     final copy-number viewer (the deliverable)
│   ├── plots/
│   └── genome_heatmap.png
└── logs/
```

### 384-well plate

Preprocessing outputs (`demux/` to `dedup/`) stay inside each subplate directory.
AneuFinder, review, CellenONE and the final viewer are plate-level:

```text
plate21/
├── plate21_1/ .. plate21_4/   demux/ filtered/ fastqc/ multiqc/ raw_bam/ bam/ dedup/
├── aneufinder/                one pass over all 384 wells
├── cellenone/
├── qc_review/
│   ├── review.html
│   ├── plots/
│   └── assets/                cell images, loaded on click
├── qc_decisions.csv
├── aneufinder_reviewed/
└── CN_review/
```

A 96-well `review.html` is one self-contained file (about 20 MB) and opens over Samba
or as a local file. A 384-well report would be about 107 MB that way, so it loads
plots and images from `plots/` and `assets/` next to it instead. Copy the whole
`qc_review/` folder if you move it. Set `qc_review.embed_assets: true` to force a
single file.

Well identity in 384 mode is `<subplate>_<well>`, for example `plate21_2_W07`. This
is the same ID the subplate gets in a standalone 96-well run.

## Running modes

Set `MODE` at submit time. Defaults in `submit_pipeline.sh` can also be edited.

| Mode | What it runs |
|------|--------------|
| `pre_review` | Everything up to the first-pass review report, then stops. Default. |
| `post_review` | Validate decisions, second AneuFinder pass, final viewer. Needs `qc_decisions.csv`. |
| `preprocessing` | Demux, filtering, alignment, dedup |
| `qc` | Preprocessing plus FastQC and MultiQC |
| `aneufinder_first` | Blacklist plus first AneuFinder pass |
| `review` | Rebuild `review.html` only |
| `validate_review` | Check `qc_decisions.csv` only |
| `aneufinder_reviewed` | Second AneuFinder pass only |
| `cn_review` | Rebuild `cn_review.html` only |
| `cellenone` | CellenONE ingest and images only |
| `blacklist` | Blacklist diagnostic plots only |

Other submit-time variables:

| Variable | Default | Purpose |
|----------|---------|---------|
| `PLATE_DIR` | set in script | Plate (or 384 plate) directory |
| `PLATE_FORMAT` | `96` | `96` or `384` |
| `SUBPLATES` | empty | 384 only: subplate directory names, when they do not follow `<plate>_<n>` |
| `SNAKEMAKE_EXTRA` | empty | Extra Snakemake flags, for example `--nolock` |

`submit_pipeline.sh` stops early with a clear message if `PLATE_FORMAT` does not
match the directory layout, or if a post-review mode runs before `qc_decisions.csv`
exists.

## Review workflow

### Automatic status (read count only)

Each well gets an automatic status from `usable_reads`, the number of mapped reads
in the deduplicated BAM. It pre-fills the decision. The reviewer can change any
decision.

| usable_reads | Status | Default decision |
|--------------|--------|------------------|
| 100,000 or more | `PASS` | PASS |
| 50,000 to 99,999 | `WARN` | REVIEW |
| below 50,000 | `FAIL` | EXCLUDE |
| missing | `UNKNOWN` | REVIEW |

Cutoffs are `qc_review.usable_reads_pass_cutoff` and `usable_reads_warn_cutoff`.
Duplication rate is shown for information only and does not affect the status.

### Making decisions

1. Open `<PLATE_DIR>/qc_review/review.html`.
2. Click a well to see its copy-number profile, read-count histogram, metrics and
   (if configured) cell images.
3. Choose PASS, EXCLUDE, REVIEW or REPEAT, add reasons and notes, and click
   **Save decision**.
4. Click **Copy terminal save command** and paste it into a cluster terminal. This
   writes `<PLATE_DIR>/qc_decisions.csv`. **Download qc_decisions.csv** and placing
   the file at that path also works.
5. Run `MODE=post_review`.

Only PASS wells go to the second pass. Set `qc_review.include_review: true` to also
include REVIEW wells. The CSV schema is in [`config/README.md`](config/README.md).

The final `CN_review/cn_review.html` shows the whole plate. Wells kept for the second
pass show their final profile; excluded wells are greyed out.

## 384-well plates

A 384-well plate is dispensed as four interleaved 96-well subplates. Subplate and
well are derived from the 384 position by formula. For row `r` (0 to 15, A to P)
and column `c` (0 to 23, 1 to 24):

```text
subplate_index = (r % 2) + 2 * (c // 12)      # 0 is SL1, 3 is SL4
well_number    = (c % 12) * 8 + (r // 2) + 1  # W01..W96
```

So SL1 is rows A, C, E... in columns 1 to 12, SL2 is rows B, D, F... in columns
1 to 12, and SL3 and SL4 are the same over columns 13 to 24. Check this against a
plate map with:

```bash
python workflow/scripts/reporting/plate384_layout.py --check-xlsx plate_mapping_384.xlsx
```

For a non-standard layout, set `plate_layout_tsv` in `config.yaml` (columns
`subplate`, `well`, `pos384`).

**Reusing subplate runs.** AneuFinder models are per cell, so a subplate already run
in 96 mode has the same models as a plate-level run. `aneufinder.plate384_models`
controls reuse:

| Value | Behaviour |
|-------|-----------|
| `auto` | Reuse if every subplate has models, otherwise compute (default) |
| `reuse` | Always reuse; error if a subplate has no models |
| `rerun` | Always recompute at plate level |

With reuse, a full 384-well review is ready in minutes.

The review report shows one 16 x 24 plate map with tabs `All | SL1 | SL2 | SL3 | SL4`.
The exported CSV has an extra `subplate` column.

## Configuration reference

All settings are in `config.yaml`. The ones most often changed:

| Key | Default | Notes |
|-----|---------|-------|
| `plate_format` | `96` | Usually set with `PLATE_FORMAT` at submit time |
| `subplates` | `[]` | Explicit subplate list for 384 mode; empty means auto-discover |
| `aneufinder.binsize` | `1000000` | Must match the GC template |
| `aneufinder.min_reads_for_model` | `100` | Wells below this are not given to AneuFinder. A near-empty well otherwise aborts the whole batch. |
| `aneufinder.plate384_models` | `auto` | See [384-well plates](#384-well-plates) |
| `qc_review.usable_reads_pass_cutoff` | `100000` | |
| `qc_review.usable_reads_warn_cutoff` | `50000` | |
| `qc_review.include_review` | `false` | Include REVIEW wells in the second pass |
| `qc_review.plot_res` / `plot_res_384` | `150` / `100` | Plot DPI |
| `qc_review.embed_assets` | `auto` | Single-file report at 96, sidecar assets at 384 |
| `cellenone.runs` | | Plate name to `.Run` folder |
| `cellenone.image.channels` | all | Channels to render; missing or empty channels are skipped |
| `resources.<rule>` | | SLURM threads, memory, time, partition per rule |

## Environments

- Snakemake itself runs from the `snakemake_scDNA` conda environment set in
  `submit_pipeline.sh`.
- Each rule uses a conda environment from `workflow/envs/`. These are built once and
  shared under `.snakemake/conda` at the repo root, so every plate reuses them.
- R analysis code uses `renv`. Restore it on a fresh checkout with
  `R -e 'renv::restore()'`.

## Repository layout

```text
scDNA_pipeline/
├── Snakefile
├── config.yaml
├── multiqc_config.yaml
├── submit_pipeline.sh       SLURM entry point
├── pre_check.sh             pre-flight checks
├── config/                  qc_decisions.csv schema and templates
├── docs/                    CellenONE and QC notes, analysis write-ups
├── resources/reference/     GC template and mappability resources
├── workflow/
│   ├── envs/                per-rule conda environments
│   └── scripts/
│       ├── preprocessing/
│       ├── analysis/
│       ├── reporting/
│       └── dev/
├── ad_hoc_checks/           one-off analyses, not part of the pipeline
├── renv/, renv.lock
└── Seqinfo/
```

## Further documentation

- [`config/README.md`](config/README.md): `qc_decisions.csv` schema and decision rules
- [`docs/CELLENONE_AND_QC.md`](docs/CELLENONE_AND_QC.md): how to read the review report, and what CellenONE measures versus what the pipeline derives
- [`FUTURE_PLANS.md`](FUTURE_PLANS.md): planned work
