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

## Plate formats: 96 and 384

`PLATE_FORMAT` (default `96`) selects how a plate directory is interpreted.

**96-well** — one plate directory = one 96-well plate = one FASTQ pair. This is the
original behaviour and is completely unchanged.

**384-well** — a 384-well plate is dispensed as four interleaved 96-well subplates,
each sequenced as its own FASTQ pair:

```text
384_well/plate21/
├── plate21_1/          # SL1  — preprocessed exactly like a standalone 96-well run
├── plate21_2/          # SL2
├── plate21_3/          # SL3
├── plate21_4/          # SL4
├── aneufinder/         # ONE plate-level pass over all 384 wells
├── cellenone/          # CellenONE cell images (optional)
├── qc_review/          # ONE 384-well review
└── CN_review/          # ONE 384-well CN viewer
```

```bash
PLATE_FORMAT=384 PLATE_DIR=/…/384_well/plate21 sbatch submit_pipeline.sh
```

Per-subplate preprocessing (demux → dedup → MultiQC) stays under each subplate
directory. Only AneuFinder and the reporting layer become plate-level.

Subplates are auto-discovered as `<plate>_<n>` directories holding a FASTQ pair. When
they do **not** carry the plate directory's name — `384_well/plate17_umi/` holds
`plate17_1`..`plate17_4` — name them explicitly instead:

```bash
PLATE_FORMAT=384 PLATE_DIR=/…/384_well/plate17_umi \
SUBPLATES="plate17_1 plate17_2 plate17_3 plate17_4" sbatch submit_pipeline.sh
```

The 384 position of a well is derived, not looked up. For 384 row `r` (0–15 = A–P)
and column `c` (0–23 = 1–24):

```text
subplate_index = (r % 2) + 2 * (c // 12)      # 0 -> SL1 … 3 -> SL4
well_number    = (c % 12) * 8 + (r // 2) + 1  # W01..W96
```

so SL1 = rows A,C,E,… × cols 1–12, SL2 = rows B,D,F,… × cols 1–12, SL3/SL4 the same
over cols 13–24, filling column-major within each subplate. Verified against all 384
rows of `plate_mapping_384.xlsx`:

```bash
python workflow/scripts/reporting/plate384_layout.py --check-xlsx plate_mapping_384.xlsx
# 384/384 OK
```

Set `plate_layout_tsv` in `config.yaml` only for a non-standard layout.

Well identity in 384 mode is `<subplate>_<well>` (e.g. `plate21_2_W07`) — the same
`sample_id` a standalone run of that subplate produces, so a cell keeps its identity
whether it is reviewed alone or as part of the plate. That one string is also the
staged BAM basename, the AneuFinder model id, the `.RData` name and the plot name,
which is why `run_aneufinder.R` needs no changes for 384 mode.

### Reusing already-processed subplates

If a plate's four subplates have already been through their own 96-well runs, the
plate-level AneuFinder pass costs nothing. AneuFinder is strictly per-cell — each BAM
is independently binned, GC-corrected and segmented — so a subplate's model is
identical to the plate-level one. `aneufinder.plate384_models` controls this:

| value   | behaviour                                                       |
|---------|-----------------------------------------------------------------|
| `auto`  | reuse when every subplate already has models, else compute (default) |
| `reuse` | always reuse; error if a subplate has no models                 |
| `rerun` | always recompute at plate level                                 |

In `reuse` mode the existing models are symlinked into the plate-level `MODELS/`
directory under their namespaced ids, and a full 384-well review becomes available in
minutes without running AneuFinder once.

## Inputs

Each (sub)plate run expects:

- paired FASTQs named after their directory: `<dir>/<dir>_R1.fastq.gz` and `_R2.fastq.gz`
- a barcode table for well demultiplexing
- reference configuration in `config.yaml`
- a mappability reference BAM for blacklist diagnostics
- a GC template RDS for the AneuFinder stage

Barcode lookup follows this order:

1. `resources/barcodes.tsv`
2. `resources/barcodes/barcodes.tsv`
3. `<plate_dir>/barcodes.tsv` — **resolved per subplate in 384 mode**, i.e.
   `<plate_dir>/<subplate>/barcodes.tsv`

This allows a shared barcode definition to be reused across multiple plate runs
without copying `barcodes.tsv` into every plate directory. 384 mode additionally
requires the `well_id` order to be identical across subplates, and says so explicitly
if it is not.

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
<PLATE_DIR>/aneufinder_reviewed/      # second-pass models (PASS wells)
<PLATE_DIR>/CN_review/                 # FINAL deliverable: cn_review.html + plots/ + genome_heatmap.png
```

The `CN_review/` folder is the clearly-named final step. The `cn_review.html` viewer
shows the **whole plate**: wells kept for the second pass are coloured by their original
automated QC status and carry their final CN profile + histogram; wells excluded at
review are greyed out. The final `cn_review.html` renders per-well plots as **SVG**
(crisp/responsive vector); the larger 96-well first-pass `review.html` uses **PNG** (kept
compact). Both are shown full-width so the genome x-axis is readable. So the normal
operator sequence is just two runs:

1. `MODE=pre_review` → review HTML, then stop
2. save `<PLATE_DIR>/qc_decisions.csv` (use the report's **Copy terminal save command**)
3. `MODE=post_review` → final `CN_review/cn_review.html`

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

## CellenONE cell images (optional, display-only)

For every well the CellenONE dispenser photographs the drop in three channels
(Transmission, Blue, Orange) before ejecting it. The pipeline overlays those,
artificially colours the fluorescence channels, and scores the transmission image
against the two ejection lines to produce an **image-based call**:

| condition                                          | call            |
|----------------------------------------------------|-----------------|
| exactly 1 object                                   | `SINGLE`        |
| ≥2 objects, rightmost right of the purple line     | `PASS`          |
| rightmost between the green and purple lines       | `CONTAMINATION` |
| all objects left of the green line                 | `FAIL`          |
| 0 objects                                          | `NO_OBJECT`     |

The nozzle is at the **left** and cells sediment leftwards toward it, so an object far
to the right is still up the capillary and will not be ejected.

**This call is display-only.** Read-count gating remains the sole driver of
`auto_status` and the default decision; the image call gets its own panel and offers a
*suggested* reason chip that the reviewer must click to apply. It introduces no new
reason tokens.

Enable it by mapping the plate to its run folder in `config.yaml` — this mapping is
not derivable from any name (CellenONE's `plate_2` is our `plate21`) and is never
guessed:

```yaml
cellenone:
  runs:
    plate21: "/…/P21_22/plate_2/K1563_plate_2_20260707_135300_812.Run"
```

Detection parameters (`DetDiaMinTrans`, the isolation window, fluorescence intensity
limits) are always read from that run's own tables, never hardcoded — they vary
between runs. Iterate on just this layer with `MODE=cellenone`.

### Plates with no cell images

Nothing needs switching off. A plate absent from `cellenone.runs` (plate17, plate19 and
every plate dispensed before the CellenONE) never builds a `cellenone/` directory, the
review rule drops `--cellenone-dir`, and `generate_qc_review.py` marks the report
`cell_images: false`. The report then omits the image panel and the droplet call on the
hover line entirely — rather than showing 384 empty frames — so the per-well evidence is
the **copy-number profile and the bin-count histogram**. Everything else (the read gate,
the diagnostics grid, decisions, the CSV) is unchanged.

The same thing happens if a run folder is configured but yields no usable image: the
report logs a warning and falls back to the no-image layout instead of failing.

`qc_review.cell_images: false` in `config.yaml` forces that layout even for a plate that
*does* have images.

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
- `CN_review/`
- `cellenone/` (only when CellenONE images are configured)
- `logs/`

In 384 mode the first six of those stay under each `<plate>/<subplate>/` directory and
the rest are plate-level.

These generated outputs are intentionally excluded from version control.

### Report size

A 96-well `review.html` is a single self-contained file with every plot base64-embedded
(~20 MB), which is what makes it open correctly over Samba or as a `file://` document.
At 384 wells that approach would produce a ~107 MB page, so 384 mode switches to
relative asset references (`qc_review.embed_assets: auto`) and the browser fetches only
the clicked well. Measured on plate21:

| | 96-well (plate21_1) | 384-well (plate21) |
|---|---|---|
| `review.html` | 20.4 MB (self-contained) | **0.46 MB** + referenced assets |
| `qc_review/plots/` | 16 MB @ 150 dpi | 46 MB @ 100 dpi (`plot_res_384`) |
| `qc_review/assets/` (cell images) | — | 13 MB (768 JPEGs, ~17 KB each) |

Both size levers are config: `qc_review.plot_res_384` (the per-well plot DPI) and
`cellenone.image.channels` — the default `[merge, trans]` is two JPEGs per well;
adding `blue` and `orange` doubles the assets directory.

`qc_review.embed_assets: true` forces the single-file behaviour at any size.

### Reviewing a 384 plate

`review.html` shows one 16×24 plate map (rows A–P, columns 1–24) with a tab strip
`All | SL1 | SL2 | SL3 | SL4` above it. Selecting a subplate dims the other 288 wells;
**dimmed wells stay clickable**, since the filter is a focus aid rather than a lockout.
At 24 columns there is no room for per-cell text, so each cell carries a tooltip and a
hover readout appears under the map. The exported CSV gains a `subplate` column.

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
