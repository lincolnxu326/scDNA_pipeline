# scDNA Pipeline Agent Context

## Branch `welldr`: read this first

This branch targets the wellDR-seq protocol (Wang et al., Cell 188, 6355-6369, 2025)
and is the default development branch. Legacy behaviour (NlaIII library, in-pipeline
demultiplexing, UMI-tools deduplication, 384 subplate mode) lives on branch `dlp+`.

### Before doing anything else

1. Read, in this order:
   - `README.md` (the target pipeline and its output contracts);
   - `docs/WELLDR_PLAN.md` (the source of truth for the migration: protocol facts,
     decisions, open items, change inventory, phases);
   - `docs/WELLDR_BRANCH_SETUP_PROMPT.md` (the original brief, for the reasoning
     behind the plan; it has already been carried out, do not run it again).
2. Work out where the migration stands: which phases in `docs/WELLDR_PLAN.md` are
   done (check `git log` and the code), and which open items in section 3.2 are still
   open.
3. Before writing any code, explain the plan to the user in a few sentences: the
   goal, the phase the branch is in, what the next phase would change, and the open
   items that block it. Then ask what they want to do. This applies whether the user
   wants to develop the pipeline or only to understand it.
4. Work one phase at a time and stop for user review at the end of each phase, as
   the plan says.

The rest of this file describes the legacy pipeline as it stood when the branch was
created. Sections on legacy-only topics apply to `dlp+` only. The AneuFinder / GC /
blacklist warnings still apply here.

## Purpose

This document is an internal handoff note for future coding agents working in
`VCAM1_GnT/TOOLS/scDNA_pipeline`.

It is intentionally detailed and operational. It captures:

- what this project is
- how the current pipeline is structured
- what has already been changed
- what is currently fragile
- what not to break
- which future plans have been drafted but intentionally parked

This file is not meant to be polished end-user documentation. The user-facing
documentation remains in `README.md` and `FUTURE_PLANS.md`.


## Project Overview

This repository contains a Snakemake-based single-cell DNA / Strand-seq style
pipeline used for per-plate processing of single-cell sequencing data.

The main workflow currently covers:

1. Preprocessing
   - demultiplexing FASTQs by well
   - deduplication
   - adapter / dimer filtering
   - alignment to hg38

2. QC
   - FastQC
   - MultiQC
   - blacklist generation for AneuFinder support
   - a currently semi-manual GC template workflow for AneuFinder

3. CNV calling
   - AneuFinder-based downstream analysis on per-well BAMs

The repo also contains ad hoc R / Quarto scripts for manual QC exploration,
Watson/Crick diagnosis, GC generation, and plotting work that is not part of
the strict production workflow.


## Current Top-Level Structure

- `Snakefile`
  - main Snakemake workflow
- `config.yaml`
  - workflow configuration
- `submit_pipeline.sh`
  - main launch script for cluster execution
- `pre_check.sh`
  - user-facing helper
- `workflow/scripts/`
  - pipeline-owned execution scripts
- `workflow/envs/`
  - conda env YAMLs used by Snakemake
- `ad_hoc_checks/`
  - manual / exploratory scripts not owned by the production workflow
- `resources/reference/`
  - shared reference outputs used across plates
- `.snakemake/conda`
  - shared env prefix used across runs
- `.conda/pkgs`
  - shared conda package cache used across runs


## Intended Logical Structure of the Workflow

The user explicitly wants the overall architecture understood as:

1. `Preprocessing`
   - demux
   - dedup
   - dimer removal
   - alignment

2. `QC`
   - FastQC
   - MultiQC
   - Ashley QC in the future
   - Strand-seq plots in the future
   - current support resources for AneuFinder such as blacklist

3. `CNV`
   - AneuFinder

This matters because future features should fit into that structure rather than
being added in an ad hoc order.


## Current Mode Semantics

The desired mode split is:

- `preprocessing`
  - preprocessing only
- `qc`
  - preprocessing + FastQC + MultiQC
- `aneufinder`
  - blacklist / GC precheck / AneuFinder only
- `full`
  - QC + AneuFinder end-to-end

Important:

- `aneufinder` mode should not force a blacklist rebuild if the shared blacklist
  already exists and is up to date.
- `aneufinder` mode should not route through the default `all` target if the
  intent is specifically AneuFinder-only behavior.
- manual GC generation remains a prerequisite for AneuFinder right now because
  of the package compatibility problem described below.


## Important Historical Changes Already Made

### 1. Shared conda env prefix

Originally, environment creation was happening in per-plate or per-run
locations. This was reworked so Snakemake environments are shared under the
pipeline root rather than recreated in sample folders.

Implication:

- same YAML should resolve to the same env across plates
- env reuse now depends on YAML content + shared conda prefix, not on plate path

### 2. User-facing shell scripts at the pipeline root

The user wants user-run helper scripts easy to find at repo root.

So these are intentionally top-level and should stay easy to discover:

- `submit_pipeline.sh`
- `pre_check.sh`

Do not bury these under a helper subdirectory unless the user explicitly asks.

### 3. `workflow/scripts` versus ad hoc scripts

Pipeline-owned execution code belongs in `workflow/scripts/...`.

Manual exploratory analysis or development-only scripts belong in
`ad_hoc_checks/`.

The user explicitly said legacy or exploratory scripts should be moved out of
the workflow tree if they are not part of production execution.

### 4. Git repo setup

This directory was turned into a Git repo and a repo-facing `README.md` was
written. Runtime outputs, rule graph outputs, and generated artifacts are
ignored.

### 5. MultiQC adjustments

Work was already done to reduce sample-name fragmentation in MultiQC and to
improve general statistics visibility, including mapped read counts.

### 6. Watson/Crick plotting work

There was iterative work on Watson/Crick / Strand-seq style plotting in
`ad_hoc_checks`, including a user request to orient chromosomes vertically and
show Watson/Crick signal on either side of the chromosome axis. Some plotting
work was moved into `QC_plots_dev.qmd`, but this area should be treated as
still exploratory and user-driven.


## Current AneuFinder / GC / Blacklist State

This is the most fragile and important current area.

### What the user wanted

Originally, the idea was to automate:

- shared mappability reference support
- shared blacklist generation
- shared GC template generation

under the pipeline root so all plates could reuse them.

### What actually happened

Blacklist generation is now shared under:

- `resources/reference/mappability/`

However, GC generation had to be backed out of the workflow.

### Why GC generation is manual right now

The user explicitly clarified that:

- the GC correlation step had been kept outside the workflow because
  `AneuFinder` and `BSgenome.Hsapiens.UCSC.hg38` conflict in the current
  project setup
- `hg38` was not reliably available in the `renv` path used by the pipeline
  execution
- therefore GC generation has to remain a manual step for now

This means:

- Snakemake may generate the blacklist
- the GC template is then produced manually from `ad_hoc_checks/generate_gc_corr.R`
- AneuFinder then consumes the manually produced shared GC RDS

### Shared resource locations currently in use

Shared blacklist:

- `resources/reference/mappability/blacklist.bed.gz`

Shared GC template:

- `resources/reference/hg38_binsize1000000_variable_bins_with_GC.rds`

### Important workflow behavior currently intended

`check_gc_rds` should verify the manual GC template exists without owning or
deleting it.

There was a bug where `check_gc_rds` declared the manual `.rds` as its output,
which caused Snakemake to treat that file as job-owned and remove or refresh it
on rerun. That behavior was explicitly identified as wrong and corrected.

Rule design principle here:

- manual artifacts must not be declared as Snakemake-owned outputs if the user
  is expected to create them outside the workflow

### Current debugging history for AneuFinder

This sequence already happened:

1. GC template initially missing from workflow path
2. GC precheck and blacklist orchestration adjusted
3. GC file accidentally treated as rule output and removed by Snakemake
4. corrected so a separate flag is used for readiness rather than owning the
   manual file
5. AneuFinder then began failing inside `correctGC()`
6. multiple hypotheses were tested:
   - wrong object class (`list` vs `GRangesList`)
   - wrong metadata type (`matrix` vs numeric column)
   - wrong metadata column name (`gc` vs `GC`)
   - wrong bin naming (`binsize_1e+06` vs `binsize_1000000` vs `0`)
   - mismatched bin structure between sample bins and GC template
7. the original legacy GC script logic was restored more closely because it was
   more likely to match the format AneuFinder expected

### Critical lessons from that debugging

- AneuFinder is extremely sensitive to the exact structure of the GC template
- class, naming, metadata field names, and bin construction path all matter
- the original legacy script uses `GC` uppercase, not lowercase
- the original legacy script may produce a `GRangesList` named `"0"`
- sample bins in `run_aneufinder.R` may surface names like:
  - `binsize_1000000`
  - `binsize_1e+06`
- `run_aneufinder.R` had to be patched multiple times to normalize or preserve
  these names

### Current best mental model

The remaining AneuFinder / GC area should be treated as unstable until the user
confirms the end-to-end run is clean.

Do not casually simplify:

- object class conversions
- `correctGC()` inputs
- GC metadata naming
- bin naming normalization

Any future edits in this area should log and compare:

- `class(bins_gc_)`
- `names(bins_gc_)`
- `class(binned_obj_)`
- `names(binned_obj_)`
- range identity between precomputed bins and GC bins


## Current `run_aneufinder.R` Status

This file has been modified during debugging to:

- source `renv/activate.R`
- avoid suppressing package startup messages
- normalize bin naming conventions such as scientific notation vs integer form
- preserve bin naming better than the earlier `GRangesList("0" = ...)` style
- deduplicate normalized bin names

This file is now an accumulation of debugging fixes. Future agents should read
it carefully before changing it again.

Do not assume it still matches upstream or legacy AneuFinder examples.


## Current `generate_gc_corr.R` Status

The ad hoc GC generation script was intentionally reverted closer to the
original working logic after newer rewrites appeared to drift away from the
format AneuFinder expected.

The user explicitly asked to:

- keep the absolute path configuration style
- keep manual top-of-file parameter settings
- not convert this into a fully CLI-driven script again

So this file currently uses hardcoded parameter values near the top and should
remain that way unless the user asks otherwise.


## Submission / SLURM Notes

### Historical job-submission issue

There was a critical bug where `submit_pipeline.sh` resolved the pipeline root
relative to the SLURM spool copy under `/tmp/slurmd/...`, causing shared state
such as `.snakemake` and `.conda` to be created in the wrong place and fail
with permissions errors.

This was fixed by basing the pipeline root on `SLURM_SUBMIT_DIR` when present.

### Desired future submission UX

The user wants future submission behavior improved significantly:

- one variable should be enough to select a plate
- launcher job names should include the plate
- SLURM logs should be named deterministically
- rule-level child jobs should have plate-aware names
- the user should not need to edit many paths per run

This is still a future-plan area, not a fully completed one.


## MultiQC Notes

The user wanted:

- proper merging of same-well metrics into one sample row
- more complete general statistics columns
- mapped read counts included

There was already work done to fix sample name normalization and to feed
samtools stats / flagstat data to MultiQC more clearly.

If future agents touch MultiQC again, they should preserve the goal that each
well appears once in the summary table.


## QC Plotting / Watson-Crick Notes

The user has a strong opinion about what a Strand-seq plot should look like.

It should not be a generic chromosome-on-x plot.

The requested visual style is:

- chromosomes arranged vertically like karyotype columns
- genomic position running vertically within each chromosome
- Watson and Crick signal shown on opposite sides of a chromosome axis

This requirement came up explicitly after an earlier attempt that was close but
not the desired layout.

There is also an active Quarto development notebook:

- `ad_hoc_checks/QC_plots_dev.qmd`

This is being used to explore:

- duplication asymmetry between R1 and R2
- read-depth unevenness across wells
- Watson/Crick imbalance
- per-plate and cross-plate QC summaries

This area is exploratory and user-driven, not yet frozen.


## QC Review Sub-stage (static HTML)

A static HTML human-review layer was added after MultiQC + AneuFinder.

What it is:

- `Snakefile` rules `render_well_profiles`, `qc_review`, `validate_qc_decisions`
  (+ alias target `qc_review_report`, and `review.html` added to `rule all`).
- `workflow/scripts/reporting/render_well_profiles.R` — renders one
  `qc_review/plots/{well}.png` per well from the existing AneuFinder
  `MODELS/method-{method}/{well}.RData` models, plus `manifest.json`.
- `workflow/scripts/reporting/generate_qc_review.py` — aggregates per-well QC
  metrics, computes an advisory PASS/WARN/FAIL from `config.yaml qc_thresholds`,
  base64-embeds the PNGs, and writes the single self-contained
  `qc_review/review.html`.
- `workflow/scripts/reporting/validate_qc_decisions.py` — validates a
  human-saved `<PLATE_DIR>/qc_decisions.csv` (per-plate, in the data dir; run on
  demand; not in `all`).
- `config/qc_decisions.csv.example`, `config/README.md`, `submit_pipeline.sh`
  `MODE=review`, `config.yaml` `qc_review:` block + `resources.qc_review`.

Key principles to preserve here:

- **The renderer is a read-only consumer of `run_aneufinder.R` output.** It must
  not modify, re-run, or take ownership of any AneuFinder file. Do not couple it
  to the fragile AneuFinder/GC code — it only `load()`s `.RData` and re-plots.
- `render_well_profiles.R` runs via `renv` (no `conda:`), exactly like
  `run_aneufinder.R` / `generate_blacklist.R`, because AneuFinder lives in renv.
- The review HTML is fully static: no server, no DB. The browser only offers a
  copy/download of `qc_decisions.csv`; it never writes into the project tree.
- **The single-file property now holds only for 96-well runs.** A 384 report would be
  ~107 MB fully embedded, so it references `plots/…` and `assets/…` relatively
  instead (`qc_review.embed_assets: auto`). See the 384 section below.
- Well IDs here are sequential `W01..W96` (not positional `A01..H12`). The plate
  map uses a sequential 8x12/grid fill; positional layout is auto-detected only
  if every well matches `A-H` + `01-12`. Verify against the barcode table if a
  truly positional plate ever appears.
- Well IDs here are sequential `W01..W96` (not positional `A01..H12`).


## Two-Pass AneuFinder (post-review second pass)

The QC review now feeds a **second AneuFinder pass**. Flow:

```
preprocessing -> run_aneufinder (FIRST pass, aneufinder/) -> qc_review/review.html
  -> [human saves <PLATE_DIR>/qc_decisions.csv] -> validate_qc_decisions
  -> derive_included_wells (qc_review/included_wells.tsv, PASS [+REVIEW if cfg])
  -> run_aneufinder_reviewed (SECOND pass, aneufinder_reviewed/) -> render_reviewed_profiles
  -> cn_review (CN_review/cn_review.html)   # final deliverable lives in CN_review/
```

Final-output layout: `qc_review/` = first-pass review (review.html, plots/, included_wells.tsv,
validated flag); `aneufinder_reviewed/` = second-pass models; **`CN_review/`** = final deliverable
(`cn_review.html`, `plots/`, `genome_heatmap.png`). Driven by `CN_REVIEW_DIR` in the Snakefile.

Critical design rules to preserve:

- **`run_aneufinder.R` is reused UNCHANGED for the second pass.** `run_aneufinder_reviewed`
  builds a symlink-only input dir (`aneufinder_reviewed/input_bams/{well}.bam ->
  ../../bam/{well}.bam`) of just the PASS wells, then calls the same `run_aneufinder.R`
  with `--input <link_dir> --output aneufinder_reviewed`. Do NOT add a well-subset
  option to the fragile R script — keep the subsetting in the Snakemake shell.
- **First/second passes are separate dirs:** `aneufinder/` (first, untouched by the
  second pass) vs `aneufinder_reviewed/` (second). Original BAMs are never copied or
  mutated — only symlinked.
- **Plots are rendered in R from `.RData`, never by parsing `profiles_*.pdf` page order.**
  `render_well_profiles.R` re-plots per-well models (id == well, reliable) and, with
  `--heatmap`, also renders the genome heatmap via `heatmapGenomewide()`. The conda QC
  env has no PDF rasterizer, and PDF page order is treated as unsafe.
- **The decisions CSV only gates the post-review targets.** `derive_included_wells` and
  everything after it depend on `qc_decisions.validated.flag`; the pre-review targets
  (`all`, `all_qc`, `qc_review_report`) never reference the CSV, so a missing CSV cannot
  block a normal run.
- **`MODE=full` was renamed to `pre_review`** (alias `full` kept) because it stops at the
  first-pass review report; the second pass is always an explicit, opt-in mode.
- **`MODE=post_review` is the normal post-review entry point** — a thin target
  (`rule post_review`, input `CN_review/cn_review.html`) that pulls the whole post-review
  chain in one run. The lower-level modes (`validate_review`, `aneufinder_reviewed`,
  `cn_review`) are kept for debugging/partial reruns. `submit_pipeline.sh` guards the
  post-review modes with an early check that `<PLATE_DIR>/qc_decisions.csv` exists.
- Inclusion semantics: `PASS` always; `REVIEW` only if `qc_review.include_review`;
  `EXCLUDE`/`REPEAT` excluded (`REPEAT` logged). `cn_review.html` is a read-only viewer
  (no decision controls / no CSV export) built by `generate_qc_review.py --report-kind cn`.
- **cn-viewer plate map** shows the WHOLE plate: wells kept for the second pass are coloured
  by their ORIGINAL auto-status (PASS/WARN/FAIL/UNKNOWN); excluded wells are greyed out
  (`included` flag per well). The legend is kind-aware (status colours in cn-mode, decision
  colours in review-mode). Do NOT restrict the cn map to only included wells.
- **Per-well plot format is hybrid** via `render_well_profiles.R --format` (svglite→
  `grDevices::svg()`→PNG fallback; wide aspect): first-pass `render_well_profiles` rule uses
  `--format png` (96 wells → keeps `review.html` ~tens of MB; SVG made it ~150 MB),
  `render_reviewed_profiles` (cn) uses `--format svg` (few PASS wells → crisp/responsive).
  The generator embeds `{well}_{profile,histogram}.svg` then `.png` (so it handles either).
  Genome heatmap stays PNG. Profile is shown full-width above the histogram.

### QC gating + plots + UI (read-count-only)

- **Gating is read-count-only.** `compute_status()` in `generate_qc_review.py` uses
  `usable_reads = metrics["mapped_reads"]` (mapped reads in the dedup BAM, from
  `bam/{well}.stats.txt`) against `qc_review.usable_reads_pass_cutoff` (100000) /
  `usable_reads_warn_cutoff` (50000). Mapping → auto_status → default decision:
  PASS→PASS, WARN→REVIEW, FAIL→EXCLUDE, UNKNOWN(missing reads)→REVIEW(+reason
  `missing_qc_metric`). **Do not** reintroduce duplication/dimer/mapping-rate gating.
- **Duplication is never gated.** The UI shows it only as a labelled diagnostic:
  "UMI-tools duplicate rate" (= `dedup_rate` from `dedup/dedup_summary.tsv`, i.e.
  umi_tools `duplicate_reads/total`) + "UMI dedup retention". Keep the source explicit.
- **Two plots per well.** `render_well_profiles.R` writes `{well}_profile.png` and
  `{well}_histogram.png` (native AneuFinder `plot(model, type="profile"/"histogram")`),
  plus a back-compat `{well}.png` copy of the profile; manifest carries
  `profile_png/histogram_png/ok_profile/ok_histogram`. Same script, both passes.
- **Reasons are a controlled vocabulary, `;`-joined** in the single `reason` CSV column
  (`low_read_count, noisy_profile, poor_bin_distribution, low_complexity,
  suspected_doublet_or_mixed_well, sample_swap_suspected, manual_exception, other,
  missing_qc_metric`). `REASONS` is duplicated in `generate_qc_review.py` and
  `validate_qc_decisions.py` — keep them in sync; the validator splits on `;` and allows
  empty. CSV schema unchanged (`sample_id,well,decision,reason,notes`).
- **UI:** decision buttons + multi-select reason chips + notes + Save; cells colour by
  decision (auto-default until Saved, then dot-marked); a per-well badge shows
  auto-default vs saved. `cn` viewer shows profile+histogram read-only.
- **Decisions file is per-plate, in the data dir:** `<PLATE_DIR>/qc_decisions.csv`
  (NOT the repo `config/`). `DECISIONS_FILE` in the Snakefile resolves
  `qc_review.decisions_file` relative to `PLATE_DIR` (absolute paths honoured). The repo
  `config/` keeps only the template (`qc_decisions.csv.example`) + schema docs. The
  review report exposes a **"Copy terminal save command"** button: it emits a
  `mkdir -p '<dir>' && cat > '<PLATE_DIR>/qc_decisions.csv' <<'QC_DECISIONS_EOF' … EOF`
  heredoc (quoted delimiter, no shell expansion) built live from the in-page decisions,
  so the user pastes it into a cluster terminal to write the file at the right path.
  `generate_qc_review.py` receives the destination via `--decisions-path` (Snakefile
  passes `DECISIONS_FILE`).


## Parked Future Plans

These plans were drafted on purpose and then explicitly parked. Do not start
implementing them just because they exist.

### 1. Ashley QC integration

Planned direction:

- integrate `ashleys-qc` on existing per-well BAMs
- do not replace current preprocessing
- default to report-only predictions
- later allow optional downstream gating

Architecture target:

- `Preprocessing` -> `QC` -> `CNV`
- Ashley belongs in the QC phase

### 2. Watson/Crick Strand-seq plots as a formal module

Planned direction:

- separate from Ashley ML
- run on existing BAMs
- generate per-cell Strand-seq karyotype-style plots
- produce plate-level strand summaries

### 3. Submission UX cleanup

Planned direction:

- better SLURM job names
- better log names
- easier plate selection
- deterministic mapping between job id, plate, and log files

### 4. Performance optimization plan

Planned direction:

- improve demultiplex speed
- parallelize dedup / filter per well
- fix align CPU allocation vs actual scheduler CPU requests

Important insight already established:

- much of the preprocessing slowdown is because several scripts are single-core
  Python loops despite requesting many threads


## Do / Do Not Guidance

### Do

- preserve user-facing root scripts at repo top level
- keep exploratory scripts in `ad_hoc_checks/`
- preserve shared conda env prefix behavior
- keep shared blacklist under `resources/reference/mappability/`
- treat the GC template as a manually managed shared resource for now
- read logs carefully before assuming the failure is in Snakemake rather than
  the R layer
- keep file-path changes conservative and explicit
- use `apply_patch` for manual edits

### Do Not

- do not reintroduce per-plate conda env creation
- do not move root helper scripts into a hidden utility folder
- do not convert manual exploratory scripts into workflow-owned steps unless the
  user explicitly asks
- do not assume the current AneuFinder / GC debug area is stable
- do not make broad “cleanup” edits in `run_aneufinder.R` or
  `generate_gc_corr.R` without checking the exact AneuFinder object expectations
- do not treat manually created GC files as Snakemake-owned outputs
- do not revert user edits or unrelated worktree changes


## 384-Well Mode (`plate_format: 384`)

A 384-well plate is four interleaved 96-well subplates, each sequenced as its own
FASTQ pair. 384 mode keeps per-subplate preprocessing exactly as-is and adds ONE
plate-level AneuFinder + QC/CN review over all 384 wells.

### The two strings the whole 96/384 branch reduces to

Defined at the top of the `Snakefile`:

```python
SUB_BASE  # per-subplate output root: str(PLATE_DIR / "{sub}")  @384
          #                           str(PLATE_DIR)            @96  <- NO wildcard
LOG_TAG   # per-subplate log tag:     "{sub}" @384, PLATE @96
```

**In 96 mode `SUB_BASE` has no wildcard**, so every 96-mode path renders
*character-identical* to what it was before 384 mode existed: same on-disk outputs,
same `.snakemake/metadata` keys, no spurious re-runs. That property is the safety
mechanism — preserve it. `assert len(SUBPLATES) == 1` guards it.

`wildcard_constraints` for `sub` and `well` are **mandatory**: Snakemake wildcards
default to `.+`, which matches `/`, so an unconstrained `{sub}` swallows path
separators.

### `wid()` — one definition of well identity

```python
def wid(sub, well): return f"{sub}_{well}" if IS_384 else well
```

This single string is: staged BAM basename → AneuFinder model id → `.RData` basename
→ plot basename → HTML well key → `included_wells.tsv` row → `sample_id`. That chain
is *why* `run_aneufinder.R` and `generate_gc_corr.R` need **zero** edits for 384 mode
— only which files sit in the staged input dir and what they are named changes.

`sample_id` is `<subplate>_<well>`, identical to what a standalone 96-well run of that
subplate produces, so a cell keeps its identity either way.

### The staging script

`workflow/scripts/analysis/stage_aneufinder_input.py` replaced two inline shell loops
and serves both AneuFinder passes and both plate formats. It symlinks only — no BAM is
ever copied or mutated.

It exists because (a) AneuFinder aborts the **entire batch** if any single well has
~no reads, so dead wells must be filtered out, and (b) the staged basename is what
namespaces plate-level output.

The `run_aneufinder_reviewed` loop it replaced was
`while IFS=$'\t' read -r sample well`, which silently read the **subplate** column of
a 3-column `included_wells.tsv` as the well name. The script detects both schemas and
errors loudly when a 2-column file is ambiguous, rather than guessing.

### Model reuse — `aneufinder.plate384_models: auto|reuse|rerun`

AneuFinder is strictly per-cell: each BAM is independently binned, GC-corrected and
segmented, and no cross-well information enters a well's model. A model computed in a
subplate's own 96-well run is therefore identical to the plate-level one. Only the
cluster PDF and heatmap are batch-wide, and `render_well_profiles.R` regenerates those
at plate level from the `.RData` anyway.

`link_subplate_models.py` symlinks the subplates' models into the plate-level
`MODELS/method-*/` under their namespaced ids. Because the plot name comes from the
`.RData` basename, that alone yields correct plate-level plots at **zero** compute.
The mode is decided at parse time so `threads`/`resources` match what the job does.

### Rules by level

- **Per-subplate** (gain `{sub}` @384): `demultiplex`, `filter_dimers`, `fastqc`,
  `align`, `deduplicate_bam`, `summarize_dedup`, `multiqc`.
- **Plate-level** (paths unchanged): `run_aneufinder`, `render_well_profiles`,
  `qc_review`, `validate_qc_decisions`, `derive_included_wells`,
  `run_aneufinder_reviewed`, `render_reviewed_profiles`, `cn_review`, `ingest_cellenone`,
  `render_cellenone_images`.
- **Untouched**: `build_index`, `get_mappability_bam`, `generate_blacklist`,
  `check_gc_rds` (the last is dead in the DAG — its shell references a `{params.gc_rds}`
  that `params:` never defines, and nothing consumes its flag. **Leave it dead.**)

Rules that referenced the `{PLATE_DIR}` / `{LOG_DIR}` globals directly in `shell:` now
use `params.workdir` / `params.logdir`, because the globals resolve to the *parent* in
384 mode.

`rule align`'s `rgid_` gained the subplate. **`rgsm_` was deliberately NOT touched** —
changing SM would invalidate every existing BAM header and force realignment of all
768 BAMs.

### `rule clean` is deliberately non-recursive

It removes plate-level outputs only and never descends into `<plate>/<subplate>/`.
Recursing would delete 768 BAMs — days of alignment — for what reads like routine
cleanup.

### Reporting

`generate_qc_review.py` stays one template with one `__DATA_BLOB__` token; all new CLI
args are additive and defaulted, so with `--subplates` empty every existing path is
unchanged. The CSS hard numbers became custom properties **whose defaults are the
96-well values**, so a 96 report is pixel-identical; the boot code overrides them when
`grid.cols > 12`.

Two long-standing declared-but-unused arguments were wired up rather than left as
traps: `render_well_profiles.R --res` (`PNG_DPI` was hardcoded at 150) and the MultiQC
link, which was a hardcoded relative `../multiqc/multiqc_report.html` and is now a
`multiqc_links` list — the old form pointed at a plate-level directory that does not
exist for a 384 plate.

The genome heatmap height is capped at 14000 px: 46 px/row is fine for 96 wells but
384 would ask for 7200×17664 in one R device and OOM the job.

## CellenONE Cell Images (display-only)

`ingest_cellenone.py` (parse-only) → `render_cellenone_images.py` (pixels) →
`cellenone/` → an image panel in both the review and CN viewers.

Hard-won facts, all verified against the real run folders:

- The run-folder → plate mapping is **not derivable** (CellenONE's `plate_2` is our
  `plate21`). Config map only; never guess.
- `*.xls` files are **tab-separated text** with CRLF and trailing empty fields, not
  Excel. `ImageFile` is wrapped as `=HYPERLINK("…")`.
- Detection/isolation criteria **vary per run** (plate21 isolates 12–25 µm,
  plate22 9–25 µm). Always read them from that run's own file.
- `EjBound` in the xls is 330 and is **not** the green line. Green = 396 (the
  `Ejection (pix)` parameter), purple = 647 (= the `BackgroundEjZone` image width).
- Match images to wells by the `(POS)` token in the filename, **never** by count or
  index: there are 385 Blue/Orange PNGs but 384 table rows, and Red exists for only
  27 of 384 wells.
- `geoprops.xls` has exactly **one** Transmission row per printed drop, so CellenONE's
  own tables cannot answer "how many objects are in this well" — detection is required.
  geoprops is for run-level stats and cross-checking only.
- Calibrate the pixel scale by **nearest centroid** to the table's (X, Y), never by
  largest component: in 5 of 30 sampled wells the largest blob is not the isolated cell.
- The median over ~60 frames beats CellenONE's own `_Background.png`, which leaves
  static artifacts at x ≈ 668 and x ≈ 692 in every well.

**The image call is display-only and must stay that way.** It never enters
`compute_status()`; it lives in `d.cellenone` and offers a *suggested* reason chip the
human clicks. It introduced no new reason tokens — that vocabulary is frozen.

In the HTML, the channel-tab listeners are wired **above** `select()`'s
`if (APP.kind !== 'review') return;`, or the tabs would be silently dead in the CN
viewer.

## High-Risk Files

If a future agent edits any of these, extra care is needed:

- `Snakefile`
- `submit_pipeline.sh`
- `workflow/scripts/analysis/run_aneufinder.R`
- `ad_hoc_checks/generate_gc_corr.R`
- `multiqc_config.yaml`
- `ad_hoc_checks/QC_plots_dev.qmd`
- `ad_hoc_checks/Watson_Crick_diagnosis.R`

Before changing the `Snakefile`, always run the **96 regression gate**: a
`MODE=pre_review` dry run against an already-complete 96-well plate must still say
"Nothing to be done". Then run a 384 dry run and confirm the job list contains **no**
`align` / `deduplicate_bam` / `fastqc` / `multiqc` / `demultiplex` jobs. The shared
profile sets `rerun-triggers: mtime`, which is what makes reuse work — never pass
`-F` / `--forceall` / `--rerun-triggers code`, and never `touch` the FASTQs.


## Recommended First Checks for Future Agents

If a future agent is asked to debug a failed run, the first checks should be:

1. Determine whether failure is:
   - Snakemake DAG / mode selection
   - cluster / SLURM
   - conda / env
   - renv / R library
   - AneuFinder internals

2. Check:
   - top-level pipeline stderr under `logs/pipeline_*.err`
   - per-plate rule log
   - child SLURM rule log under `.snakemake/slurm_logs`

3. For AneuFinder failures specifically, log:
   - class and names of the loaded GC object
   - class and names of the sample binned object
   - whether bin ranges match exactly


## If the User Later Asks for a Commit

At various points the user asked for a commit, then redirected work before the
commit was made.

So if asked again:

- review the current working tree carefully
- do not assume the pending edits are all ready
- summarize what will be committed before actually committing if the worktree is
  mixed


## Final Note

This project has a relatively clear biological / workflow intent, but the
current AneuFinder / GC implementation details are still mid-debug. Future
agents should optimize for preserving working assumptions and producing precise
diagnostics over cleanup for elegance.
