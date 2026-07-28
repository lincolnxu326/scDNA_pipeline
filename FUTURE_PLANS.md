# Future Plans

This file tracks agreed future work that has been scoped but deliberately parked until the current pipeline is more stable.

## 1. Advanced QC Package

Bundle the next QC upgrade as one coordinated package:

- integrate Ashley QC on existing per-well BAM outputs
- add Watson/Crick Strand-seq plots
- keep both under one future advanced QC block
- keep predictions report-only by default at first
- avoid replacing the current FASTQ-to-BAM path during the first integration

Planned direction:

- use `ashleys-qc` directly on the BAMs already produced by this pipeline
- keep strand plots as a separate but coordinated post-alignment module
- write outputs into a plate-level QC area
- delay any downstream filtering or gating until the base workflow is stable

## 2. Submission UX And Logging Cleanup

Future submission improvements should make plate execution easier and logs easier to trace.

Planned direction:

- make one variable enough to select which plate to run
- derive the rest of the run paths automatically
- improve top-level `sbatch` job naming so plate and mode are visible
- improve rule-level submitted job names so they are not opaque random strings
- standardize scheduler `out` and `err` log paths and filenames
- make it easier to match SLURM job IDs to pipeline runs and per-rule logs

## 3. Performance Improvements And SLURM Resource-Allocation Audit

The current preprocessing path is slower than it should be, especially in demultiplexing and plate-level preprocessing. More generally, most stages across the whole pipeline feel slow to run, and the `resources:` block in `config.yaml` (per-rule `threads`/`mem_mb`/`time`/`partition`) was set generously rather than measured, so it's not clear how much of the observed runtime is real compute time versus time spent queued waiting for an oversized SLURM allocation.

Planned direction:

- optimize demultiplexing, which is currently dominated by single-process Python and gzip I/O
- review UMI-tools deduplication metrics after test plate runs
- split adapter filtering into per-well jobs for the same reason
- review alignment resource mapping so allocated SLURM CPUs match requested tool threads
- **audit actual vs. requested resources per rule**, using SLURM accounting (`sacct`/`seff`) or
  Snakemake `benchmark:` directives, before changing any numbers:
  - several per-well rules (`align`, `deduplicate_bam`) currently request the same large
    allocation (32 threads, 256–512 GB) as plate-level rules (`demultiplex`, AneuFinder passes);
    confirm whether each rule's tool actually uses that many threads/that much memory, or
    whether the request is oversized for what runs per well
  - check whether single-process scripts (e.g. the demux script) are requesting thread
    counts they cannot use, which wastes allocation without speeding up the job
  - check whether large requested allocations are causing queue delays on the `ncpu`
    partition that look like the pipeline "hanging" but are actually scheduling wait time
  - right-size `threads`/`mem_mb`/`time` per rule from the measured data, separately for
    plate-level vs. per-well rules, rather than reusing one large default across very
    different job sizes
  - only add tunable performance knobs once the current pipeline behavior is stable and
    the measured baseline exists to compare against
- add benchmarking and controlled performance knobs only after the current pipeline behavior is stable

## 4. Downstream Gating From QC Decisions (IMPLEMENTED)

Implemented as the **two-pass AneuFinder workflow**: the first pass (`aneufinder/`)
feeds `qc_review/review.html`; a validated `<PLATE_DIR>/qc_decisions.csv` is turned into
`qc_review/included_wells.tsv` and the second pass (`aneufinder_reviewed/`) reruns
AneuFinder on the PASS wells only, with a final `qc_review/cn_review.html` viewer.
First-pass outputs and original BAMs are never mutated. See `README.md` (Two-Pass
AneuFinder Workflow) and `AGENT_CONTEXT.md`.

Remaining parked ideas in this area:

- **REPEAT auto-rerun loop:** wells marked `REPEAT` are currently excluded and only
  logged; a future loop could re-sequence / reprocess them and re-review.
- **Multi-round review:** support more than two passes (review → reprocess → review)
  with versioned decision files.
- **Cohort assembly:** combine PASS wells across plates into a cohort-level CN matrix.

## 6. Higher-Throughput, Multi-Subplate Runs (IMPLEMENTED)

**Implemented as `plate_format: 384`.** Run it with:

```bash
PLATE_FORMAT=384 PLATE_DIR=/…/384_well/plate21 sbatch submit_pipeline.sh
```

What landed against the plan below:

- Per-subplate FASTQ→BAM processing kept exactly as-is, under
  `<plate>/<subplate>/`, with subplates tracked as members of one parent plate.
- One plate-level AneuFinder pass and one 384-well `qc_review/` + `CN_review/`,
  with each well still labelled by its subplate and a tab strip
  (`All | SL1 | SL2 | SL3 | SL4`) giving the per-subplate view over the same map.
- Additive, not a rewrite: 96-mode paths render character-identically (`SUB_BASE`
  carries no wildcard at 96), so existing plates neither move nor re-run. The
  decisions CSV gained an optional `subplate` column; the 5-column form stays valid.
- Because AneuFinder is per-cell, already-processed subplates are reused rather
  than recomputed (`aneufinder.plate384_models: auto`), so a plate whose subplates
  are already done gets a full 384 review in minutes.

Not carried over: a *separate* per-subplate review artefact. The subplate view is a
filter on the plate map rather than four extra HTML files, and each subplate's own
standalone `qc_review/` still exists if it was ever run on its own.

See `README.md` ("Plate formats: 96 and 384") and `AGENT_CONTEXT.md` ("384-Well Mode")
for the design and its invariants.

---

*Original plan, kept for context:*

The pipeline currently assumes one run = one 96-well subplate = one pair of FASTQs.
When a plate has more wells (e.g. 384), sequencing still runs at 96 wells per lane,
so a single plate produces multiple FASTQ pairs (subplates, e.g. `plate17_1`..`plate17_4`),
and each is currently analyzed as a fully separate pipeline run with its own
intermediate files and its own `qc_review`/`cn_review` outputs. There is no
plate-level view today.

As throughput needs grow, we want a run defined at the plate level, not the
subplate level, whenever the plate is a multiple of 96 wells:

- user-driven: user sets the expected well count for a run
- auto-detected: pipeline discovers the set of FASTQ pairs that belong to the same
  plate/subplate group and treats them as one logical run

Planned direction:

- keep per-subplate FASTQ-to-BAM processing as-is (that boundary matches how
  sequencing actually delivers data), but track subplates as members of one
  parent plate/run instead of independent runs
- roll intermediate files and review outputs up to the plate level so all
  subplates of a run are visible together, while still supporting drill-down to
  a single subplate
- `qc_review` and `cn_review` should offer both a full-plate overview (all wells
  sequenced in the run) and a per-subplate view (how that subplate performed on
  its own), not just one or the other
- preserve current per-subplate output contracts where possible so this is
  additive rather than a rewrite of the existing 96-well path

## 7. Clustering And Tree-Building In `cn_review`

Once a run routinely covers more wells (e.g. a full 384-well plate across its
subplates), `cn_review` should help assess whether cells cluster well as a
plate, not just list per-well copy-number calls.

Planned direction:

- build a tree/clustering view over the reviewed (PASS) wells in a run, using
  their copy-number profiles
- exact tree format (e.g. hierarchical clustering dendrogram vs. a phylogenetic
  tree) is still undecided and needs more thought once real multi-subplate data
  is available
- this depended on [Higher-Throughput, Multi-Subplate Runs](#6-higher-throughput-multi-subplate-runs),
  which is now **implemented** — a full plate's wells are reviewed together, and the
  second AneuFinder pass already produces one plate-level `MODELS/` directory over the
  PASS wells of all four subplates. That is exactly the input a clustering/tree view
  needs, so this item is now unblocked.
- the genome-wide heatmap in `CN_review` is the natural place to hang it; note its
  height is capped at 14000 px (see `render_well_profiles.R`), so a 384-well tree view
  should not assume one readable row per cell at full plate scale

## 8. Principles For Future Changes

When these future upgrades are implemented:

- preserve current output contracts where possible
- prefer additive modules over disruptive rewrites
- keep user-facing launch steps simple
- keep generated outputs and runtime state out of Git
- document new behavior in the README at the same time as implementation
