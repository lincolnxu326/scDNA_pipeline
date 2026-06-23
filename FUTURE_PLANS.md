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

## 3. Performance Improvements

The current preprocessing path is slower than it should be, especially in demultiplexing and plate-level preprocessing.

Planned direction:

- optimize demultiplexing, which is currently dominated by single-process Python and gzip I/O
- review UMI-tools deduplication metrics after test plate runs
- split adapter filtering into per-well jobs for the same reason
- review alignment resource mapping so allocated SLURM CPUs match requested tool threads
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

## 5. Principles For Future Changes

When these future upgrades are implemented:

- preserve current output contracts where possible
- prefer additive modules over disruptive rewrites
- keep user-facing launch steps simple
- keep generated outputs and runtime state out of Git
- document new behavior in the README at the same time as implementation
