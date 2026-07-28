# `config/` — QC decisions template & schema

This directory is **documentation only**: it holds the template
(`qc_decisions.csv.example`) and the schema for the human QC-decisions file.

**The live decisions file lives WITH the plate data, per plate:**

```
<PLATE_DIR>/qc_decisions.csv
```

(i.e. in the data directory next to that plate's `barcodes.tsv` and outputs — not
in this pipeline repo). The path is configurable via `qc_review.decisions_file` in
`config.yaml` (relative paths resolve against `<PLATE_DIR>`). Keeping it per-plate
means each plate's decisions stay with its data and never collide between plates.

It records the manual keep / exclude call for each well of a plate after a
reviewer has inspected `<PLATE_DIR>/qc_review/review.html` (plate map + per-well
AneuFinder profile + read-count histogram + automated read-count status).

## How it is produced

1. Run the pipeline through QC + AneuFinder + the first-pass review report:
   ```bash
   sbatch submit_pipeline.sh          # MODE="pre_review"
   ```
2. Open `<PLATE_DIR>/qc_review/review.html` (over Samba or as a local file).
3. Click wells, set a **decision** + **reason(s)** (+ optional notes).
4. Click **Copy terminal save command** and paste it into a terminal on the
   cluster — it writes your decisions straight to `<PLATE_DIR>/qc_decisions.csv`
   (a `cat > … <<'EOF'` heredoc). *(Alternatively use **Download qc_decisions.csv**
   / **Copy CSV only** and place the file at that path yourself.)*
5. Run the post-review phase in one go:
   ```bash
   sbatch submit_pipeline.sh          # MODE="post_review"
   ```
   This validates the CSV, derives the PASS wells, runs the second AneuFinder pass,
   and builds the final deliverable `<PLATE_DIR>/CN_review/cn_review.html`. (The
   individual steps are also
   available as `MODE=validate_review`, `aneufinder_reviewed`, `cn_review` for
   debugging.)

The browser never writes into the project directory — it only offers copyable text
(CSV or the save command) or a download. `<PLATE_DIR>/qc_decisions.csv` is the
reproducible artefact; in this repo only `qc_decisions.csv.example` is tracked.

## Schema

`<PLATE_DIR>/qc_decisions.csv` has a header row and one row per well. **Two schemas
are valid**, and both are accepted indefinitely — a 96-well plate keeps producing the
5-column form:

**96-well (5 columns)** — `qc_decisions.csv.example`

| column      | meaning                                                        |
|-------------|----------------------------------------------------------------|
| `sample_id` | `{plate}_{well}` (e.g. `plate17_4_W01`)                         |
| `well`      | well ID, must be one of the plate's wells (e.g. `W01`)         |
| `decision`  | one of `PASS`, `EXCLUDE`, `REVIEW`, `REPEAT`                    |
| `reason`    | empty, or a `;`-separated list of controlled tokens (see below) |
| `notes`     | free text (quote if it contains a comma)                       |

**384-well (6 columns)** — `qc_decisions_384.csv.example`

Identical, plus a `subplate` column immediately after `sample_id`. A well ID like
`W07` is only unique *within* a subplate, so 384 mode keys on the `(subplate, well)`
pair:

| column      | meaning                                                        |
|-------------|----------------------------------------------------------------|
| `sample_id` | `{subplate}_{well}` (e.g. `plate21_2_W07`)                      |
| `subplate`  | subplate directory name (e.g. `plate21_2`)                      |
| `well`      | well ID within that subplate (e.g. `W07`)                       |
| `decision`  | as above                                                        |
| `reason`    | as above                                                        |
| `notes`     | as above                                                        |

`sample_id` is the same string whether a cell is reviewed in a standalone 96-well run
of its subplate or as part of the whole 384 plate.

**Allowed `reason` tokens** (the `reason` cell may be empty or hold several joined
with `;`): `low_read_count`, `noisy_profile`, `poor_bin_distribution`,
`low_complexity`, `suspected_doublet_or_mixed_well`, `sample_swap_suspected`,
`manual_exception`, `other`, `missing_qc_metric`.

This vocabulary is **frozen** — the CellenONE image QC layer deliberately introduced
no new tokens; it only *suggests* an existing one (see below).

Example with multiple reasons:

```csv
sample_id,well,decision,reason,notes
SAMPLE_002,A02,EXCLUDE,low_read_count;noisy_profile,Weak CN signal
```

The review report writes whichever schema matches the plate it was built for, and
`validate_qc_decisions.py` / `derive_included_wells.py` detect it from the header.
`included_wells.tsv` mirrors the input schema: `sample_id\twell` or
`sample_id\tsubplate\twell`.

## Automatic gating (read-count only)

The review report pre-fills each well's decision from an automatic **read-count-only**
status. `usable_reads` = **mapped reads in the dedup BAM** (`reads mapped` from
`bam/{well}.stats.txt`). Duplication is **not** used for gating (it is shown only as a
labelled diagnostic). Cutoffs live in `config.yaml` under `qc_review`:

| usable_reads | auto_status | default decision | default reason |
|--------------|-------------|------------------|----------------|
| `>= usable_reads_pass_cutoff` (100000) | PASS | PASS | (empty) |
| `>= usable_reads_warn_cutoff` (50000) and `<` pass | WARN | REVIEW | (empty) |
| `< usable_reads_warn_cutoff` | FAIL | EXCLUDE | low_read_count |
| missing | UNKNOWN | REVIEW | missing_qc_metric |

These are only *defaults* in the UI — the human decision always wins.

## How decisions gate the second AneuFinder pass

This CSV is the gate for the **second (post-review) AneuFinder pass**
(`MODE=aneufinder_reviewed`). After it is validated, `derive_included_wells.py`
turns it into `qc_review/included_wells.tsv`, and AneuFinder is rerun on just those
wells into `aneufinder_reviewed/`:

| decision | included in second pass? |
|----------|--------------------------|
| `PASS`   | yes (always)             |
| `REVIEW` | only if `qc_review.include_review: true` in `config.yaml` |
| `EXCLUDE`| no                       |
| `REPEAT` | no — and logged as "flagged for rerun"                    |

The first-pass outputs (`aneufinder/`) and the original BAMs are never modified by
the second pass. A missing `<PLATE_DIR>/qc_decisions.csv` only blocks the post-review
targets — it never blocks normal pre-review pipeline execution.
