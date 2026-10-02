# scDNA Pipeline: wellDR-seq branch

> **Status: under development.**
>
> This branch (`welldr`) targets the wellDR-seq protocol (Wang et al., Cell 188,
> 6355-6369, 2025). This README describes the **target** pipeline. Until the phases in
> [`docs/WELLDR_PLAN.md`](docs/WELLDR_PLAN.md) are implemented, the code on this branch
> is still the `dlp+` pipeline and will not run wellDR-seq data.
>
> The legacy protocol (NlaIII digestion, Y-adapter ligation, in-pipeline demultiplexing,
> UMI-tools deduplication, 384-well subplate mode) lives on branch `dlp+`.

Snakemake workflow for one wellDR-seq plate. Each plate produces one DNA library and
one RNA library from the same cells. The DNA branch aligns per-well FASTQs, marks
duplicates by position and calls copy number with AneuFinder. The RNA branch
demultiplexes wells from Read 1 and counts genes with STARsolo. A human review step
sits between a first and a second AneuFinder pass. After review, copy number is
joined to the RNA layer of the same cells.

The repository holds code only. All plate data, intermediate files and results live
in the plate directory you point the pipeline at, outside version control.

## Contents

1. [Quick start](#1-quick-start)
2. [What the pipeline does](#2-what-the-pipeline-does)
3. [Inputs](#3-inputs)
4. [CellenONE run folder](#4-cellenone-run-folder)
5. [Outputs](#5-outputs)
6. [Running modes](#6-running-modes)
7. [Review workflow](#7-review-workflow)
8. [Post-review and integration](#8-post-review-and-integration)
9. [Sequencing run requirements](#9-sequencing-run-requirements)
10. [Configuration reference](#10-configuration-reference)
11. [Environments](#11-environments)
12. [Repository layout](#12-repository-layout)
13. [Open items](#13-open-items)
14. [Further documentation](#14-further-documentation)

## 1. Quick start

```bash
# 1. Check the setup
bash pre_check.sh

# 2. First pass: DNA and RNA processing, AneuFinder on all wells, review report. Stops here.
MODE=pre_review PLATE_DIR=/path/to/plateNN sbatch submit_pipeline.sh

# 3. Open <PLATE_DIR>/qc_review/review.html, make decisions, save qc_decisions.csv

# 4. Second pass: AneuFinder on PASS wells, final copy-number viewer, RNA + CN integration
MODE=post_review PLATE_DIR=/path/to/plateNN sbatch submit_pipeline.sh
```

Plate size (384 or 5,184 wells) is read from `well_map.tsv`. There is no
`PLATE_FORMAT` setting.

To run several plates at the same time from this directory, add
`SNAKEMAKE_EXTRA=--nolock`. Each plate writes only inside its own `PLATE_DIR`.

## 2. What the pipeline does

### 2.1 DNA branch

| # | Stage | Tool | Output folder |
|---|-------|------|---------------|
| D1 | Check DNA inputs against `well_map.tsv` (missing, empty and unexpected wells) | Python | `qc/` |
| D2 | Adapter trimming (Nextera) | fastp | `trimmed/` |
| D3 | Read QC | FastQC, MultiQC | `fastqc/`, `multiqc/` |
| D4 | Alignment | bowtie2, samtools | `raw_bam/` |
| D5 | Position-based duplicate marking | samtools markdup | `markdup/`, `bam/` |
| D6 | Library complexity per well | Python | `markdup/` |
| D7 | Index hopping QC: reads in empty wells, by row and by column | Python | `qc/` |
| D8 | First AneuFinder pass, all wells | AneuFinder | `aneufinder/` |

Duplicates are marked by alignment position, without UMIs. This is valid because
Tn5 fragment ends are random (paper, STAR Methods, "Data preprocessing"; the authors
used `sambamba markdup`).

### 2.2 RNA branch

| # | Stage | Tool | Output folder |
|---|-------|------|---------------|
| R1 | Demultiplex wells from Read 1: check fixed adapter and polyT, match RNA-CB1 and RNA-CB2 with up to 1 mismatch each | Python | `rna/demux/` |
| R2 | Trimming (adapters, polyA, Nextera N7 and S5) | Trimmomatic or fastp | `rna/trimmed/` |
| R3 | Gene counting, SmartSeq mode | STARsolo | `rna/starsolo/` |
| R4 | RNA QC per well: reads, genes, counts, mitochondrial fraction | Python or R | `rna/qc/` |

The RNA library has **no UMI**. The 10 bp random sequence in Read 1 is added during
the cDNA enrichment PCR, not at reverse transcription, so it does not identify
original molecules. STARsolo runs with `--soloUMIdedup Exact NoDedup`. The random
sequence is kept in the read name for QC only.

### 2.3 Shared stages

| # | Stage | Output folder |
|---|-------|---------------|
| S1 | CellenONE image ingest and rendering (optional, display-only) | `cellenone/` |
| S2 | Static HTML review report | `qc_review/` |
| S3 | Human review, saved as a CSV | `qc_decisions.csv` |
| S4 | Second AneuFinder pass, PASS wells only | `aneufinder_reviewed/` |
| S5 | Final copy-number viewer | `CN_review/` |
| S6 | RNA + copy-number integration | `integrated/` |

`MODE=pre_review` runs D1 to D8, R1 to R4, S1 and S2. `MODE=post_review` runs S4 to S6.

## 3. Inputs

### 3.1 Plate directory

One directory per plate. One plate is one DNA library and one RNA library.

```text
plateNN/
├── well_map.tsv                         required
├── dna/
│   └── fastq/
│       ├── <well_id>_R1.fastq.gz        one pair per well, from bcl-convert
│       └── <well_id>_R2.fastq.gz
├── rna/
│   ├── plateNN_RNA_R1.fastq.gz          one pair for the whole plate
│   └── plateNN_RNA_R2.fastq.gz
└── <CellenONE .Run folder>              optional, may also live elsewhere
```

DNA wells are already separated by the sequencing facility's demultiplexing on the
i7 and i5 index reads, so the pipeline does not demultiplex DNA. Name the per-well
FASTQs (or symlinks to them) `<well_id>_R1.fastq.gz` and `<well_id>_R2.fastq.gz`.

The RNA FASTQ is split by library only. The pipeline demultiplexes RNA wells from
Read 1 (stage R1).

### 3.2 well_map.tsv

`well_map.tsv` is the only file that knows the plate geometry. Every per-well rule,
the review plate map and the CellenONE mapping read it. Tab-separated, with a header,
one row per well.

| Column | Meaning |
|--------|---------|
| `well_id` | Well identifier used in every file name and in `sample_id` |
| `row` | Row index, 1-based |
| `col` | Column index, 1-based |
| `dna_i7` | DNA-CB1, 8 bp, as read in the i7 index read |
| `dna_i5` | DNA-CB2, 8 bp, as read in the i5 index read (orientation: see section 9) |
| `rna_cb1` | RNA-CB1, 8 bp, as it appears in Read 1 |
| `rna_cb2` | RNA-CB2, 8 bp, as it appears in Read 1 |
| `cellenone_pos` | CellenONE position token for this well. **Open item, see section 13.** |

**Proposal for `well_id`:** `R<row>C<col>` with two-digit zero-padded row and column,
for example `R01C01` to `R72C72`. Letter rows (A to P) do not extend to 72 rows. The
`sample_id` used in review files is `<plate>_<well_id>`, for example `plate30_R05C17`.
This format is a proposal and is not yet fixed.

Example (barcodes are placeholders):

```text
well_id	row	col	dna_i7	dna_i5	rna_cb1	rna_cb2	cellenone_pos
R01C01	1	1	NNNNNNNN	NNNNNNNN	NNNNNNNN	NNNNNNNN	(open)
R01C02	1	2	NNNNNNNN	NNNNNNNN	NNNNNNNN	NNNNNNNN	(open)
```

In the published 72 x 72 design the i7 index (DNA-CB1) and RNA-CB2 encode the column,
and the i5 index (DNA-CB2) and RNA-CB1 encode the row. The paper text states only
"72 (row index) by 72 (column index)" barcodes; the axis assignment is to be verified
against Table S5 in Phase 0. The pipeline does not rely on it: it reads every well's
barcodes from `well_map.tsv`.

`workflow/scripts/make_samplesheet.py` (planned) writes the bcl-convert sample sheet
for the DNA library from `well_map.tsv`.

### 3.3 Read structure

**DNA library.** Paired-end genomic reads. Well identity is in the index reads only.

| Read | Content |
|------|---------|
| i7 index (8 bp) | DNA-CB1 |
| i5 index (8 bp) | DNA-CB2 |
| R1, R2 | Tn5-tagmented genomic DNA, Nextera adapters at the 3' end of short inserts |

**RNA library.** Read 1 uses the custom sequencing primer WDR_Read1. Positions are
0-based.

| Read 1 bases | Content |
|--------------|---------|
| 0 to 7 | RNA-CB2 |
| 8 to 17 | 10 bp random sequence (not a UMI) |
| 18 to 31 | fixed adapter `GAGGCGTAGTGGCT` |
| 32 to 39 | RNA-CB1 |
| 40 onward | polyT |

Read 2 is cDNA. The library is 3'-anchored.

### 3.4 Reference files

Set in `config.yaml`. All must exist before a run.

| Key | What it is |
|-----|------------|
| `genome.fasta` | Reference FASTA, hg38. Chromosomes named `chr1` to `chr22`, `chrX`, `chrY`. |
| `genome.index_prefix` | bowtie2 index for that FASTA |
| `rna.star_index` | STAR index for hg38 with a matching GTF (planned) |
| `rna.gtf` | Gene annotation used by STARsolo and by the gene-to-bin map (planned) |
| `mappability.reference_bam` | Euploid reference BAM for blacklist diagnostics. **To be rebuilt from wellDR-seq euploid cells, see section 13.** |
| GC template RDS | `resources/reference/hg38_binsize1000000_variable_bins_with_GC.rds`, or set `aneufinder.gc_rds`. Must match `aneufinder.binsize`. |

The authors aligned to hg19. This pipeline uses hg38.

## 4. CellenONE run folder

Cells are dispensed on CellenONE. The CellenONE layer is optional and display-only:
it never changes a well's automatic status or default decision.

Map the plate name to its `.Run` folder in `config.yaml`:

```yaml
cellenone:
  enable: true
  runs:
    plateNN: "/path/to/<run name>.Run"
```

The folder contents, the image call (`SINGLE`, `PASS`, `CONTAMINATION`, `FAIL`,
`NO_OBJECT`, `NO_IMAGE`) and the `cellenone/` output files are unchanged from the
`dlp+` branch and are described in
[`docs/CELLENONE_AND_QC.md`](docs/CELLENONE_AND_QC.md). In short, the folder needs
`Reordered_*_isolated.xls` and the `*_Printed_*_(<POS>)_Trans_*.png` images; the
fluorescence images, `geoprops.xls`, `BackgroundEjZone` and `cellenREPORT/` are
optional.

**Open item:** how CellenONE names well positions for 384-well and 5,184-well
targets, and how those names map to `well_map.tsv` (`cellenone_pos`). Until this is
settled, images are matched only for layouts where the mapping has been verified.

## 5. Outputs

Everything is written inside `PLATE_DIR`. Nothing is written into the repo except
logs under `logs/`.

```text
plateNN/
├── well_map.tsv, dna/, rna/            inputs (section 3)
├── qc/                                 input check, index hopping QC
├── trimmed/                            per-well trimmed DNA FASTQs
├── fastqc/                             FastQC per well
├── multiqc/                            MultiQC over fastp, markdup and STARsolo logs
├── raw_bam/                            aligned DNA BAMs before duplicate marking
├── markdup/                            duplicate-marking metrics, complexity per well
├── bam/                                <well_id>.bam (+ .bai, .stats.txt, .flagstat.txt)
├── mappability/                        blacklist diagnostic plots
├── rna/
│   ├── demux/                          per-well RNA FASTQs and demux stats
│   ├── trimmed/
│   ├── starsolo/                       STARsolo output: Gene and GeneFull matrices
│   └── qc/                             per-well RNA metrics
├── aneufinder/                         first pass, all wells (MODELS/, profiles, heatmap)
├── cellenone/                          only if CellenONE is configured
├── qc_review/
│   ├── review.html                     review report (open this)
│   ├── plots/
│   └── assets/
├── qc_decisions.csv                    written by you after review
├── aneufinder_reviewed/                second pass, PASS wells only
├── CN_review/
│   ├── cn_review.html                  final copy-number viewer
│   ├── plots/
│   └── genome_heatmap.png
├── integrated/                         gene-to-bin map, per-cell gene CN, RNA + CN object
└── logs/
```

**AneuFinder handoff.** `bam/<well_id>.bam` is the file AneuFinder reads. It is the
coordinate-sorted, duplicate-marked BAM. This contract is the same as on `dlp+`
(where it was the UMI-deduplicated BAM), so `run_aneufinder.R` does not change.

## 6. Running modes

Set `MODE` at submit time.

| Mode | What it runs |
|------|--------------|
| `pre_review` | DNA and RNA processing, first AneuFinder pass, review report, then stops. Default. |
| `post_review` | Validate decisions, second AneuFinder pass, `CN_review`, integration. Needs `qc_decisions.csv`. |
| `dna` | DNA branch D1 to D7 only (planned) |
| `rna` | RNA branch R1 to R4 only (planned) |
| `review` | Rebuild `review.html` only |
| `validate_review` | Check `qc_decisions.csv` only |
| `aneufinder_reviewed` | Second AneuFinder pass only |
| `cn_review` | Rebuild `cn_review.html` only |
| `integrate` | Integration stage only (planned) |
| `cellenone` | CellenONE ingest and images only |
| `blacklist` | Blacklist diagnostic plots only |

Submit-time variables:

| Variable | Default | Purpose |
|----------|---------|---------|
| `PLATE_DIR` | set in script | Plate directory |
| `SNAKEMAKE_EXTRA` | empty | Extra Snakemake flags, for example `--nolock` |

`PLATE_FORMAT` and `SUBPLATES` are removed on this branch.

## 7. Review workflow

The review flow is the same as on `dlp+`:

1. Open `<PLATE_DIR>/qc_review/review.html`.
2. Click a well to see its copy-number profile, read-count histogram, metrics,
   RNA panel and (if configured) cell images.
3. Choose PASS, EXCLUDE, REVIEW or REPEAT, add reasons and notes, and click
   **Save decision**.
4. Click **Copy terminal save command** and paste it into a cluster terminal. This
   writes `<PLATE_DIR>/qc_decisions.csv`.
5. Run `MODE=post_review`.

Changes on this branch:

- **Plate map.** A generic rows x columns grid built from `well_map.tsv`, replacing
  the fixed 96-well grid and the 384 subplate tabs.
- **DNA metrics.** Duplication rate and library complexity from `markdup/` replace
  UMI-tools retention.
- **RNA panel.** Per-well genes detected, RNA reads, mitochondrial fraction.
- **Scale.** At 5,184 wells the report always uses sidecar `plots/` and `assets/`
  rather than embedding images. Copy the whole `qc_review/` folder if you move it.

### 7.1 Automatic status

Each well gets an automatic status from `usable_reads` that pre-fills the decision.
The reviewer can change any decision.

| usable_reads | Status | Default decision |
|--------------|--------|------------------|
| at or above `usable_reads_pass_cutoff` | `PASS` | PASS |
| between the two cutoffs | `WARN` | REVIEW |
| below `usable_reads_warn_cutoff` | `FAIL` | EXCLUDE |
| missing | `UNKNOWN` | REVIEW |

On this branch `usable_reads` is mapped, non-duplicate reads with MAPQ of at least 1
in `bam/<well_id>.bam` (proposal). The cutoffs carried over from `dlp+` (100,000 and
50,000) are placeholders until they are re-tuned on wellDR-seq data (section 13). For
reference, the authors excluded cells under 100K reads when the plate mean was about
500K, and cells with more than 10% empty bins.

RNA metrics are shown in the review but do not set the automatic status.

## 8. Post-review and integration

`MODE=post_review` runs:

1. Validation of `qc_decisions.csv` and derivation of the PASS wells
   (`qc_review.include_review: true` also includes REVIEW wells).
2. Second AneuFinder pass on the PASS wells, into `aneufinder_reviewed/`.
3. The final viewer `CN_review/cn_review.html`.
4. Integration, into `integrated/` (file names are proposals):

| File | Contents |
|------|----------|
| `gene_bin_map.tsv` | Each gene from the GTF assigned to its AneuFinder bin |
| `cell_gene_cn.tsv.gz` | Per-cell integer copy number for each gene, from the second-pass models |
| `<plate>_rna_cn.rds` | Combined object: RNA counts from STARsolo plus per-cell copy number, joined on `well_id` |

Downstream analyses planned on top of `integrated/`, following the paper (see
`docs/WELLDR_PLAN.md`): DNA subclones from UMAP and dbscan on segment values, gene
dosage correlation between subclone pseudobulk expression and subclone copy number,
cis and trans differential expression between subclones, and RNA-inferred copy
number (CopyKAT, inferCNV) as a concordance check against the DNA calls.

## 9. Sequencing run requirements

The run configuration (instrument, read lengths) is not yet fixed. These are the
minimum requirements the pipeline depends on.

**DNA library**

- Index reads: 8 cycles i7 and 8 cycles i5.
- Demultiplex with 0 index mismatches, giving one FASTQ pair per well.
- The orientation of the i5 sequence in the sample sheet depends on the instrument.
  Record it explicitly with each run and verify it from the Undetermined index counts
  before running the pipeline.
- Indexing is combinatorial (non-unique dual): wells share i7 values along one axis
  and i5 values along the other. Index hopping therefore moves reads between wells
  that share a row or a column. Stage D7 measures this.

**RNA library**

- Read 1 of at least 48 cycles with the custom WDR_Read1 primer, so that both cell
  barcodes and the start of polyT are read (8 + 10 + 14 + 8 = 40 bp, plus 8 bp polyT).
- The authors also spiked in a custom index 2 primer (WDR_Idx5) and sequenced on a
  NextSeq 2000 (paper, STAR Methods step 7).

## 10. Configuration reference

All settings are in `config.yaml`. Keys marked (planned) do not exist yet.

| Key | Notes |
|-----|-------|
| `well_map` (planned) | Path to `well_map.tsv`, default `<PLATE_DIR>/well_map.tsv` |
| `dna.trim` (planned) | fastp adapter settings |
| `dna.min_mapq` (planned) | MAPQ threshold for `usable_reads`, proposed 1 |
| `rna.demux` (planned) | Adapter, polyT check, barcode mismatches (1) |
| `rna.star_index`, `rna.gtf` (planned) | STAR index and annotation |
| `rna.qc` (planned) | Gene-count and mitochondrial-fraction thresholds |
| `aneufinder.binsize` | Must match the GC template |
| `aneufinder.min_reads_for_model` | Wells below this are not given to AneuFinder. To be re-tuned. |
| `qc_review.usable_reads_pass_cutoff`, `usable_reads_warn_cutoff` | To be re-tuned |
| `qc_review.include_review` | Include REVIEW wells in the second pass |
| `qc_review.embed_assets` | Always sidecar at 5,184 wells |
| `cellenone.runs` | Plate name to `.Run` folder |
| `resources.<rule>` | SLURM threads, memory, time, partition per rule |

The `preprocessing:` block (UMI, barcode and trim offsets) and
`aneufinder.plate384_models` are removed on this branch.

## 11. Environments

- Snakemake itself runs from the `snakemake_scDNA` conda environment set in
  `submit_pipeline.sh`.
- Each rule uses a conda environment from `workflow/envs/`, shared under
  `.snakemake/conda` at the repo root. The RNA branch adds `workflow/envs/rna.yaml`
  (planned).
- R analysis code uses `renv`. Restore it with `R -e 'renv::restore()'`.

## 12. Repository layout

```text
scDNA_pipeline/
├── Snakefile
├── config.yaml
├── multiqc_config.yaml
├── submit_pipeline.sh       SLURM entry point
├── pre_check.sh             pre-flight checks
├── config/                  qc_decisions.csv schema and templates
├── docs/                    protocol notes, migration plan, CellenONE notes
├── resources/reference/     GC template and mappability resources
├── workflow/
│   ├── envs/
│   └── scripts/
│       ├── analysis/
│       ├── reporting/
│       └── rna/             (planned)
├── ad_hoc_checks/           one-off analyses, not part of the pipeline
└── renv/, renv.lock
```

## 13. Open items

These are not decided. Do not fill them in without the user.

1. **CellenONE position mapping.** How CellenONE names well positions for 384-well
   and 5,184-well targets, and how they map to `well_map.tsv` (`cellenone_pos`).
2. **Sequencing run configuration.** Instrument, read lengths and i5 orientation.
   Only the minimum requirements in section 9 are known.
3. **`well_id` format.** `R01C01` is a proposal.
4. **Mappability reference and blacklist.** To be rebuilt from euploid cells
   sequenced with wellDR-seq, because the current reference comes from the old
   protocol.
5. **Read cutoffs.** `usable_reads_pass_cutoff`, `usable_reads_warn_cutoff` and
   `min_reads_for_model` to be re-tuned on wellDR-seq data.
6. **AneuFinder bin size.** The pipeline uses 1 Mb; the authors used about 220 kb
   variable bins. A change needs a matching GC template.

## 14. Further documentation

- [`docs/WELLDR_PLAN.md`](docs/WELLDR_PLAN.md): phased implementation plan for this branch
- [`docs/WELLDR_BRANCH_SETUP_PROMPT.md`](docs/WELLDR_BRANCH_SETUP_PROMPT.md): the original brief behind the plan (record only)
- [`AGENT_CONTEXT.md`](AGENT_CONTEXT.md): start here if you use a coding agent on this repo
- [`docs/CELLENONE_AND_QC.md`](docs/CELLENONE_AND_QC.md): how to read the review report and the CellenONE layer
- [`config/README.md`](config/README.md): `qc_decisions.csv` schema and decision rules
- Wang et al. (2025), Cell 188, 6355-6369; author code at https://github.com/navinlabcode/wellDR-seq
