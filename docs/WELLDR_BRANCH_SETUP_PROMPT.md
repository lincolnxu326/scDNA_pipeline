> **Record only.** This is the original brief used to create the `welldr` branch
> (branch rename, target README, `docs/WELLDR_PLAN.md`). It has already been carried
> out; do not run it again. Read it for the reasoning and the decisions behind the plan.
> `docs/WELLDR_PLAN.md` is the current source of truth and wins where the two differ.

# Agent prompt: set up the `welldr` branch (docs only)

You are working in the scDNA pipeline repo at `VCAM1_GnT/TOOLS/scDNA_pipeline` (GitHub: `lincolnxu326/scDNA_pipeline`). Your job in this session is limited to git housekeeping and documentation. **Do not change any pipeline code, config, environment or script, and do not run the pipeline.** The output of this session is a new branch that carries (1) a README describing the target pipeline and (2) a phased implementation plan that a later agent will execute one phase at a time.

## 1. Background

The lab is moving from its current single-cell DNA protocol to the wellDR-seq protocol (Wang et al., "Coalescing single-cell genomes and transcriptomes to decode breast cancer progression", Cell 188, 6355-6369, 2025; code at https://github.com/navinlabcode/wellDR-seq). If `docs/references/wellDR-seq_Wang2025_Cell.pdf` exists, read its STAR Methods ("wellDR-seq procedure" and "Data preprocessing") before writing anything.

The current pipeline (branch `main`, to be renamed `dlp+`) is described in `README.md`, `AGENT_CONTEXT.md`, `FUTURE_PLANS.md` and `docs/01_METHODS.md`. Read all four and skim `Snakefile` and `config.yaml` first.

### 1.1 Decisions already made by the user

1. The old protocol and the new protocol live on **separate branches**. The new protocol is **not** a config switch or mode inside one pipeline. Code that only serves the old protocol is removed on the new branch (in later phases, not now).
2. Branch names: the current pipeline becomes `dlp+`; the new pipeline is `welldr`.
3. One run = one plate = one DNA library + one RNA library. No multiplexing of several plates/chips in one library.
4. Plates will be 384 wells at first and 5,184 wells (72 x 72) later. Cells are still dispensed on CellenONE.
5. The human-in-the-loop QC review (`review.html`, `qc_decisions.csv`, two-pass AneuFinder, `CN_review`) is kept.
6. CNV calling stays with AneuFinder. CNV results are later added to the RNA layer of the same cells.
7. **Open, do not decide:** how CellenONE names well positions for 384 and 5,184 formats and how they map to `well_map.tsv`. Record it as an open item.
8. **Open, do not decide:** sequencing run configuration (instrument, read lengths). Record the minimum requirements below as requirements, not as known values.

### 1.2 Protocol facts the plan must use

DNA library
- Fragmentation by Tn5 tagmentation (Illumina TDE1) of protease-lysed cells. Fragment ends are random, so position-based duplicate marking is valid without UMIs. The authors used bowtie2 and `sambamba markdup`.
- Well identity is in the index reads: DNA-CB1 is the i7 index, DNA-CB2 is the i5 index (8 bp each). In the published 72 x 72 chip, i7 encodes the column and i5 the row. `bcl2fastq` with 0 index mismatches produces one FASTQ pair per well, so no in-pipeline DNA demultiplexing is needed.
- Combinatorial (non-unique) dual indexing means index hopping moves reads between wells that share a row or a column.
- Author DNA QC: exclude MAPQ < 1; exclude cells under 100K reads when the plate mean is about 500K; exclude cells with more than 10% empty bins. CNA profiles are stable from about 250K reads per cell at 220 kb bins (Spearman 0.92) and about 500K unique reads at 50 kb.

RNA library
- Oligo-dT RT with RNA-CB1, template switching (biotinylated TSO), enrichment PCR with RNA-CB2, tagmentation of amplified cDNA. 3'-anchored library.
- Read 1 uses a custom sequencing primer (WDR_Read1). Read 1 structure (0-based):
  - bp 0-7: RNA-CB2 (column in the published chip)
  - bp 8-17: 10 bp random sequence
  - bp 18-31: fixed adapter `GAGGCGTAGTGGCT`
  - bp 32-39: RNA-CB1 (row in the published chip)
  - bp 40 onwards: polyT
- Read 2 is cDNA.
- The RNA FASTQ is split by library only, so the pipeline must demultiplex RNA wells from Read 1.
- **There is no UMI.** The 10 bp random sequence is added during the cDNA enrichment PCR, not at RT, so it does not identify original molecules and must not be used as a UMI. It may be kept in the read name for QC.
- Author processing: Trimmomatic (adapters: `GTACTCTGCGTTGATACCACTGCTT`, polyA, Nextera N7 `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`, Nextera S5 `CTGTCTCTTATACACATCTGACGCTGCCGACGA`), then STARsolo `--soloType SmartSeq --soloUMIdedup Exact NoDedup --soloStrand Unstranded --soloFeatures Gene GeneFull`.
- Author RNA QC (tumour): keep cells with 200 to 10,000 genes and mitochondrial fraction at most 30%.

Integration in the paper (to be reproduced with AneuFinder output instead of varbin/CopyKit)
- DNA subclones from UMAP + dbscan on segment ratios.
- Gene dosage: subclone pseudobulk expression (Seurat `AggregateExpression`, DESeq2 `vst`) correlated (Pearson) with subclone integer copy number per gene.
- cis/trans DE between subclones (fold change > 1.4, adjusted p < 0.05).
- RNA-inferred CNV (CopyKAT, inferCNV) as a concordance check against DNA calls.

Run requirements (values not yet known)
- DNA: index reads 8 + 8 cycles. The i5 orientation in the sample sheet depends on the instrument and must be recorded explicitly and verified from the Undetermined index counts.
- RNA: Read 1 at least 48 cycles with the custom WDR_Read1 primer, so that both barcodes and the start of polyT are read.

## 2. Tasks

Do these in order. Stop and report to the user at each point marked **STOP**.

### 2.1 Git housekeeping

1. Run `git status` and `git branch -a`. `main` currently has uncommitted modifications (at least `Snakefile`, `README.md`, `FUTURE_PLANS.md`, `.gitignore`, `.Renviron`, `.Rprofile`, `.renvignore`, `Seqinfo/*`). **STOP**: show the user the list with a one-line summary of each diff and ask whether to commit them to `main` before renaming, stash them, or leave them. Never discard them.
2. Rename the local branch `main` to `dlp+` (`git branch -m main 'dlp+'`; quote the name in shell commands because of the `+`).
3. Create `welldr` from `dlp+` and switch to it.
4. Leave `feature/umi-tools-dedup` and `worktree-no-cell-image-mode` untouched.
5. Do not push yet (see 2.5).

### 2.2 Rewrite `README.md` on `welldr`

Describe the **target** pipeline. Put a status block at the top stating that the branch is under development, that the code in it is still the `dlp+` pipeline until the phases in `docs/WELLDR_PLAN.md` land, and that the legacy protocol lives on branch `dlp+`.

Keep the structure and tone of the current README (contents list, quick start, stage table, inputs, outputs, modes, review workflow, configuration, environments, layout, further documentation). Content to cover:

1. Inputs per plate: `well_map.tsv` (the only file that knows plate geometry), `dna/fastq/{well}_R1/R2.fastq.gz` from bcl-convert, `rna/{plate}_RNA_R1/R2.fastq.gz`.
2. `well_map.tsv` schema: `well_id, row, col, dna_i7, dna_i5, rna_cb1, rna_cb2, cellenone_pos`. `cellenone_pos` is marked as open (see 1.1 item 7). `well_id` format must work for 72 rows; propose a format and mark it as a proposal.
3. Stage table, DNA branch: check DNA inputs, adapter trimming, FastQC, bowtie2 alignment, position-based duplicate marking, library complexity, index hopping QC (reads in empty wells by row and column), AneuFinder pass 1.
4. Stage table, RNA branch: Read 1 demultiplexing (adapter and polyT check, up to 1 mismatch per barcode), trimming, STARsolo SmartSeq (Exact and NoDedup), RNA QC.
5. Review: unchanged flow. Additions: duplication and complexity metrics replace UMI retention, an RNA panel per well, and a generic rows x cols plate map driven by `well_map.tsv`. Sidecar assets are always used at 5,184 wells.
6. Post-review: AneuFinder pass 2, `CN_review`, and integration (gene to bin map, per-cell gene copy number, combined RNA + CN object joined on `well_id`).
7. Output tree per plate (`trimmed/`, `raw_bam/`, `bam/`, `markdup/`, `rna/`, `aneufinder/`, `cellenone/`, `qc_review/`, `aneufinder_reviewed/`, `CN_review/`, `integrated/`). `bam/{well}.bam` remains the AneuFinder handoff.
8. Modes: `pre_review` (DNA + RNA processing + review report), `post_review` (pass 2 + CN_review + integration). `PLATE_FORMAT` is removed.
9. Run requirements from 1.2.
10. Open items: CellenONE position mapping, run configuration, `well_id` format, rebuilt mappability reference and blacklist from new-protocol euploid cells, re-tuned read cutoffs.

### 2.3 Write `docs/WELLDR_PLAN.md`

This file is the implementation brief for future agents. Write it as instructions to an agent. Structure:

1. Objective
2. Protocol facts (copy section 1.2 of this prompt, tidied)
3. Decisions and open items (section 1.1)
4. Change inventory, a table with columns Item / Status (Remove, New, Change, Keep) / Reason, covering at least:
   - Remove: `demultiplex` rule and `extract_umi_barcode.py`; `filter_dimers` and `filter_adapter_dimers.py`; `deduplicate_bam` (umi_tools) and `summarize_umi_tools_dedup.py`; the 384 subplate machinery (subplate formula in `plate384_layout.py`, `link_subplate_models.py`, `aneufinder.plate384_models`, SL1-SL4 tabs, `PLATE_FORMAT`); the `preprocessing:` UMI/barcode/trim block in `config.yaml`.
   - New: `check_dna_inputs`, `trim_dna`, `markdup`, `dna_complexity`, `index_hopping_qc`, `rna_demux`, `trim_rna`, `starsolo`, `rna_qc`, `gene_bin_map`, `cn_to_genes`, `build_object`, helper `make_samplesheet.py`, env `rna.yaml`, STAR hg38 index.
   - Change: `align` (trimmed input, simpler read groups), well wildcard (`W\d+` to `well_map.tsv` IDs), `generate_qc_review.py` (generic grid, new metrics, RNA panel, 5,184 scale), `ingest_cellenone.py` (positions from `well_map.tsv`, pending the open item), `render_well_profiles.R` and `cn_review` (5,184 scale, heatmap height cap), `multiqc` (fastp, markdup, STARsolo logs), mappability reference and blacklist inputs, `usable_reads` and `min_reads_for_model` cutoffs, `submit_pipeline.sh` modes.
   - Keep: `run_aneufinder.R` (do not simplify; see the GC and bin-naming warnings in `AGENT_CONTEXT.md`), `check_gc_rds`, `generate_blacklist`, the two-pass flow, `validate_qc_decisions`, `derive_included_wells`, `fastqc`, CellenONE image calls.
5. Phases. For each phase give: scope, files touched, tests or checks, what to show the user (figures first, tables as backing data), and a hard **stop for review** before the next phase.
   - Phase 0: protocol spec (`docs/PROTOCOL_wellDR.md`), `well_map.tsv` schema and validator, `make_samplesheet.py`, a pilot plate. CellenONE mapping stays open.
   - Phase 1: remove legacy-only code; DNA input layer driven by `well_map.tsv`. Check: reads per well on a plate map, Undetermined fraction, reads in empty wells.
   - Phase 2: trimming, alignment, duplicate marking, complexity, index hopping QC. Check: duplication rate against depth, agreement of `samtools markdup` and `sambamba markdup` on a subset of wells.
   - Phase 3: AneuFinder on new data, rebuilt mappability reference and blacklist, GC template check, re-tuned cutoffs, review report on the new grid and at 5,184 scale. Check: profiles of known euploid and aneuploid cells, review report renders.
   - Phase 4: RNA branch. Check: genes and counts per well, mitochondrial fraction, DNA and RNA recovery per well on one plate map.
   - Phase 5: integration. Check: per-gene CN matrix sanity, gene dosage correlations, RNA-inferred CNV against DNA calls.
6. Principles: additive within the branch where possible, output contracts documented in the README at the same time as implementation, generated data out of git, one phase per session.

### 2.4 Small pointer edits

1. On `welldr`, add a short section at the top of `AGENT_CONTEXT.md` saying this branch targets wellDR-seq, that `docs/WELLDR_PLAN.md` is the source of truth for the migration, and that legacy behaviour lives on `dlp+`. Do not rewrite the rest of the file.
2. Do not edit any file on `dlp+` beyond what the user approves in 2.1.

### 2.5 Commit and report

1. Commit on `welldr` with a clear message. Changed files must be limited to `README.md`, `docs/WELLDR_PLAN.md` and `AGENT_CONTEXT.md`. Verify with `git diff --stat 'dlp+'..welldr`.
2. **STOP**: report the branch layout, the diff stat, and the open items. Ask before pushing. On approval, push `dlp+` and `welldr` to `origin`. **Do not** delete `origin/main` or change the GitHub default branch unless the user explicitly asks; tell them that renaming the default branch on GitHub is a separate step they control.

## 3. Writing rules

- Plain scientific register. Numbered sections. No em dashes. No promotional or "AI style" phrasing.
- State facts with their source (paper section, author repo file, or this repo's file). Mark proposals as proposals and open items as open.
- Do not invent values for the open items.
