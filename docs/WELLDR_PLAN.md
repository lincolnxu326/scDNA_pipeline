# wellDR-seq migration plan

This file is the implementation brief for agents working on branch `welldr`. Read it
in full before starting any phase. Work on one phase per session and stop for user
review at the end of each phase.

## 1. Objective

Convert this pipeline from the legacy single-cell DNA protocol (branch `dlp+`) to the
wellDR-seq protocol (Wang et al., "Coalescing single-cell genomes and transcriptomes
to decode breast cancer progression", Cell 188, 6355-6369, 2025; author code at
https://github.com/navinlabcode/wellDR-seq).

The finished pipeline processes one plate (one DNA library and one RNA library) into:

1. per-well duplicate-marked DNA BAMs and AneuFinder copy-number calls,
2. a per-well RNA count matrix from STARsolo,
3. the existing human review report and two-pass AneuFinder flow,
4. a combined RNA + copy-number object for the cells that passed review.

The target behaviour and output contracts are described in `README.md` on this
branch. Keep the README in step with the code as each phase lands.

Read before starting: `README.md`, `AGENT_CONTEXT.md`, `docs/01_METHODS.md`,
`docs/CELLENONE_AND_QC.md`, and the paper at
`docs/references/wellDR-seq_Wang2025_Cell.pdf` (STAR Methods, "wellDR-seq procedure"
and "Data preprocessing"). The PDF is not tracked in git; ask the user for it if it
is missing.

## 2. Protocol facts

Sources are given in brackets. "Paper" means Wang et al. 2025. "Brief" means the
user's instructions for this migration, recorded here; items marked "verify" have
not been confirmed against a primary source.

### 2.1 DNA library

1. Cells are lysed with protease, then the naked DNA is fragmented by Tn5
   tagmentation (Illumina TDE1). [Paper, STAR Methods, procedure step 2]
2. Fragment ends are random, so position-based duplicate marking is valid without
   UMIs. The authors used bowtie2 and `sambamba markdup`. [Paper, "Data
   preprocessing"]
3. Well identity is in the index reads. DNA-CB1 is the i7 index
   (`CAAGCAGAAGACGGCATACGAGAT[8]GTCTCGTGGGCTCGG`), DNA-CB2 is the i5 index
   (`AATGATACGGCGACCACCGAGATCTACAC[8]TCGTCGGCAGCGTC`), 8 bp each. [Paper, procedure
   steps 3 and 4] This holds for the design without chip multiplexing, which is the
   design used here (section 3, decision 3).
4. In the published 72 x 72 chip, i7 encodes the column and i5 the row. [Brief;
   the paper text says only "72 (row index) by 72 (column index)". Verify against
   Table S5.]
5. `bcl2fastq` with 0 index mismatches gives one FASTQ pair per well, so the pipeline
   does no DNA demultiplexing. [Paper, "Data preprocessing"]
6. Indexing is combinatorial (non-unique dual). Index hopping moves reads between
   wells that share a row or a column.
7. Author DNA QC [Paper, "Data preprocessing"]:
   - exclude reads with MAPQ < 1;
   - exclude cells with fewer than 100K reads when the plate mean is about 500K;
   - exclude cells with more than 10% of bins empty.
8. CNA profiles are stable from about 250K reads per cell at 220 kb bins (Spearman
   rho 0.92) and from about 500K unique mapped reads at 50 kb bins (rho 0.92).
   [Paper, Results, and Figures S2C and S2D]
9. The authors aligned to hg19 and called copy number with variable bins (about
   220 kb), lowess GC correction and CBS. This pipeline uses hg38 and AneuFinder.
   [Paper, "Breast tumor DNA data analysis"]

### 2.2 RNA library

1. Oligo-dT reverse transcription with RNA-CB1
   (`GAGGCGTAGTGGCT[8]T30VN`), template switching with a biotinylated TSO
   (`/Bio/TCTCCGACTCAGTACATrGrGrG`), enrichment PCR with RNA-CB2
   (`AAGCAGTGGTATCAACGCAGAGTAC[8][10 random]GAGGCGTAGTGGCT`), then tagmentation of the
   amplified cDNA. The library is 3'-anchored. [Paper, procedure steps 2, 3 and 7]
2. Read 1 uses the custom primer WDR_Read1
   (`GCCTGTCCGCGGAAGCAGTGGTATCAACGCAGAGTAC`). A custom index 2 primer (WDR_Idx5) is
   also spiked in. The authors sequenced on a NextSeq 2000. [Paper, step 7]
3. Read 1 structure, 0-based [derived from the primer sequences in 2.2.1 and 2.2.2;
   Brief]:

   | Bases | Content |
   |-------|---------|
   | 0-7 | RNA-CB2 (column in the published chip; verify) |
   | 8-17 | 10 bp random sequence |
   | 18-31 | fixed adapter `GAGGCGTAGTGGCT` |
   | 32-39 | RNA-CB1 (row in the published chip; verify) |
   | 40 onward | polyT |

4. Read 2 is cDNA.
5. The RNA FASTQ is split by library only. The pipeline must demultiplex RNA wells
   from Read 1. [Paper, "Data preprocessing": "split into FASTQ files of individual
   cells using ... RNA-CB1 and RNA-CB2 in Read 1"]
6. **There is no UMI.** The 10 bp random sequence is part of the RNA-CB2 primer and
   is added in the cDNA enrichment PCR, not at RT, so it does not identify original
   molecules. Do not use it as a UMI. It may be kept in the read name for QC.
   [Paper, step 3 and Discussion: "does not include unique molecular identifiers"]
7. Author RNA processing [Paper, "Data preprocessing"]:
   - Trimmomatic with adapters `GTACTCTGCGTTGATACCACTGCTT` (adaptor_smt2), polyA,
     Nextera N7 `CTGTCTCTTATACACATCTCCGAGCCCACGAGAC`, Nextera S5
     `CTGTCTCTTATACACATCTGACGCTGCCGACGA`;
   - STARsolo `--outSAMtype BAM SortedByCoordinate --soloType SmartSeq
     --soloUMIdedup Exact NoDedup --soloStrand Unstranded --soloFeatures Gene
     GeneFull`.
8. Author RNA QC [Paper]:
   - tumour: keep cells with 200 to 10,000 genes and mitochondrial fraction of at
     most 30% ("Breast tumor single cell RNA data processing");
   - cell lines: 500 to 10,000 genes, mitochondrial fraction of at most 10% ("RNA
     technical metrics comparisons").

### 2.3 Integration in the paper

Reproduce these with AneuFinder output in place of the authors' variable binning and
CopyKit:

1. DNA subclones from UMAP (`uwot`) and dbscan on segment ratios; subclones of fewer
   than 4 cells removed. [Paper, "DNA copy number profiles clustering"]
2. Gene dosage: subclone pseudobulk expression (Seurat `AggregateExpression`, DESeq2
   `vst`) correlated (Pearson) with subclone integer copy number per gene; subclones
   with at least 20 cells with RNA. [Paper, "Global gene dosage analysis"]
3. cis and trans DE between subclones: Seurat `FindMarkers`, fold change > 1.4,
   adjusted p < 0.05; cis if the DE gene lies in a subclonal bin and changes in the
   same direction as copy number. [Paper, "Classification of subclonal CNA cis and
   trans DE genes"]
4. RNA-inferred CNV (CopyKAT, inferCNV) as a concordance check against DNA calls.
   [Paper, "CNA inference from single cell RNA data"]

### 2.4 Run requirements (values not yet known)

1. DNA: index reads of 8 + 8 cycles. The i5 orientation in the sample sheet depends
   on the instrument; record it explicitly per run and verify it from the
   Undetermined index counts.
2. RNA: Read 1 of at least 48 cycles with WDR_Read1, so both barcodes and the start
   of polyT are read (8 + 10 + 14 + 8 = 40 bp, plus 8 bp polyT).

### 2.5 Differences from the paper that are intended

1. Cells are dispensed on CellenONE into plates, not on the ICELL8 cx nanowell
   system used in the paper. The nanolitre volumes in the paper are for ICELL8.
2. Reference is hg38, not hg19.
3. Copy number is called with AneuFinder, not variable binning + CBS.

## 3. Decisions and open items

### 3.1 Decisions made by the user

1. The old and new protocols live on separate branches. The new protocol is not a
   config switch or mode. Code that only serves the old protocol is removed on
   `welldr` (from Phase 1 onward).
2. Branch names: `dlp+` for the legacy pipeline, `welldr` for this one.
3. One run is one plate: one DNA library and one RNA library. No multiplexing of
   several plates or chips in one library.
4. Plates are 384 wells first and 5,184 wells (72 x 72) later. Cells are dispensed on
   CellenONE.
5. The human review (`review.html`, `qc_decisions.csv`, two-pass AneuFinder,
   `CN_review`) is kept.
6. Copy number stays with AneuFinder. Copy-number results are added to the RNA layer
   of the same cells.

### 3.2 Open items (do not decide; ask the user)

1. How CellenONE names well positions for 384-well and 5,184-well targets, and how
   they map to `well_map.tsv` (`cellenone_pos`).
2. Sequencing run configuration: instrument and read lengths. Only the minimum
   requirements in 2.4 are known.
3. `well_id` format. Proposal: `R<row>C<col>`, two-digit zero-padded, for example
   `R01C01`.
4. Mappability reference and blacklist, to be rebuilt from euploid cells sequenced
   with wellDR-seq.
5. Re-tuned read cutoffs: `usable_reads_pass_cutoff`, `usable_reads_warn_cutoff`,
   `min_reads_for_model`.
6. AneuFinder bin size (currently 1 Mb; the authors used about 220 kb). A change needs
   a matching GC template.

## 4. Change inventory

Paths are relative to the repo root. "Rule" means a Snakefile rule.

| Item | Status | Reason |
|------|--------|--------|
| rule `demultiplex`, `workflow/scripts/preprocessing/extract_umi_barcode.py` | Remove | DNA wells are demultiplexed by bcl-convert on the index reads |
| rule `filter_dimers`, `workflow/scripts/preprocessing/filter_adapter_dimers.py` | Remove | Specific to the Y-adapter ligation library; replaced by adapter trimming |
| rule `deduplicate_bam` (umi_tools), rule `summarize_dedup`, `workflow/scripts/preprocessing/summarize_umi_tools_dedup.py` | Remove | The DNA library has no UMI; replaced by position-based marking |
| 384 subplate machinery: subplate formula in `workflow/scripts/reporting/plate384_layout.py`, `workflow/scripts/analysis/link_subplate_models.py`, `aneufinder.plate384_models`, `SUB_BASE` / `LOG_TAG` / `wid()` subplate logic, SL1-SL4 tabs in the review, `PLATE_FORMAT` | Remove | One plate is one library; geometry comes from `well_map.tsv` |
| `preprocessing:` block in `config.yaml` (UMI, barcode, trim offsets, dimer adapter) | Remove | Legacy read structure |
| rule `check_dna_inputs` | New | Compare `dna/fastq/` against `well_map.tsv`: missing, empty, unexpected wells |
| rule `trim_dna` | New | Nextera adapter trimming of short inserts (fastp proposed) |
| rule `markdup` | New | Position-based duplicate marking (`samtools markdup`; compare with `sambamba markdup` in Phase 2) |
| rule `dna_complexity` | New | Duplication rate and estimated library size per well |
| rule `index_hopping_qc` | New | Reads in empty wells summarised by row and by column |
| rule `rna_demux` | New | Read 1 parsing: adapter and polyT check, RNA-CB1 and RNA-CB2 with up to 1 mismatch each |
| rule `trim_rna` | New | Adapters, polyA, Nextera N7 and S5 (section 2.2.7) |
| rule `starsolo` | New | SmartSeq mode, `--soloUMIdedup Exact NoDedup`, Gene and GeneFull |
| rule `rna_qc` | New | Genes, counts, mitochondrial fraction per well |
| rule `gene_bin_map` | New | Assign each GTF gene to an AneuFinder bin |
| rule `cn_to_genes` | New | Per-cell integer copy number per gene from second-pass models |
| rule `build_object` | New | Combined RNA + CN object joined on `well_id` |
| `workflow/scripts/make_samplesheet.py` | New | bcl-convert sample sheet for the DNA library from `well_map.tsv` |
| `workflow/envs/rna.yaml` | New | STAR, Trimmomatic or fastp, R packages for integration |
| STAR hg38 index | New | Built once, shared under `resources/reference/` or a configured path |
| `well_map.tsv` schema and validator | New | Single source of plate geometry and barcodes |
| rule `align` | Change | Input from `trimmed/`; read groups simplified to plate and well |
| `well` wildcard | Change | `W\d+` replaced by the `well_id` format from `well_map.tsv` |
| `workflow/scripts/reporting/generate_qc_review.py` and `assets/` | Change | Generic rows x cols grid, duplication and complexity metrics, RNA panel, 5,184-well scale |
| `workflow/scripts/reporting/ingest_cellenone.py` | Change | Positions from `well_map.tsv` (`cellenone_pos`), pending open item 3.2.1 |
| `workflow/scripts/reporting/render_well_profiles.R`, rule `cn_review` | Change | 5,184-well scale; keep the heatmap height cap |
| rule `multiqc`, `multiqc_config.yaml` | Change | Add fastp, markdup and STARsolo logs; drop umi_tools |
| mappability reference and blacklist inputs | Change | Rebuilt from wellDR-seq euploid cells (open item 3.2.4) |
| `usable_reads`, `min_reads_for_model` | Change | Definition (non-duplicate, MAPQ >= 1, proposal) and cutoffs (open item 3.2.5) |
| `submit_pipeline.sh`, `pre_check.sh` | Change | New modes; remove `PLATE_FORMAT` and 384 checks |
| `workflow/scripts/analysis/run_aneufinder.R` | Keep | Do not simplify. See "Current AneuFinder / GC / Blacklist State" in `AGENT_CONTEXT.md` (GC object class, `GC` column name, bin naming) |
| rule `check_gc_rds` | Keep | Leave as is; it is currently dead in the DAG (`docs/01_METHODS.md` B5) and must never own the GC RDS |
| rule `generate_blacklist`, `get_mappability_bam` | Keep | Inputs change (above), logic stays |
| two-pass flow: `qc_review`, `validate_qc_decisions`, `derive_included_wells`, `run_aneufinder_reviewed`, `render_reviewed_profiles` | Keep | Decision 3.1.5 |
| rule `fastqc` | Keep | Input path changes only |
| CellenONE image calls (`render_cellenone_images.py`) | Keep | Display-only layer, unchanged |

## 5. Phases

Each phase ends with a hard stop. Show the user the listed figures first, with tables
as backing data, and wait for approval before starting the next phase.

### Phase 0: protocol specification and plate map

**Scope**

1. Write `docs/PROTOCOL_wellDR.md`: library structures, primer sequences, read
   layouts, run requirements, all with sources (section 2).
2. Verify the axis assignment of i7, i5, RNA-CB1 and RNA-CB2 against the paper's
   Table S5 and the author repository. Record the result.
3. Define the `well_map.tsv` schema and write a validator (unique `well_id`; unique
   `(dna_i7, dna_i5)` and `(rna_cb1, rna_cb2)` pairs; barcode lengths; minimum
   Hamming distance within each barcode set; row and column ranges).
4. Write `workflow/scripts/make_samplesheet.py` (bcl-convert v2 sample sheet, i5
   orientation as an explicit argument, `Sample_ID` = `well_id`).
5. Prepare a `well_map.tsv` for the pilot plate with the user. Leave `cellenone_pos`
   empty if open item 3.2.1 is still open.

**Files touched:** `docs/PROTOCOL_wellDR.md`, `workflow/scripts/make_samplesheet.py`,
`workflow/scripts/validate_well_map.py` (name proposed), tests under `tests/` or
`ad_hoc_checks/`.

**Checks:** validator passes on the pilot map and fails on deliberately broken maps;
the sample sheet round-trips to the same barcodes.

**Show the user:** the pilot plate map coloured by barcode, the Hamming distance
distribution of each barcode set, the generated sample sheet.

**Stop for review.**

### Phase 1: remove legacy-only code; DNA input layer

**Scope**

1. Remove the items marked Remove in section 4.
2. Add `well_map.tsv` loading to the Snakefile; replace the well wildcard.
3. Add rule `check_dna_inputs` and point `fastqc` at `dna/fastq/`.

**Files touched:** `Snakefile`, `config.yaml`, `submit_pipeline.sh`, `pre_check.sh`,
the removed scripts, `workflow/scripts/reporting/plate384_layout.py`.

**Checks**

1. `snakemake -n` resolves a DAG for the pilot plate.
2. Reads per well from the DNA FASTQs, on the plate map.
3. Undetermined fraction of the DNA run (from bcl-convert reports).
4. Reads in wells marked empty.

**Show the user:** reads-per-well plate map; Undetermined fraction; empty-well reads.

**Stop for review.**

### Phase 2: trimming, alignment, duplicate marking, complexity, index hopping

**Scope:** rules `trim_dna`, `align` (changed), `markdup`, `dna_complexity`,
`index_hopping_qc`; MultiQC updated.

**Checks**

1. Duplication rate against sequencing depth per well.
2. Agreement of `samtools markdup` and `sambamba markdup` on a subset of wells
   (duplicate fraction per well, and read-level agreement).
3. Index hopping: reads in empty wells against the row and column totals of the
   wells sharing their i7 or i5.

**Show the user:** duplication against depth scatter; markdup tool agreement plot;
index hopping heatmap by row and column.

**Stop for review.**

### Phase 3: AneuFinder on wellDR-seq data and the review report

**Scope**

1. Run AneuFinder pass 1 on `bam/<well_id>.bam`. Confirm which duplicate and MAPQ
   filters apply inside AneuFinder (`run_aneufinder.R` does not set
   `remove.duplicate` or `min.mapq`, so AneuFinder defaults apply; verify).
2. Rebuild the mappability reference BAM and blacklist from wellDR-seq euploid cells
   (open item 3.2.4).
3. Check the GC template against the chosen bin size (open item 3.2.6).
4. Propose re-tuned cutoffs from the data (open item 3.2.5); the user decides.
5. Review report on the generic grid, with duplication and complexity metrics, at
   384 wells and at a simulated 5,184-well scale.

**Checks:** profiles of known euploid and aneuploid cells; empty-bin fraction per
cell; the review report renders and stays responsive at 5,184 wells.

**Show the user:** profiles of control cells; reads against empty-bin fraction with
proposed cutoffs; screenshots of the review report at both scales.

**Stop for review.**

### Phase 4: RNA branch

**Scope:** rules `rna_demux`, `trim_rna`, `starsolo`, `rna_qc`; `workflow/envs/rna.yaml`;
STAR index; RNA panel in the review report.

**Checks**

1. Read 1 structure: fraction of reads with the fixed adapter at bases 18-31 and
   polyT from base 40; barcode match rates at 0 and 1 mismatch.
2. Genes and counts per well; mitochondrial fraction.
3. DNA and RNA recovery per well on one plate map.

**Show the user:** Read 1 structure summary; genes per well plate map; DNA reads
against RNA genes per well.

**Stop for review.**

### Phase 5: integration

**Scope:** rules `gene_bin_map`, `cn_to_genes`, `build_object`; analyses in section
2.3.

**Checks**

1. Per-gene CN matrix sanity: genes per bin, agreement of gene CN with bin CN, chrX
   and chrY behaviour.
2. Gene dosage correlations between subclone pseudobulk expression and copy number.
3. RNA-inferred CNV (CopyKAT, inferCNV) against DNA calls.

**Show the user:** subclone UMAP; gene dosage plot; RNA-inferred against DNA copy
number heatmaps.

**Stop for review.**

## 6. Principles

1. Work additively within the branch where possible; remove legacy code only where
   section 4 says so.
2. Document output contracts in `README.md` in the same change that implements them.
3. Keep generated data, runtime state and the paper PDF out of git.
4. One phase per session. Do not start the next phase without user approval.
5. Do not fill in open items (section 3.2) without the user.
6. Do not simplify `run_aneufinder.R` or the GC template handling (see
   `AGENT_CONTEXT.md`).
7. Commit to `welldr` only. Do not edit branch `dlp+`.
