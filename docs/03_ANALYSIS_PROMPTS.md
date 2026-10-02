# Analysis prompts — outstanding and completed

Part of the scNLAIII-seq document set. Read `00_START_HERE.md` and `02_FINDINGS.md` first — several
hypotheses below are already closed, and re-deriving them is the commonest way to waste a cycle.

Paste one prompt at a time into a session rooted at
`/nemo/project/proj-tracerX/working/VCAM1_GnT/TOOLS/scDNA_pipeline`.

**Conventions assumed by all prompts.** Data root `/nemo/project/proj-tracerX/working/VCAM1_GnT/DATA/384_well/`, plates `plate17_umi_test`, `plate21`, `plate22`, each holding `<plate>_1..4` subplate dirs with `demux/`, `filtered/`, `raw_bam/`, `bam/`, `dedup/`, `aneufinder/`. Reference `…/DATA/ref/GRCh38_canonical.fa`. Read structure: R1 = `[8nt well barcode][14nt linker incl. CATG]`; R2 = `[12nt UMI][8nt well barcode][14nt linker incl. CATG]`; `r1_trim_5prime: 22`, `r2_trim_5prime: 34` (see `config.yaml`). Write new scripts under `ad_hoc_checks/<name>/`, never into `workflow/`. Run R from `/tmp` with the native module `R/4.4.1`, **not** from the repo root — the repo `.Rprofile` activates an incomplete renv and hides system ggplot2 (see `ad_hoc_checks/plate_qc_summary/README.md`).

**Status, 29 July 2026.** A1, A2 and A4 are **done** — see the READMEs in `ad_hoc_checks/nlaIII_digest_check/`, `ad_hoc_checks/insert_size_vs_spri/` and `ad_hoc_checks/background_floor/`. Their prompts are kept below for the record. **A6 and the two new prompts, A11 and A12, are now the highest-value next steps.** A11 measures the realised gain from 3′ adapter trimming (the biggest single win available); A12 re-demuxes with `BC1 == BC2` enforced, which A4 predicts removes 33 % of contaminating molecules against 2.3 % of genuine ones.

A9 is the one that directly produces a defensible QC threshold, and it should now include percent-at-CATG as a second axis.

---

## A1 — ✅ DONE — Is the NlaIII digestion actually complete?

> **Result:** digestion is complete and ligation is specific (96–97 % of PASS read ends at a genuine CATG; 2–3 % of inserts skip an uncut site; zero soft-clips in 1.5 bn read ends). Two new findings: FAIL wells are library-like in 91 % of cases on plate17 but only 31 % / 23 % on plate21 / 22; and 40 % of the digest is below the 79 bp read length and structurally unalignable.

```
In this scNLAIII-seq pipeline, genomic DNA is digested with NlaIII (recognition site CATG,
leaving a 4-nt 3' CATG overhang) and a forked Y-adapter is ligated to the resulting ends.
If digestion were complete and ligation specific, essentially every sequenced fragment should
begin exactly at a genuine CATG site in hg38.

Measure how true that is, for plate21, plate22 and plate17_umi_test.

Method:
1. For each subplate, take the aligned BAMs. Use the RAW bams (raw_bam/) as well as the
   deduplicated ones (bam/) so we can see whether dedup changes the picture.
2. For each properly-paired primary alignment, take the 5' end of the read in genomic
   orientation (leftmost coordinate for forward reads, rightmost for reverse) and, allowing
   for the 4 nt of CATG that the trimming step removes, ask whether the reference sequence
   immediately 5' of that position is CATG. Be explicit in the code about the exact offset
   you use and justify it against config.yaml's r1_trim_5prime: 22 / r2_trim_5prime: 34 and
   the 14-nt linker that ends in CATG. Sanity-check the offset on a handful of reads by
   printing the reference context.
3. Report, per subplate and per well:
     - % of read starts at a genuine CATG
     - % at CATG within +/-1, +/-2 bp (to absorb trimming/soft-clip slop)
     - the distribution of offsets from the nearest CATG
4. Separately, build the in-silico NlaIII digest of GRCh38_canonical.fa: every CATG position,
   and the resulting fragment length distribution. Report the number of sites and fragments.
5. Compare observed capture against that in-silico map: what fraction of all NlaIII sites in
   the genome is hit at least once (a) in the best well, (b) in a median PASS well, (c) in a
   median FAIL well, (d) pooled across a whole subplate.

Then answer directly:
- Is digestion complete, or are we losing ends?
- Are the reads in FAIL wells (usable_reads < 50000) at CATG sites at the same rate as reads
  in PASS wells (>= 100000)? If FAIL-well reads are LESS often at CATG, their reads are not
  genuine library and that changes what "FAIL" means.
- What is the empirically usable fragment ceiling per cell, i.e. how many NlaIII fragments
  are in a size range we actually observe?

Write results to ad_hoc_checks/nlaIII_digest_check/ with a short README.md stating the
conclusion in plain language, plus the tidy tables and plots.
```

---

## A2 — ✅ DONE — Is the 0.8× SPRI cleanup throwing away the library?

> **Result:** no. Short fragments are *enriched* (O/E 1.4–2.2 renormalised), and the depletion runs against long fragments. More decisively, a cleanup performed on the pool cannot generate between-well differences at all. But the ~188 bp dimers were never removed: plate21 FAIL wells are a median 25.6 % dimer and 111 of 339 are over 50 %. Recovering the sub-100 bp classes is worth ×2.1–2.4 in molecules per well.

```
Hypothesis to test: the first 0.8x SPRI cleanup (performed on the POOLED ligation, before any
amplification, when the molecules ARE the library complexity) is size-selecting away the
genuine short NlaIII fragments while retaining the ~130 bp adapter dimers, and is losing a
large fraction of molecules outright at picogram input.

Using plate21, plate22 and plate17_umi_test:

1. From the BAMs, extract the insert size (TLEN) distribution per well and per subplate.
   Use properly-paired primary alignments only. Report the full distribution, not just the
   mean - I want the density from 0 to 1000 bp with the low tail resolved.
2. Overlay that on the in-silico NlaIII fragment length distribution of GRCh38_canonical.fa
   (build it if A1 has not already). Normalise both to density so they are comparable.
3. Quantify the depletion: for size bins 0-100, 100-150, 150-200, 200-300, 300-500, 500-1000,
   >1000 bp, report observed/expected ratio.
4. Ask whether the depletion profile differs between PASS wells (usable_reads >= 100000) and
   FAIL wells (< 50000). If PASS wells simply have MORE of the same size distribution, the
   size selection is uniform across wells and cannot explain the between-well inequality.
   If FAIL wells are additionally depleted of short fragments, the cleanup interacts with
   per-well yield and is a bigger problem than it looks.
5. Estimate: if we recovered the fragments in the 0-200 bp range at the same efficiency as
   the 200-500 bp range, how many more unique molecules per well would we have? Give a
   number with its assumptions stated.

Deliverable: ad_hoc_checks/insert_size_vs_spri/ with a README.md that answers, in one
sentence, whether the bead cleanup is a first-order complexity killer here.
```

---

## A3 — Are the failing wells exhausted, or just under-sequenced?

```
For every well of plate21, plate22 and plate17_umi_test, estimate library complexity and
saturation properly, so we can state definitively whether more sequencing would rescue any
wells.

1. Install/use preseq (c_curve and lc_extrap) on the per-well raw BAMs. Where preseq fails on
   low-count wells, fall back to a Good-Toulmin or Chao1-style estimator and say which was
   used for each well.
2. Per well report: observed unique molecules, estimated total library size, the fraction of
   the library already observed, and the extra reads needed to reach 100,000 unique molecules
   (report "unreachable" where the estimated library size is below 100,000).
3. Produce a saturation curve panel for: the 5 best wells, 5 median PASS wells, 5 WARN wells,
   and 20 random FAIL wells.
4. Aggregate to the subplate and plate level: what is the total estimated library complexity
   per subplate, and how much of it have we already sequenced?
5. Directly answer: how many wells could cross the 100k threshold with more reads, and at
   what cost? Express it as reads-per-rescued-well.

Also do a separate check on deduplication itself, because the previous benchmark
(ad_hoc_checks/plate17_umi_tools_benchmark.qmd) was never actually run and contains no
results:
6. Re-run dedup on a sample of wells with umi_tools directional, unique, and cluster, and
   with position-only dedup (no UMI). Compare unique_reads across methods. Because NlaIII
   fragments all start at fixed CATG sites, positions are heavily reused; report
   mean_umi_per_pos and max_umi_per_pos and state whether the 12-nt UMI is anywhere near
   saturation at any position. Conclude whether directional dedup is over-collapsing.

Deliverable: ad_hoc_checks/complexity_saturation/ with README.md giving a one-line answer to
"would more reads help?" and a one-line answer to "is dedup over-collapsing?".
```

---

## A4 — ✅ DONE — What actually is the ~3,000-read floor in dead wells?

> **Result:** both hopping and chimera, separably. 27.6 % of FAIL-well molecules are demonstrably another well's — ~17 pp index hopping, ~11 pp PCR chimera. 18 % of FAIL wells are empty; 54 % hold a real but minuscule library, so the headline is "library made and lost". Because the UMI is attached during that well's own ligation, shared exact UMIs can only post-date ligation, which excludes ambient DNA and the dispense step. Also found: the demux assigns on BC1 alone and discards BC2 unused.

```
Plates 17-22 have NO cell-free control wells, so we have never calibrated the background.
Failing wells still receive a median of ~50,000 raw reads and ~3,300 unique molecules. I need
to know whether that floor is (a) a genuine but tiny single-cell library, (b) barcode hopping
/ index misassignment from the handful of high-yield wells, or (c) adapter-mediated chimeras
formed during the 26-28 cycle PCR.

For plate21 and plate22:

1. Barcode integrity. The demux requires the well barcode in BOTH R1 and R2. Recover, from
   the raw FASTQs, the reads that were DISCARDED because R1 and R2 barcodes disagreed. Report
   the rate, and build the R1-barcode x R2-barcode confusion matrix. A hopping/chimera process
   produces an asymmetric matrix concentrated on the high-yield wells as donors; sequencing
   error produces a symmetric matrix concentrated on Hamming-distance-1 pairs. Say which it
   looks like.
2. Molecule sharing. For each FAIL well, compute what fraction of its unique molecules
   (defined as genomic position + strand, ignoring UMI) is also present in the top-5 wells of
   the same subplate. Compare against the null expectation given the number of NlaIII sites
   and the number of molecules involved - compute that null explicitly, do not eyeball it.
   Genuine independent single cells from the same tissue will share sites at the null rate;
   hopping-derived reads will share far above it.
3. UMI sharing. Repeat (2) but requiring position + strand + exact UMI. Sharing an exact
   12-nt UMI at the same position is essentially impossible between two genuinely independent
   molecules, so this is the sharp test.
4. Report a per-well "estimated contamination fraction" and re-derive each well's
   contamination-corrected usable_reads.
5. State how many currently-FAIL wells have essentially zero genuine molecules once
   contamination is subtracted.

Deliverable: ad_hoc_checks/background_floor/ with README.md, a per-well contamination table,
and a clear statement of what fraction of the "FAIL" population is empty vs genuinely
low-yield. This determines whether the problem is "no library made" or "library made and lost".
```

---

## A5 — Are the wells that survive a biased sample of cells?

```
Concern: plate23's qPCR shows cell diameter strongly predicts whether a well yields
amplifiable material (AUC 0.72-0.75, p ~ 1e-17; 8% positive at <=13um rising to 56% at >19um).
If cell size tracks DNA content, then the wells that survive to copy-number calling are
enriched for S-phase and polyploid cells - exactly the cells that produce spurious CN calls.
For a project about copy number in NORMAL kidney, that is a threat to the biological
conclusion, not just to yield.

For the wells that passed on plate17, plate21 and plate22:

1. From the AneuFinder models in aneufinder/MODELS/method-edivisive/*.RData, extract per-well:
   inferred ploidy / mean copy-number state, the fraction of the genome not at the modal
   state, segment count, and the AneuFinder quality metrics (spikiness, bhattacharyya,
   entropy, num.segments).
2. Derive an S-phase indicator per well independent of AneuFinder: within-genome coverage
   variance in large bins after GC correction, and replication-timing correlation (correlate
   per-bin coverage against a public hg38 Repli-seq / replication timing track; fetch one if
   not present and record the source). S-phase cells show a characteristic positive
   correlation with early replication timing.
3. Ask whether ploidy, S-phase score, segment count or genome-fraction-altered correlate with
   usable_reads across PASS wells. Use Spearman and report effect sizes, not just p-values.
4. For plates 21 and 22, join to the CellenONE isolated cell diameter
   (cellenone/wells_raw.tsv) and test whether diameter predicts ploidy or S-phase score among
   the PASS wells specifically. This is the direct test of the selection-bias worry.
5. Report the fraction of PASS wells that would be excluded by a reasonable S-phase / ploidy
   filter, and what the CN conclusions look like with and without them.

Deliverable: ad_hoc_checks/selection_bias/ with README.md stating whether the surviving cells
are a biased sample and, if so, how much of the current CN signal survives correction.
```

---

## A6 — Full adapter dimer and free-adapter accounting

```
The optimisation history shows adapter dimers dominating early plates (>95%) and the current
0.08 nM adapter concentration was chosen to suppress them. I want a complete accounting of
what the sequenced molecules actually are, per well, for plate21, plate22 and plate17_umi_test.

Working from the RAW FASTQs (before the pipeline's filter_dimers step):

1. Classify every read pair into: (i) genuine insert - adapter/linker present at the expected
   5' position and genomic sequence after it; (ii) adapter dimer - linker sequence
   (CAGTCAGCGT and its reverse complement) present again downstream, i.e. adapter ligated to
   adapter; (iii) adapter-only / no insert; (iv) linker present but insert too short to align;
   (v) unclassifiable. Report counts and percentages per well.
2. Plot dimer fraction against usable_reads per well (log scale) for all 1152 wells across the
   three plates, coloured by plate. Fit and report the relationship. We already know from
   well_metrics.csv that dimer rate rises from ~2% in the top yield decile to ~60% in the
   bottom two deciles - I want that confirmed from raw reads and broken down by dimer subtype.
3. Critically: distinguish "this well made dimers instead of library" from "this well made
   nothing and dimers are all that is left". Compute absolute dimer counts per well, not just
   the fraction. If absolute dimer counts are CONSTANT across wells while genuine inserts vary
   160-fold, then dimer formation is a fixed background and the variance is entirely in the
   genuine-ligation channel - which would strongly support the ligation-efficiency model.
   If absolute dimer counts are HIGHER in failing wells, adapter is being consumed by
   self-ligation in those wells and that is a different (and more fixable) problem.
4. Estimate, per well, the implied number of adapter molecules consumed by dimers versus by
   genuine ligation, and compare against the 0.4 uL x 2 nM = ~4.8e8 adapter molecules dispensed.

Deliverable: ad_hoc_checks/dimer_accounting/ with README.md answering point 3 explicitly -
it is the discriminating result.
```

---

## A7 — How many PCR founder molecules does each well actually have?

```
The protocol runs 26-28 PCR cycles in four blocks with a bead cleanup between each, from a
pooled sub-picogram library. I want to know whether that produced jackpotting, and to
estimate the number of independent founder molecules per well.

For plate21, plate22 and plate17_umi_test:

1. From the raw (pre-dedup) BAMs, build the duplicate family size distribution per well: for
   each unique molecule (position + strand + UMI, using the same grouping umi_tools applies),
   how many reads support it? Report the full distribution per well, not summary statistics.
2. Fit the family-size distribution. Under clean PCR, family sizes are roughly geometric /
   negative-binomial around the mean amplification factor. Jackpotting produces a heavy right
   tail - a small number of enormous families. Quantify the tail: what fraction of all reads
   sits in the top 1% of families? Report per well and per subplate.
3. Compare PASS wells against FAIL wells. If FAIL wells have the SAME mean family size but
   fewer families, they had fewer founder molecules and the loss is upstream of PCR (ligation
   or cleanup). If FAIL wells have SMALLER families, they were out-competed during the pooled
   PCR and the loss is at amplification. This is the key discriminating result - state it
   plainly.
4. Estimate founder molecules per well = number of distinct families, and cross-check against
   the preseq estimate from A3.
5. Test for a PCR-block signature: do the family sizes cluster around powers of 2 consistent
   with amplification starting at different cycles?

Deliverable: ad_hoc_checks/pcr_founders/ with README.md answering point 3 in one sentence.
```

---

## A8 — Model the plate-position effect properly

```
There is a real but secondary spatial gradient in well yield. In plates 21+22 pooled,
Kruskal-Wallis across rows gives p = 0.0039 and across columns p = 4.2e-7; Spearman of well
number against usable reads is rho = -0.14, p = 8e-5. Separately, plate23's qPCR shows rows
A-H at 20.8% positive versus rows I-P at 36.5%, and that row effect is INDEPENDENT of cell
diameter (row index vs diameter rho = -0.013). I want this characterised properly rather than
as a marginal test.

1. Map every well of plate17, plate21, plate22 to its true physical 384-well position
   (subplate interleaving - see the Snakefile's plate layout logic, and plate_layout_tsv in
   config.yaml). Do NOT assume sequential W01..W96 corresponds to A01..H12; verify against
   barcodes.tsv and the layout formula, and state which mapping you used.
2. Fit a model of log10(usable_reads + 1) with terms for: row, column, edge-vs-interior,
   distance from plate centre, subplate, and plate. Use a mixed model with subplate as a
   random effect. Report variance explained by each term.
3. Test specifically for an evaporation signature: edge wells (row A/P, column 1/24) versus
   interior. Report the effect size in fold-change of usable reads.
4. Test for a dispensing-order signature: if the CellenONE dispense order is recoverable from
   the .Run folder (DropNo in Reordered_*_isolated.xls), test whether yield correlates with
   dispense order independently of physical position. This separates "how long the well sat
   open" from "where the well is in the thermal block".
5. Repeat (2)-(4) on the plate23 qPCR positivity outcome so we have the same model on both
   readouts, and report whether the spatial patterns agree between the two plates.

Deliverable: ad_hoc_checks/spatial_model/ with README.md quantifying how much of the total
variance in yield the spatial terms explain. If it is under ~10%, we deprioritise it; if it is
over ~25%, evaporation/thermal control becomes a primary target.
```

---

## A9 — Replace the 50k/100k thresholds with empirically derived ones

```
The current QC gate (usable_reads >= 100000 PASS, 50000-99999 WARN, < 50000 FAIL) is inherited
rather than derived. Derive it from the data.

Using the PASS wells across plate17, plate21 and plate22 (the wells with enough depth to
subsample from):

1. For each such well, subsample the deduplicated BAM to a ladder of unique-molecule counts:
   10k, 20k, 30k, 50k, 75k, 100k, 150k, 200k, 300k, 500k, and full depth. Do 5 independent
   subsampling replicates at each level, with a fixed seed recorded.
2. Run AneuFinder (edivisive, 1 Mb bins, the same config the pipeline uses) on every
   subsample. Reuse run_aneufinder.R unchanged - stage subsampled BAMs into a temp input dir
   exactly as stage_aneufinder_input.py does. Do NOT modify run_aneufinder.R.
3. For each subsample, compare the CN profile against that well's own full-depth profile as
   ground truth. Report: per-bin CN concordance, breakpoint concordance (within a tolerance
   you state), and the AneuFinder quality metrics (spikiness, bhattacharyya, num.segments).
4. Plot concordance against unique-molecule count, with a curve per well and a pooled loess.
   Identify the knee: the depth at which concordance stops improving materially, and the depth
   below which between-replicate variability explodes.
5. Do this separately for 1 Mb and 500 kb bins so we know what depth each resolution needs.
6. Report the recommended PASS and WARN thresholds with a stated tolerance (e.g. "the depth at
   which 95% of replicates reproduce the full-depth call in >=95% of bins"), and say how the
   current 100k/50k cutoffs compare.

Deliverable: ad_hoc_checks/depth_vs_cn_quality/ with README.md giving the recommended
thresholds, the evidence, and a proposed patch to config.yaml's qc_review cutoffs. Do NOT
change config.yaml yourself - just state the proposed values.
```

---

## A11 — NEW — measure the gain from 3′ adapter trimming

```
A1 and A2 established that there is no 3' adapter trimming anywhere in this pipeline and that
bowtie2 runs end-to-end (bowtie2_params: "--very-sensitive"). R1 is 79 bp after r1_trim_5prime: 22,
so any NlaIII fragment shorter than the read runs into the adapter and cannot align. That makes
5,585,477 fragments - 40.4% of the digest - structurally invisible, and the counterfactual in
A2 predicts x2.39 / x2.35 / x1.92 more unique molecules per median PASS well on plate21 / plate22 /
plate17 if they were recovered.

Test that prediction on real data before anyone commits to re-aligning 1152 BAMs.

1. Pick ONE subplate with a good spread of well yields - plate21_1 is the obvious choice
   (12 PASS wells, 8.4 M unique molecules).
2. Re-run preprocessing for that subplate only, into a SEPARATE output directory, with 3' adapter
   trimming added. Do it two ways so we can compare:
     (a) cutadapt / fastp trimming of the read-through adapter before alignment, keeping bowtie2
         end-to-end. State the exact adapter sequence you trim and justify it from the read
         structure (R1 = 8 nt barcode + 14 nt linker + insert; the read-through is the reverse
         complement of the opposite adapter arm).
     (b) bowtie2 --local instead of end-to-end, no trimming.
   Do NOT modify the production plate directories or the Snakefile - build this as a
   side-by-side run under ad_hoc_checks/.
3. Also lower or remove bowtie2's -X 500 in one arm and report what it changes.
4. For each arm, recompute per well: total aligned reads, unique molecules after umi_tools dedup,
   distinct NlaIII fragments recovered as complete pairs (reuse the scan from
   ad_hoc_checks/nlaIII_digest_check/scripts/scan_bams.py), and the insert-size distribution.
5. Report the REALISED fold change in unique molecules per well against the predicted x2.39, and
   break it down by well gate status. State clearly whether the rescued short fragments map
   uniquely (report their MAPQ distribution) - the concern is that after trimming they carry too
   little unique sequence.
6. Report the cost: wall-clock and CPU-hours to re-run one subplate, extrapolated to 1152 wells.
7. Check the rescued fragments are not GC- or repeat-biased in a way that would distort CNV
   calling: GC distribution and blacklist overlap of the newly recovered fragments vs the
   previously recovered ones.

Deliverable: ad_hoc_checks/trim_benchmark/ with README.md leading on the realised fold change and
a clear recommendation on whether to re-run the whole project.
```

---

## A12 — NEW — re-demux with BC1 == BC2 enforced

```
A4 established that workflow/scripts/preprocessing/extract_umi_barcode.py assigns a read to a
well on BC1 ALONE (line ~181): BC2 is parsed, written into the read header, and never used. So
no read has ever been discarded for barcode disagreement, and every barcode-discordant read is
already inside demux/, inside the BAMs, and inside usable_reads.

A4 measured what enforcing BC1 == BC2 would do: it removes 33% of a FAIL well's contaminating
molecules and only 2.3% of its genuine ones (43x directional; PASS wells 28% vs 0.59%). It
reaches the PCR-chimera route only - an index-hopped read carries both inline barcodes intact.

Implement and validate it. Do NOT change the production pipeline until the numbers are in.

1. Add an optional strict mode to extract_umi_barcode.py that requires BC1 == BC2 (allowing a
   configurable mismatch tolerance of 0, and also test 1). Keep the current behaviour as the
   default and gate the new behaviour behind a config flag, so no existing output changes.
2. Re-demux ONE sub-plate of plate21 and ONE of plate22 both ways, into separate directories,
   and run them through to deduplicated BAMs.
3. Report per well: reads assigned, reads dropped by the strict filter, unique molecules,
   and the gate status under each mode. Quantify how many wells change PASS/WARN/FAIL class.
4. Cross-check against A4's per-well contamination table
   (ad_hoc_checks/background_floor/tables/per_well_contamination.csv): does the strict filter
   preferentially remove the molecules A4 independently flagged as belonging to another well?
   Report the overlap and the realised removal rate against A4's predicted 33% / 2.3%.
5. Re-run the A1 percent-at-CATG metric on the strict-mode BAMs. If the filter is removing
   contamination rather than real library, the CATG rate of FAIL wells should RISE.
6. Report the cost in reads: A4 measured that strict demux would have dropped 4.50-8.18% of
   assigned reads per sub-plate. Confirm that and say where the loss falls (PASS vs FAIL wells).

Deliverable: ad_hoc_checks/bc2_strict_demux/ with README.md giving a clear recommendation on
whether to make strict mode the default, and what it costs.
```

---

## A10 — After plate23 is sequenced (run this once, not before)

```
plate23 was assayed by per-well qPCR BEFORE library prep and has now been sequenced. This is
the only plate where we can link a per-well pre-sequencing measurement to a per-well read
yield, so treat it as the key experiment.

Inputs: the plate23 sequencing output processed through the standard 384 pipeline, plus
"plate23 - qPCR QC/plate23_well_classification.csv" (columns: well, amp_status, cq, tm1,
classification in {pos, amb, neg}, and the CellenONE metrics diam, elong, circ, inten, red,
orange, blue).

1. Join qPCR call to per-well usable_reads. Report the read-yield distribution for pos, amb
   and neg wells separately, and the PASS rate within each.
2. Compute how well the qPCR call predicts sequencing success: AUC of Cq against usable_reads,
   and the PASS rate in qPCR-positive versus qPCR-negative wells. Report the positive and
   negative predictive value of the qPCR call at the 100k threshold.
3. THE CRITICAL QUESTION. Among qPCR-POSITIVE wells - wells we know made amplifiable material -
   what fraction still fails to reach 100,000 usable reads? That number is the size of the
   loss that happens AFTER ligation, i.e. in pooling, bead cleanup and PCR. Everything else in
   this investigation hinges on it. State it prominently.
4. Conversely, among qPCR-NEGATIVE wells, how many nonetheless sequenced well? A non-trivial
   number would mean the qPCR is insensitive (consistent with it running at its single-copy
   detection limit around Cq 35) rather than the wells being empty.
5. Re-test cell diameter against usable_reads on this plate, and against usable_reads
   CONDITIONAL on qPCR-positive. If diameter predicts qPCR positivity but not read yield among
   positives, then cell size acts entirely at the ligation step - which is what the
   concentration-limited-ligation model predicts.
6. Decompose the overall success rate into the two sequential bottlenecks:
   P(sequencing success) = P(qPCR positive) x P(PASS | qPCR positive). Report both terms with
   confidence intervals.

Deliverable: ad_hoc_checks/plate23_qpcr_vs_seq/ with README.md leading on the answer to
point 3.
```
