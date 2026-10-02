# Findings — what is established, what is open

Read `00_START_HERE.md` first. This is the knowledge base: every result that has been established, with its evidence, its source and its confidence. **Read §3 before proposing any analysis** — eight hypotheses have already been rejected on data.

**Last updated:** 29 July 2026, after A1 (`ad_hoc_checks/nlaIII_digest_check/`), A2 (`ad_hoc_checks/insert_size_vs_spri/`) and A4 (`ad_hoc_checks/background_floor/`). Each has its own README with full method and caveats; this document summarises and cross-links them.

---

## 1. The headline

Only **8–18 % of wells** reach the 100,000-usable-read gate. The cause is a **molecular-complexity ceiling fixed before the library reaches the sequencer**, distributed extremely unequally between wells.

| plate | dispensing | PASS (≥100k) | aligned reads | **total unique molecules** | dup rate | median unique/well | top-13 of 96 share | Gini |
|---|---|---|---|---|---|---|---|---|
| plate17 | FACS | **69 / 384 (18.0 %)** | 903 M | **29.7 M** | 96.7 % | 12,341 | 64 % | 0.78 |
| plate21 | CellenONE | 41 / 384 (10.7 %) | 341 M | **31.9 M** | 90.7 % | 3,676 | 92 % | 0.90 |
| plate22 | CellenONE | 32 / 384 (8.3 %) | 365 M | **18.4 M** | 95.0 % | 4,019 | 90 % | 0.89 |

Read that row-wise: **plate17 was sequenced 2.6× deeper than plate21 and ended up with less total complexity.** The ceiling is ~20–32 M molecules per plate and is conserved across different patients, sorters and months.

**Operational consequence:** every well is already sequenced ~7× deeper than its own complexity. More reads cannot help. ~40–60 M reads per sub-plate would give essentially the same data as the ~100 M currently delivered. TapeStation, Qubit and library qPCR measure **mass**, not **complexity**, which is why pre-sequencing QC did not predict this.

---

## 2. Where the reads go

Per-sub-plate funnel, plates 21+22 (`ad_hoc_checks/plate_qc_summary/output_plate21_22/stage_summary.csv`):

| stage | retention |
|---|---|
| sequenced → demultiplexed | 94–97 % |
| → adapter/dimer filtered | 88–94 % |
| → aligned | 83–92 % |
| → **deduplicated (unique)** | **2.3 – 11.7 %** |

Every upstream step is healthy. **88–98 % of mapped reads are PCR duplicates.**

Per-well, plates 21+22 (`well_metrics.csv`, n = 768):

| | median demux reads | median unique | median dup rate |
|---|---|---|---|
| PASS wells (n = 73) | 6,971,298 | 520,727 | 90.3 % |
| FAIL wells (n = 686) | 49,834 | 3,268 | 88.9 % |

Two things follow. **(a)** The ~140× spread exists in *raw reads assigned per barcode*, i.e. before any bioinformatics. **(b)** The duplication rate is the same in failing and passing wells, so every well was sequenced to saturation and the read distribution is a faithful readout of the *molecular* distribution — it is not "some wells were under-sequenced".

The distribution of log10(usable reads) is **bimodal** with a valley around 10⁴·⁵–10⁵.

---

## 3. Hypothesis status — check this before proposing work

| # | hypothesis | status | evidence |
|---|---|---|---|
| **H1** | Adapter ligation is concentration-limited at 80 pM; yield ∝ p² | **SUPPORTED** | §4 |
| **H1b** | Adapter self-ligation is structural (self-complementary CATG + PNK), so dimers scale as [adapter]² | **SUPPORTED** | §4.3, §6.2 |
| **H5** | ~50 % of dispensed objects never yield amplifiable DNA (cell quality / DNA content) | **SUPPORTED — second bottleneck** | §7 |
| **H11** | The ~3k-read floor in failing wells is partly not real library | **CONFIRMED and quantified** | §5 |
| **H12** | No 3′ adapter trimming hides 40 % of the digest | **CONFIRMED — new, and fixable** | §6.3 |
| **H13** | Index hopping on a 3 × 3 combinatorial grid in one lane | **CONFIRMED — the larger contamination route** | §5.3 |
| **H6** | Evaporation / protease carryover / thermal gradient | **OPEN — per-well, still standing** | §7.1 |
| **H3** | 26–28 cycle four-block PCR jackpots the pool | **OPEN — plate-level, cannot be per-well** | §6.1 |
| **H10** | Over-aggressive UMI deduplication | **OPEN, unlikely** | the old benchmark was never run |
| **H2** | First 0.8× SPRI destroys / size-biases complexity | **REJECTED** | §6.2 |
| **H2b** | …and a pooled step could explain the between-well spread | **REJECTED — arithmetically** | §4.1 |
| **H4** | Incomplete NlaIII digestion | **REJECTED** | §6.4 |
| **H7** | Barcode-specific adapter failure | **REJECTED** | §3.1 |
| **H8** | Cell dispensing / sorting | **REJECTED** | §3.2 |
| **H9** | Sequencing depth | **REJECTED** | §1 |
| **H14** | Ambient DNA / dirty nozzle contaminating wells before library prep | **REJECTED** | §5.2 |

### 3.1 H7 — not the barcodes

Per-well success across the 8 sub-plates: observed **43 / 37 / 12 / 4** well positions succeeding in 0 / 1 / 2 / 3 of them, against a binomial expectation of **43.2 / 36.3 / 13.3 / 2.8**. Indistinguishable from independence — no systematically dead well, no defective oligo. The barcode set is also clean by design (96 × 8 nt, all 50 % GC, minimum Hamming distance 3).

### 3.2 H8 — not the dispensing

Four independent lines:

1. **The controlled experiment already ran.** FACS (plate17) = 18.0 % PASS; CellenONE (plates 21/22), which images every droplet and verifies one cell per well, = 10.7 % / 8.3 %. The technology that fixes dispensing made it worse.
2. **Nothing CellenONE measured predicts read yield.** 11 metrics, 0 significant against gate group, 0 against raw read count, every |ρ| < 0.075, every BH-corrected p = 1.00 (`docs/CELLENONE_VS_GATE.md`). The cleanest illustration: 8.2 % of PASS vs 8.5 % of FAIL wells had a second object in the ejection zone, Fisher p = 1.00.
3. **CellenONE isolated 384/384 wells** on both plates, cross-checked against the vendor PDF.
4. **The variance sits between sub-plates that share one cell suspension.** plate17_1 gave 17.1 M molecules and 34/96 PASS wells; plate17_2 gave 3.4 M and 8/96 — a 5.1× swing between four aliquots of the same sort (plate21 4.6×, plate22 3.6×). Cells cannot do that; library prep can.

**Keep the CellenONE.** It is not the problem, and the images are the only forensic tool for "was there even a cell in there". It is just not the lever.

---

## 4. The surviving explanation

### 4.1 The elimination is structural, and it is the strongest argument in the whole investigation

**A step applied to the pool cannot create between-well differences.** The 0.8× SPRI is performed on the *pooled* ligation; so is the 24–28 cycle PCR; so is sequencing. Every well's molecules pass through the same beads in the same tube. These steps set how much library the plate makes; they are arithmetically incapable of deciding which wells get it.

| step | shared or per-well | status |
|---|---|---|
| pooling, 0.8× SPRI, 24–28 cycle PCR, sequencing | shared | **cannot create between-well variance** |
| contamination (chimera + hopping) | shared, post-ligation | **consequence, not cause** (§5.2) |
| cell dispensing | per-well | rejected (§3.2) |
| barcode / adapter identity | per-well | rejected (§3.1) |
| NlaIII digestion | per-well | rejected (§6.4) |
| ligation *specificity* | per-well | rejected (§6.4) |
| ambient DNA / dirty nozzle | per-well, pre-ligation | rejected (§5.2) |
| **lysis / protease inactivation** | **per-well** | **STILL STANDING** |
| **adapter ligation *efficiency*** | **per-well** | **STILL STANDING** |

### 4.2 Yield scales as the square of ligation efficiency

A forked Y-adapter must ligate to **both** ends of an NlaIII fragment for that fragment to be exponentially amplifiable by the P5/P7 primer pair. (Confirmed by the read structure: R1 and R2 both begin with a well barcode plus linker.) So per-well yield ∝ p², where p is the per-end efficiency.

A1 measured the ceiling: **13,815,697 CATG sites** in GRCh38_canonical, 13,815,673 fragments, of which **7,132,080 (51.6 %)** fall in the 79–399 bp window actually recovered. Counting distinct fragments recovered as complete pairs (both ends on consecutive CATG sites), plate21:

| | distinct fragments | % of the 6.79 M window | implied p |
|---|---|---|---|
| best well | 432,905 | 6.4 % | ~18 % |
| median PASS well | 157,032 | 2.3 % | ~11 % |
| median FAIL well | 537 | 0.008 % | ~0.6 % |

plate22 is the same shape (best 451,455 → p ≈ 18 %; median PASS 82,204 → p ≈ 8 %; median FAIL 473 → p ≈ 0.6 %).

**The spread in p is 13–17×.** The pre-A1 prediction was that ~12× was what the square law needed to produce the observed ~160× spread in yield.

> **Honest caveat, do not overstate this.** Taking the square root of a yield ratio is a *definitional* transform, not an independent test. What A1 genuinely added is that the ceiling and the recovery are measured rather than assumed, and that the absolute efficiency is now pinned down. A plate-average p of ~6–11 % is a ligation running badly across the board, not one failing catastrophically in some wells and working in others.

### 4.3 Why p is low: the adapter is at 80 pM and self-ligates by design

**0.4 µL of 2 nM into ~10 µL = 0.08 nM = 80 pM.** That is four to five orders of magnitude below a standard Illumina ligation (1–15 µM), with a 4-nt CATG overhang whose annealing is transient at 16 °C. A sticky-end ligation is bimolecular — its rate depends on absolute concentration, not only on the molar ratio.

The titration history that led there (`scDNAseq pilot.pptx`):

| plate | adapter | gDNA : adapter | outcome |
|---|---|---|---|
| 3 & 4 | 6.7 nM | 1 : 447 | ">95 % self-ligation adapter dimer" |
| 5 | 0.67 nM | 1 : 45 | "more fragments of interest" |
| 6 | 0.133 nM | 1 : 9 | "DNA lost somehow in the process" |
| 7 → present | **0.08 nM** | 1 : 9 | 20/96 decent wells; now 8–18 % |

**And the reason there is no good operating point is structural.** CATG is self-complementary, so two adapters anneal to one another exactly as an adapter anneals to a genomic end — and because T4 PNK is in the annealing mix, phosphorylating both oligos, those self-ligations are covalently sealed. Therefore:

- dimer formation ∝ **[adapter]²**
- genuine ligation ∝ **[adapter] × [fragment ends]**

One cell supplies ~10⁷–10⁸ NlaIII ends against 4.8 × 10⁸ adapter molecules, so adapter is always in large excess and the squared term always wins. Lowering the adapter suppresses dimers faster than it suppresses library — which is what the titration observed — but it drives the genuine reaction into the concentration-limited regime.

**Predicted dimer size, computed from the oligo sequences:** 55 (TYA) + 64 (BYA) + 37 (Ad1 tail) + 32 (Ad2 tail) = **188 bp**. A genuine molecule is 188 bp + insert. The 0.8× SPRI cut is ~200 bp with a soft 150–300 bp transition, so the dimer sits inside the transition and no single-sided bead ratio can cleanly separate 188 bp from ~240 bp. **A2 confirmed the dimer was never removed** (§6.2).

### 4.4 Two things worth testing that follow directly

1. **PEG-6000/8000 at 5–10 % in the ligation.** Molecular crowding raises the effective local concentration by orders of magnitude without adding a single adapter molecule, so it should raise genuine ligation without raising dimers proportionally. No record of it having been tried.
2. **An unphosphorylated BYA oligo** (order without 5′-P, drop PNK). In an adapter–adapter junction *both* nicks need an adapter 5′-phosphate; in an adapter–genome junction one nick can be sealed using the genomic 5′-phosphate NlaIII leaves. **This is a hypothesis to test, not an established fact** — with an unphosphorylated bottom strand the adapter is held on by an 18 bp stem (Tm ≈ 52–56 °C), which may melt during the 72 °C 5-min fill-in that opens the PCR. If it does, lengthen the stem or lower the fill-in temperature.

Also noted: pilot deck slide 22 records *"Increase reaction volume to 4-6-8-10."* Raising the volume lowers both the DNA and the adapter concentration, which slows a bimolecular ligation further. It helps pipetting accuracy; it works against the reaction.

---

## 5. A4 — what the failing wells actually contain

Full method: `ad_hoc_checks/background_floor/README.md`. Measured on plates 21 and 22: 8 sub-plates, 768 wells, 389 M barcode-assigned read pairs, 24.5 M deduplicated molecules.

### 5.1 Composition and the corrected read counts

Of **686 FAIL wells**:

| | wells | share |
|---|---|---|
| **genuine low yield** — a real, tiny library | **367** | **53.5 %** |
| mostly contamination, detectable genuine remnant | 195 | 28.4 % |
| **empty** — upper 95 % bound on genuine molecules < 50 | **124** | **18.1 %** |

Aggregated, **27.6 %** of all FAIL-well molecules are demonstrably another well's.

| plate | status | median contamination | median usable reads | **corrected** |
|---|---|---|---|---|
| plate21 | PASS | 4.2 % | 614,056 | 584,122 |
| plate21 | FAIL | **48.5 %** | 3,083 | **925** |
| plate22 | PASS | 4.6 % | 395,918 | 377,603 |
| plate22 | FAIL | **36.2 %** | 3,477 | **1,979** |

No well changes gate class under correction. **But note the direction: correcting makes the failing wells look worse**, so the true PASS-to-FAIL spread widens to roughly **630×**.

`c_corrected` is a **lower bound** — it cannot see contamination whose source is outside this run (reagent DNA, carry-over from an earlier plate), because there is no second well to match against. Everything it counts is positively identified.

### 5.2 The UMI dates the contamination — which is what rejects H14

The UMI is carried on the adapter and attached **during that well's own ligation**. So free-floating tissue DNA that landed in a well before library prep — ambient DNA, a lysed neighbour, a dirty nozzle — would be digested and ligated *in that well* and would pick up *that well's* fresh, independent UMIs. It could produce shared **positions**, but never a shared **exact UMI at the same base**.

Two wells can hold the same (position, strand, UMI) only if the same physical adapter-ligated molecule reached both, which can only happen **during PCR or on the flowcell**.

For the median FAIL well (1,592 molecules):

| donor set | shared molecules | chance expectation |
|---|---|---|
| top-5 wells of its sub-plate | **54** | 0.002 (analytic) |
| every other well of its sub-plate | **92** | 0.11 (permutation) |
| any other well in the run | **262** | 0.50 (permutation) |

Three to five orders of magnitude above chance. **The exact-UMI test does not merely detect the contamination; it dates it, and excludes everything before ligation.**

Two independent corroborations:

- For molecules held by exactly two wells, the BC2 tag the demux discards **names the partner well 48.3–48.6 % of the time against 1.05 % by chance** — adapter-mediated chimera formation observed molecule by molecule.
- The barcode-derived and coordinate-derived well-pair matrices, built from completely different inputs, agree on *which* wells exchange material at Spearman **r = 0.79–0.88**.

### 5.3 Two routes, separated — and the larger one is a submission choice

The 3 × 3 combinatorial index grid in one lane (see `01_METHODS.md` §B6) separates the mechanisms cleanly:

- **Within a sub-plate all 96 wells share one i7+i5 pair**, so index hopping cannot move a read between them. Any within-sub-plate transfer must be inline-barcode recombination during PCR.
- **Across libraries**, a single hop preserves the inline barcode, so it must land on the **same well id**. Measured: libraries differing in one index share **0.807 %** of molecules at the same well id vs **0.0008 %** at a different one (**975×**); pairs differing in two indices share **50× less**.

| route | molecules | % of FAIL molecules | mechanism |
|---|---|---|---|
| within own sub-plate | 188,410 | **10.97 %** | PCR chimera |
| same well id, other sub-plate | 158,994 | 9.25 % | index hopping |
| other plate | 139,862 | 8.14 % | index hopping |
| other well id, other sub-plate | 15,825 | 0.92 % | mixed |
| **union** | **474,576** | **27.62 %** | |

**Index hopping is the larger route (~17 pp), PCR chimera the smaller (~11 pp)** — which inverts the priority of the two available fixes.

### 5.4 The barcode swaps are molecular, not sequencing error

Of 389 M BC1-assigned reads, **701,358 carry a BC2 that is a different well's valid barcode** — **25.7× above** a sequencing-error null measured from the observed distance spectrum rather than assumed.

| Hamming d | observed | expected if error | obs/exp |
|---|---|---|---|
| 1 | **0** | – | – |
| 2 | **0** | – | – |
| 3 | 94,728 | 10,500 | 9.0× |
| 5 | 191,682 | 6,471 | 30× |
| 7 | 86,395 | 901 | 96× |
| all | 701,358 | 27,341 | **25.7×** |

`d = 1` and `d = 2` being exactly zero is the internal check that barcode parsing is correct (the code's minimum distance is 3). **obs/exp rising with distance is impossible for an error process** — its mass must pile up at the smallest reachable distance. The matrix is symmetric (slope of log(M[i,j]/M[j,i]) on log(a_i/a_j) = −0.06 to +0.03; error predicts 1) and its rate scales with the partner's depth (slope 0.50–0.75; error predicts 0). That is the product signature `κ·a_i·a_j` of a swap, not the donor signature `a_i·e_ij` of an error.

### 5.5 The fixes, in order

1. **Unique dual indices.** Removes the ~17 pp hopping route outright. A `BC1 == BC2` filter **cannot** touch it — a hopped read carries both inline barcodes intact and correct.
2. **Enforce `BC1 == BC2` in the demux.** One line at `extract_umi_barcode.py:181`; BC2 is already parsed. Removes **33 %** of a FAIL well's contaminating molecules and **2.3 %** of its genuine ones (43× directional; PASS wells 28 % vs 0.59 %). Catches the chimera route only — 49.2 % of chimera-route molecules run-wide, 2.2 % of hop-route.
3. **Add a background-aware QC axis** — `c_corrected` from A4, or the cheaper percent-at-CATG from A1.
4. **Put cell-free control wells on future plates.**

---

## 6. A1 and A2 — what the chemistry and the cleanup actually do

### 6.1 The pooled steps cost you, but not selectively

**26–28 PCR cycles in four blocks** with a cleanup between each, from a pooled sub-picogram template. Close to a worst-case design for jackpotting, and mechanically why duplication is 90–97 %. Plates 9&10 recorded **more cycles giving less product** (18 cyc → 30 ng/µL; 26 cyc → 8 ng/µL) — a classic sign of running past plateau. **Still untested (A7).** Note it cannot explain between-well spread (§4.1).

### 6.2 A2 — the bead cleanup is not eating the library, but it never removed the dimers

Renormalised inside the measurable 79–499 bp window, observed/expected by insert size:

| insert bin | plate21 | plate22 | plate17 |
|---|---|---|---|
| 79–99 | 1.53 | 1.51 | 0.58 |
| 100–149 | **2.15** | **2.03** | **1.37** |
| 150–199 | 1.26 | 1.36 | 1.51 |
| 200–299 | 0.38 | 0.43 | 1.00 |
| 300–499 | 0.15 | 0.15 | 0.51 |

**Short fragments are enriched, not depleted; the depletion runs against long fragments.** A 0.8× size selection would produce the opposite gradient. H2 rejected.

But the 188 bp dimer predicted in §4.3 sat on the cut-off and survived:

| plate | PASS median dimer % | FAIL median dimer % | FAIL wells > 50 % dimer |
|---|---|---|---|
| plate21 | 1.9 | **25.6** | **111 of 339** |
| plate22 | 2.6 | **12.9** | **103 of 347** |
| plate17 | 0.3 | 2.1 | 7 of 266 |

The dimers are not *causing* the low yield — they are what fills the sequencing when there is little real library to compete with them — but they consume a real share of reads. Note the plate ordering matches the CATG-concordance result exactly: plate17's failing wells are clean on both measures, plate21/22's are dirty on both.

### 6.3 A1/A2 — 40 % of the digest is structurally invisible, and it is recoverable

There is **no 3′ adapter trimming anywhere in the pipeline** and bowtie2 runs end-to-end. R1 is 79 bp after the fixed 22 nt trim, so a fragment shorter than the read runs into the adapter and cannot align.

| fragment class | fragments | % of digest |
|---|---|---|
| < 79 bp — unalignable end-to-end | 4,595,696 | 33.3 % |
| < 100 bp — 0 % recovery below 50 bp, 9–23 % for 50–99 bp | 5,585,477 | **40.4 %** |
| 79–399 bp — the window actually recovered | 7,132,080 | 51.6 % |

Recovering the sub-100 bp classes at the efficiency of the best-observed class (100–149 bp), for the median PASS well:

| plate | recovered now | extra | fold |
|---|---|---|---|
| plate21 | 157,032 | +221,045 | **×2.39** |
| plate22 | 82,204 | +110,352 | **×2.35** |
| plate17 | 49,519 | +48,547 | **×1.92** |

**Two caveats stated by the analysis, not to be dropped:** this assumes the rescued short class would be recovered as efficiently as the 100–149 bp class (an optimistic ceiling — very short fragments may ligate, amplify or cluster differently, and after trimming carry less unique sequence to map), and it holds per-fragment efficiency fixed, which implicitly assumes proportionally more sequencing. **It is a gain in recoverable complexity, not a free gain at fixed depth.** It also changes every existing BAM — a re-run decision, not a patch. **A11 is designed to measure the realised gain.**

### 6.4 A1 — the enzyme chemistry is clean

| evidence | plate21 | plate22 | plate17 |
|---|---|---|---|
| read ends exactly at a genuine CATG (raw, PASS wells) | **96.4 %** | **96.3 %** | **97.2 %** |
| pairs spanning exactly one NlaIII fragment (median PASS well) | 95.8 % | 95.7 % | 96.2 % |
| inserts skipping ≥1 uncut CATG | 2.0 % | 2.2 % | 3.2 % |
| soft-clipped alignments, of 1.5 bn read ends | **0** | **0** | **0** |

Digestion is complete to within a few percent; ligation is specific. **H4 rejected.**

The residual off-target ends are not random shearing either: **79–94 % of read ends that miss the site map sit on 4-mers one substitution from CATG**, against 4.7 % expected by chance — a 15–20× enrichment, which is what star activity or partial-site cleavage looks like. Off-map ends are 4.7 % of PASS reads and 15.6 % of FAIL reads.

*Caveat:* the "skips an uncut site" measure is conditioned on what was sequenced, and uncut molecules are longer and recovered less efficiently, so 2–3 % is a **lower bound** on the per-site uncut rate in the tube.

### 6.5 A1 — "FAIL" means different things on different plates

Median % of read ends at a genuine CATG, PASS vs FAIL, dedup BAMs:

| plate | median PASS | median FAIL | p | FAIL wells that still look like NlaIII library |
|---|---|---|---|---|
| plate17 | 92.1 % | 91.4 % | 0.37, n.s. | **242 of 265 (91 %)** |
| plate21 | 93.9 % | 83.8 % | 5 × 10⁻¹³ | 104 of 333 (31 %) |
| plate22 | 92.7 % | 69.9 % | 9 × 10⁻¹³ | 74 of 318 (23 %) |

On plate17 essentially every failing well is a real cell that was simply too shallow. On plates 21 and 22 most contain reads that are largely not NlaIII product. **Two different failure modes have been counted as one, and read count alone cannot separate them.** Percent-at-CATG costs one pass over the dedup BAM. A4 reproduced the same plate ordering by a completely independent measure.

### 6.6 A methodological note that has bitten three times

**Compute your null; do not assume it.** Three separate results in this project turned on it:

- A4's position-sharing test: the naive null (uniform over the digest) gives 7.0 %, but the correct held-out-real-library null is **22.3 %** — so an eyeballed "50 % of positions shared!" is mostly ordinary coincidence, and the position test is nearly uninformative here.
- A4's UMI collision probability: measured `q̂ = 2.715 × 10⁻⁷`, **4.55× the uniform 4⁻¹²**, because UMI position 1 is 60 % A. Using 4⁻¹² would have understated chance collisions fourfold.
- A4's barcode-error null: measured from the observed distance spectrum rather than assumed, so every correlated instrument artefact is already inside it.

Also: **poly-G UMIs must be removed.** On a two-colour NovaSeq X, G is "no signal", so a dark cluster spells poly-G. `GGGGGGGGGGGG` alone occurs 12,681 times in plate21_1/W86 and appears in every well. Dropping UMIs with ≥10 of 12 identical bases removes 2.94 % of read 1s.

---

## 7. Cell quality — the second bottleneck

### 7.1 A plate-position gradient exists, in both sequencing and qPCR

- Plates 21/22 sequencing: Kruskal-Wallis across rows p = 0.0039, across columns **p = 4.2 × 10⁻⁷**; well number vs usable reads Spearman ρ = −0.14 (p = 8 × 10⁻⁵). Best-to-worst column ≈ 7× in mean log10 yield.
- plate23 qPCR: rows A–H **20.8 %** positive vs rows I–P **36.5 %**.

Crucially, in plate23 the row effect is **independent of cell size** (row index vs diameter ρ = −0.013, p = 0.79; in a joint logistic model both terms stay significant). So it is a **handling / thermal / evaporation gradient, not a cell gradient**. Candidates: evaporation across three long incubations under only 5 µL of oil; multichannel dispensing order; incomplete protease inactivation varying with block position; thermocycler non-uniformity. Effect size (~7×) is real but an order of magnitude smaller than the well-to-well lottery. **Open — see A8.**

### 7.2 Cell size predicts library formation (plate23)

`DATA/384_well/plate23 - qPCR QC/plate23_well_classification.csv`, n = 384: **110 positive / 75 ambiguous / 199 negative**.

| metric | AUC (pos) | p | AUC (pos+amb) | p |
|---|---|---|---|---|
| **diameter** | **0.716** | 3.8 × 10⁻¹¹ | **0.752** | 1.4 × 10⁻¹⁷ |
| elongation | 0.481 | 0.56 | 0.462 | 0.20 |
| circularity | 0.481 | 0.57 | 0.459 | 0.16 |
| transmission intensity | 0.445 | 0.09 | 0.427 | 0.013 |
| Orange (**CD13**) | 0.431 | 0.014 | 0.443 | 0.026 |
| Blue (**CD24**) | 0.463 | 0.22 | 0.480 | 0.46 |
| Red (**VCAM1**) | 0.500 | 1.00 | 0.500 | 1.00 |

Monotone across the whole range, not a debris threshold:

| diameter (µm) | n | % positive | % pos+amb |
|---|---|---|---|
| ≤ 11 | 25 | 8 % | 20 % |
| 11–13 | 62 | 8 % | 18 % |
| 13–15 | 81 | 19 % | 32 % |
| 15–17 | 107 | 34 % | 59 % |
| 17–19 | 57 | 40 % | 70 % |
| > 19 | 52 | **56 %** | 77 % |

**Interpretation:** a continuous, monotone diameter effect is what you expect when yield is proportional to input DNA rather than saturating — the signature of a concentration-limited ligation. Independent corroboration of §4, from a different assay.

**Red/VCAM1 has AUC exactly 0.500 with p = 1.00 because the channel is empty** (signal in 27 of 384 wells; the pipeline drops it run-wide). That row is a null result about the *channel*, not about VCAM1 biology. Do not report it as evidence that VCAM1 status does not matter.

**Why plates 21/22 hid this:** `CELLENONE_VS_GATE.md` found diameter unrelated to reads (p = 0.71). Not a contradiction — when ~90 % of wells fail for a downstream reason, the downstream failure dominates and a modest upstream signal cannot show through. plate23 measures a step *before* that noise is added.

### 7.3 ⚠ Selection-bias risk for the biological result

If larger cells preferentially succeed and size tracks DNA content, the wells surviving to copy-number calling are **enriched for S-phase and polyploid cells** — exactly the cells that generate spurious CN profiles. Slide 1 of the pilot deck already flagged triploid/tetraploid cells needing correction. **For a project about copy number in *normal* kidney this threatens the conclusion, not just the yield.** Untested — see A5.

### 7.4 The two bottlenecks, decomposed

```
P(sequencing success) = P(makes a library at all) × P(reaches 100k | makes a library)
       8–18 %         =        ~29–48 %           ×          20–35 %
```

The first term is plate23's qPCR positivity; the second is the gap that pooling, cleanup and PCR impose. **A10 measures both directly once plate23 is sequenced** — that is the single highest-value experiment available.

---

## 8. Open questions

Ranked by how much they would change the picture.

1. **What is plate23's qPCR actually measuring?** Primer sequences, target, and the exact protocol step at which the plate was sampled. Whole well or an aliquot? At ~2 genome copies per well, sampling 1 µL of 6 µL means most true positives read negative by Poisson chance alone, and the true positivity rate would be far above 28.6 %.
   - **A hard constraint:** the protocol appendix's NEBNext library quant runs **20 cycles**, but plate23 reports Cq values up to 38. **Plate23 therefore did not use that assay as written.**
   - **A diagnostic to run first:** the positive wells have Tm **91.60 °C ± 0.22** — extraordinarily tight. A real multiplexed library of thousands of different inserts should melt *broadly*; a single sharp peak is what one defined species gives, and the obvious candidate is the 188 bp adapter dimer. Check the derivative-curve shape (`plate23/Comparative Ct with Melt_01_Melt Curve Raw_*.csv`) and run a few positives on a gel before treating "pos" as "made library".
2. **Was a per-well index PCR performed on plate23?** Decides whether the assay is adapter-facing or P5/P7.
3. **Adapter plate age and dilution history for plate17 vs plates 21/22.** Four serial dilutions, frozen and thawed. A single 2 nM plate shared across sub-plates prepped on different days would explain the between-sub-plate variance directly.
4. **Which sub-plates were prepped on which days, by whom, in which thermocycler?** After §4.1 this is the strongest unassigned lead in the project. It is not currently recorded anywhere.
5. **The plate21/22 TapeStation trace** (`Plate21&22/P21&22_post_cleanup.pdf`) — does a 188 bp dimer peak dominate after cleanup?
6. **Why is T4 PNK in the adapter annealing mix?** If the intent was simply "ligation needs a phosphate", the genomic 5′-phosphate NlaIII leaves may be sufficient on its own (§4.4).
7. **plate19's status** — only `K1184` appears processed under `DATA/384_well/plate19`. Was it ever sequenced as a 384 plate?

---

## 9. What to do

### Free, and now

- **Unique dual indices** on every future submission; do not put eight libraries on a 3 × 3 combinatorial grid in one lane. Removes the single largest contamination route.
- **Enforce `BC1 == BC2`** in the demux — one line, 43× directional.
- **Re-run one sub-plate with 3′ adapter trimming** and measure the realised gain against the predicted ×2.1–2.4 before committing to re-aligning 1152 BAMs (A11).
- **Cut sequencing depth** to ~40–60 M reads per sub-plate. You are paying ~7× for saturated libraries.
- **Add a background-aware second QC axis** — percent-at-CATG or `c_corrected`.
- **Add control wells** to every future plate: ≥4 no-cell and ≥4 with 6 pg bulk gDNA, spread across positions, not all in column 24.
- **Record per-sub-plate prep metadata** — date, operator, thermocycler, adapter plate lot.
- **Make `-X 500` a conscious choice.**

### Bench

- **W1 — sequence plate23.** The only dataset linking a per-well pre-library measurement to a per-well read yield. Change nothing else about how it is processed or the link is lost.
- **W2 — ligation matrix.** One plate, four or five arms, all with control wells: (i) current; (ii) + 5 % PEG-8000; (iii) volume back to ~4 µL; (iv) PEG + reduced volume; (v) unphosphorylated BYA. Read out by per-well qPCR **before pooling** — an answer in days, no sequencing.
- **W3 — cleanup test.** Split one pooled ligation: 0.8× vs 1.2× SPRI vs column, each ± carrier.
- **W4 — cycle reduction.** Replace the four-block 24–28 cycle scheme with qPCR-monitored amplification stopped in exponential phase, single cleanup.

### Analysis

See `03_ANALYSIS_PROMPTS.md`. Highest value outstanding: **A6** (absolute dimer counts per well — the discriminating test for whether dimer formation is a fixed background), **A11** (realised gain from 3′ trimming), **A12** (`BC1 == BC2` re-demux).

---

## 10. Confidence and provenance

| claim class | confidence | basis |
|---|---|---|
| The funnel, the complexity ceiling, the per-well spread | **High** | direct measurement, three plates, 1152 wells |
| H4, H2, H7, H8, H9, H14 rejected | **High** | each measured directly, effect sizes far from the boundary |
| H11, H12, H13 confirmed | **High** | A1/A2/A4, with internal cross-checks between independent measures |
| p ≈ 6–11 % and the 13–17× spread | **Medium-high** | measured ceiling and recovery; the square root is definitional, and p is an end-to-end efficiency, not purely ligation |
| The 188 bp dimer size | **Medium-high** | computed from oligo sequences; corroborated by A2's dimer rates; not yet confirmed on a TapeStation trace |
| The unphosphorylated-BYA proposal | **Hypothesis only** | reasoning from sequence and protocol; the 18 bp stem may not survive the 72 °C fill-in |
| Everything about plate23's assay | **Low** | the assay is not identified; see §8.1 |
| Lysis vs ligation as the surviving per-well cause | **Not yet separated** | both survive the elimination; W2 is designed to separate them |
