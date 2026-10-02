# Methods — wet lab and dry lab

Read `00_START_HERE.md` first for the path map and conventions. This document is the reference for **how the material and the data were made**. Findings live in `02_FINDINGS.md`.

Source of truth for the bench protocol: `scDNAseq Template.docx` (and its dated per-plate copies, e.g. `scDNAseq - 20260211 - FF.docx`). Source of truth for the pipeline: `PIPELINE/config.yaml` and `PIPELINE/Snakefile`. Where this document states a number, it has been checked against those files; where something is uncertain it says so.

---

# PART A — WET LAB

## A1. What the chemistry is

A home-made single-cell WGS method in the scNLAIII-seq / scKaryo-seq lineage. The logic:

1. One cell per well, lysed in place under oil.
2. **NlaIII** digests the genome at every `CATG`, leaving 4-nt **3′ CATG overhangs**.
3. A **forked (Y) adapter carrying a well-specific 8-nt barcode and a 12-nt UMI** is ligated to those ends. This is where well identity enters the molecule.
4. All 96 wells are **pooled**, cleaned up, and amplified together with Nextera P5/P7 primers.
5. Sequenced paired-end; wells are recovered computationally from the inline barcode.

Because barcoding happens in-well and amplification happens in-pool, **every step from pooling onwards is shared by all 96 wells**. That fact carries a lot of weight in the failure analysis — see `02_FINDINGS.md` §4.

## A2. The Y-adapter — structure and sequences

Source: `Y adapter oligos.xlsx`, sheets `Oligos`, `Well barcode`, `Nextera adapters`.

Two oligos per well, ordered as 96 pairs, 100 nmol scale, standard desalting.

```
TYA_W##  (top,    55 nt)  5'-[Nextera R1 + ME  33 nt]-[WB 8 nt]-[linker ACGCTGACTGCATG 14 nt]-3'
BYA_W##  (bottom, 64 nt)  5'-[CAGTCAGCGT 10 nt]-[WB' 8 nt]-[N12 UMI]-[ME' + Nextera R2  34 nt]-3'
```

Worked example, W01 (`WB = TCTCATCG`):

```
TYA_W01  TCGTCGGCAGCGTCAGATGTGTATAAGAGACAG TCTCATCG ACGCTGACTGCATG
BYA_W01  CAGTCAGCGT CGATGAGA NNNNNNNNNNNN CTGTCTCTTATACACATCTCCGAGCCCACGAGAC
```

`reverse_complement(BYA[0:18])` equals `TYA[-22:-4]`, so the annealed adapter is:

- an **18 bp duplex stem**
- a **4-nt 3′ CATG overhang** on the top strand — the ligating end
- two single-stranded fork arms: Nextera Read 1 / P5 side on the top strand, Nextera Read 2 / P7 side on the bottom strand

**The barcode is on both strands** (`WB` on top, its reverse complement on the bottom), so it appears in both R1 and R2 after sequencing. **The UMI is on the bottom oligo only**, so it appears in R2 only.

### Barcode set

96 barcodes, all **8 nt**, all **exactly 50 % GC**, **minimum pairwise Hamming distance 3**. Distance spectrum: `d3:411 d4:836 d5:1335 d6:1334 d7:593 d8:51`. No single base-call error can convert one well's barcode into another's — a property that A4 used as an internal control (`02_FINDINGS.md` §5).

### Nextera index primers

```
Ad1 (i5)  AATGATACGGCGACCACCGAGATCTACAC [i5 8nt] TCGTCGGCAGCGTCAGATGTGTAT     61 nt
Ad2 (i7)  CAAGCAGAAGACGGCATACGAGAT [i7 8nt] GTCTCGTGGGCTCGGAGATGTG           54 nt
```

Eight Ad1 (N501–N508) and twelve Ad2 (N701–N712) are on file. **Which pair each library used is recorded per plate in the protocol document, step 33.** This matters — see the index-hopping finding.

### Adapter annealing and dilution

| component | volume | note |
|---|---|---|
| TYA_W## (100 µM) | 1 µL | |
| BYA_W## (100 µM) | 1 µL | |
| 10× T4 Ligation Buffer (with ATP) | 1 µL | |
| NF water | 6.5 µL | |
| **T4 PNK** | **0.5 µL** | **phosphorylates the 5′ end of both oligos** |
| total | 10 µL | annealed product ≈ 10 µM |

Program: 37 °C 30 min → 95 °C 5 min → ramp to 4 °C at ~2 °C/min.

Then serially diluted **10 µM → 1 µM → 100 nM → 10 nM**, stored at −20 °C, and on the day diluted to the **2 nM working plate** (10 µL of 10 nM into 40 µL ligation buffer).

> **Two things about this to keep in mind.** (i) PNK phosphorylates *both* oligos, and the CATG overhang is self-complementary, so the adapter can form covalently sealed self-ligation products — the protocol's own note at step 20 ("Add the Y-adapter before the ligase to prevent adapter self-ligation") shows this was known and mitigated by order of addition rather than by design. (ii) The working concentration is reached by **four serial dilutions** of a hand-made plate, frozen and thawed between uses.

## A3. Bench protocol, step by step

Reagents of note: QIAGEN Protease 19155, NlaIII NEB R0125S (10 U/µL), T4 DNA Ligase NEB M0202S (400 U/µL), NEBNext Ultra II Q5 Master Mix, SPRIselect, silicone oil (Merck 317667), DNA LoBind plasticware throughout. All work in a PCR hood.

### Day 0 — sorting

| step | action | volume in well |
|---|---|---|
| | plate pre-loaded with 3 µL lysis buffer + 5 µL silicone oil | |
| 1–3 | sort one cell per well, spin 1000 × g 3 min, store −20 °C | **4 µL aqueous** under 5 µL oil |

Lysis buffer = 450 µL reconstituted Protease + 50 µL 10× NEBuffer 4 ⇒ **~1 mAU/µL protease in 1× NEBuffer 4**.

### Day 1 — lysis

| step | action | volume |
|---|---|---|
| 7 | add 2 µL lysis buffer per well | **6 µL** |
| 9 | 55 °C **overnight** → 75 °C 20 min → 80 °C 5 min | |

> The protocol flags the 75/80 °C inactivation as **essential**, because residual protease digests the restriction enzyme and the ligase downstream.

### Day 2 — digestion and ligation

| step | action | volume | working amount |
|---|---|---|---|
| 10–13 | add 2 µL fragmentation mix (per 250 µL: 25 µL 10× rCutSmart, 25 µL NlaIII 10 U/µL, 25 µL 3 mg/mL BSA, 175 µL water) | **8 µL** | **2 U NlaIII per well** |
| 15 | 37 °C **3 h** → 65 °C 20 min | | |
| 20 | add **0.4 µL of 2 nM Y-adapter** | 8.4 µL | **0.08 nM = 80 pM final** |
| 22 | add 1.6 µL ligation mix (per 200 µL: 20 µL T4 ligase 400 U/µL, 20 µL 10× ligase buffer, 80 µL 10 mM ATP, 80 µL water) | **10 µL** | **64 U T4 ligase per well** |
| 24 | 16 °C **overnight** → 65 °C 10 min | | |

> The reaction buffer at digestion is mostly **NEBuffer 4** carried over from the lysis step, with only ~0.25× rCutSmart contributed by the fragmentation mix. NlaIII is specified for rCutSmart.

### Day 3 — pooling, cleanup, amplification

| step | action |
|---|---|
| 25 | pool **all wells** (aqueous + oil) into one LoBind tube, discard the oil phase, transfer the aqueous phase, record the volume |
| 26–32 | **0.8× SPRI #1** — make up to 500 µL, add 400 µL beads, 85 % EtOH wash, elute in **21 µL** |
| 33–34 | **PCR block 1**: 20 µL template + 25 µL Q5 2× ReadyMix + 2.5 µL Ad1 + 2.5 µL Ad2. Program: **72 °C 5 min**, 98 °C 30 s, then **3 cycles** of 98 °C 10 s / 65 °C 30 s / 65 °C 45 s |
| 35–39 | **0.8× SPRI #2**, elute 21 µL, quantify |
| 41–42 | **PCR block 2** — same primer pair, **5 cycles** |
| 43–44 | **0.8× SPRI #3**, quantify |
| 45 | **PCR block 3** — **10 cycles**. *Total so far: 18 cycles.* **0.8× SPRI #4** |
| 46–47 | **PCR block 4** — additional cycles calculated for the sequencing scale (6–10 observed), aiming for 20 µL at 10 ng/µL. **Total 24–28 cycles.** |
| 48–49 | **0.8× SPRI #5**, quantify |
| 50 | Nanopore shallow long-read sequencing for QC |

**Five 0.8× SPRI cleanups in total.** The first is the important one: it runs on the pooled ligation *before any amplification*, at which point the molecules **are** the library complexity.

### Library quantification qPCR (protocol appendix)

NEBNext Illumina Library Quant: 10 µL standard or sample + 40 µL master mix. Program **95 °C 1 min, 20 × (95 °C 15 s / 63 °C 1 min), melt.**

> Note the cycle count. Plate23's per-well qPCR reports Cq values up to 38, which cannot come from a 20-cycle program — so **plate23 did not use this assay as written**. What it did use is still unresolved (`02_FINDINGS.md` §8).

## A4. Cell dispensing

Two methods have been used on 384-well plates:

- **FACS** (plates 13–20). Cells sorted directly into the pre-loaded plate.
- **CellenONE** (plates 21–23). A piezo-acoustic single-cell dispenser that images every droplet before ejection and records per-cell measurements.

### CellenONE channel assignment — supplied by the experimenter, not derivable from the run folder

The `.Run` folder records only the hardware channel name and its threshold. It contains no fluorophore, filter or wavelength information. **The mapping is:**

| CellenONE channel | marker | status |
|---|---|---|
| **Transmission** | — | the only channel set to `Positive`, i.e. the only one that gated cell selection |
| **Blue** | **CD24** | recorded, `Inspect` only |
| **Orange** | **CD13** | recorded, `Inspect` only |
| **Red** | **VCAM1** | **did not work well** — signal present in only 27 of 384 wells on plates 21/22, and the pipeline drops the channel from rendering because it is empty run-wide |

> Because all three fluorescence channels were set to `Inspect` rather than `Positive`, **fluorescence never influenced which cells the instrument picked**. Selection was on transmission (size/shape) alone. Any analysis that treats a CellenONE well as marker-sorted is wrong: the marker gate was applied upstream at FACS, and the CellenONE fluorescence is an observational record only — and an unreliable one for VCAM1.

### CellenONE geometry and gates

| quantity | value | source |
|---|---|---|
| ejection line ("green") | 396 px | the `.par` file's `Ejection (pix)` |
| purple line | 647 px | width of `BackgroundEjZone` |
| `EjBound` in the `.xls` | 330 | **NOT the ejection line — do not use it** |
| detection gate `DetDiaMinTrans` | 5 µm | plates 21 and 22 |
| isolation gate, plate21 | 12–25 µm, constant | |
| isolation gate, plate22 | **varied mid-run**: 9–25 µm (20 wells), 9–30 µm (72 wells), 12–30 µm (292 wells) | judge each well against its own row's criteria |
| pixel scale | plate21 1.407 px/µm, plate22 1.392 px/µm | fitted per run against the table's `Diameter` |

The nozzle is at the **left**; cells sediment leftwards toward it, so an object far to the right is still up the capillary and will not be ejected.

Full provenance of every CellenONE-derived metric is in `PIPELINE/docs/CELLENONE_AND_QC.md`. The short version: `cellenone/wells_raw.tsv` holds the instrument's own measurements plus that well's own criteria; `cellenone/objects.tsv` holds *our* pixel-detected objects; the "droplet call" is **display-only and must stay that way** — it never enters the QC gate.

---

# PART B — DRY LAB

## B1. Read structure — the single most useful thing to know

Straight off a raw FASTQ (`plate21_1_R1.fastq.gz`):

```
R1   ACCGTGAA ACGCTGACTGCATG AGTGAACTCCCATTCACA...
     |8 nt bc| 14 nt linker | genomic                8 + 14 = 22 = r1_trim_5prime

R2   ATGTTATAGTTA ACCGTGAA ACGCTGACTGCATG GAATGTTCTTC...
     | 12 nt UMI | 8 nt bc | 14 nt linker | genomic  12 + 8 + 14 = 34 = r2_trim_5prime
```

**The last 4 nt of the 14-nt linker are the NlaIII site itself.** The enzyme leaves a 3′ CATG overhang, the adapter carries the complementary overhang, and after ligation that CATG sits once in the molecule. The fixed 22/34 trim removes it, so **a perfect read starts at the base immediately 3′ of a genomic CATG**.

Consequences you can rely on (all verified in A1):

```
forward read (5′ end = leftmost base)   cut position = reference_start - 4
reverse read (5′ end = rightmost base)  cut position = reference_end
```

CATG is palindromic, so it reads `CATG` on the + strand either way. Reads are 101 bp; after the 22 nt R1 trim, **79 bp are aligned**. bowtie2 runs end-to-end, so **soft clipping is exactly zero** across 1.5 billion read ends — `reference_start` / `reference_end` are the true read ends and need no clip correction.

The linker also contains the reverse complement of `CAGTCAGCGT`, which is what the dimer filter keys on.

## B2. Pipeline stages

```
demultiplex → filter_dimers → align → deduplicate → (FastQC/MultiQC)
    → run_aneufinder → render_well_profiles → qc_review (review.html)
    → [human saves qc_decisions.csv] → validate_qc_decisions → derive_included_wells
    → run_aneufinder_reviewed → render_reviewed_profiles → cn_review (final)
```

Modes, selected by `MODE=…` in `submit_pipeline.sh`:

| mode | target |
|---|---|
| `preprocessing` | demux, dimer filter, align, dedup |
| `qc` | the above + FastQC + MultiQC |
| `aneufinder` | blacklist / GC precheck / AneuFinder only |
| `pre_review` (alias `full`) | everything up to `qc_review/review.html` |
| `review` | render the review report |
| `post_review` | the whole post-review chain in one run |
| `validate_review`, `aneufinder_reviewed`, `cn_review` | lower-level, for partial re-runs |

## B3. Key parameters (from `config.yaml` — verify, do not assume)

| parameter | value | note |
|---|---|---|
| `umi_length` | 12 | from the start of R2 |
| `barcode1_length` / `barcode1_offset` | 8 / 0 | in R1 |
| `barcode2_length` / `barcode2_offset` | 8 / 12 | in R2, after the UMI |
| `r1_trim_5prime` / `r2_trim_5prime` | 22 / 34 | |
| `adapter_sequence` | `CAGTCAGCGT` | dimer filter key |
| `filter_both_reads` | true | a pair is a dimer if the sequence is in both reads |
| `umi_tools_method` | `directional` | error-aware UMI grouping |
| `bowtie2_params` | `--very-sensitive` | end-to-end; **no 3′ adapter trimming anywhere**; `-X 500` left at default |
| `aneufinder.binsize` / `method` | 1,000,000 / `edivisive` | |
| `aneufinder.min_mappability` | 0.85 | |
| `aneufinder.min_reads_for_model` | 100 | wells below this are not staged — AneuFinder aborts the whole batch on a near-empty well |
| `qc_review.usable_reads_pass_cutoff` | 100,000 | |
| `qc_review.usable_reads_warn_cutoff` | 50,000 | |
| `plate_format` | 96 or 384 | |

`usable_reads` = `SN reads mapped:` from `bam/{well}.stats.txt`, i.e. **mapped reads in the deduplicated BAM, counting both mates**. A "molecule" in the A4 analysis is one read-1 alignment, so `molecules ≈ usable_reads / 2`. Watch this factor of two when comparing documents.

Gate → status → default decision: `≥100k → PASS → PASS`; `50k–100k → WARN → REVIEW`; `<50k → FAIL → EXCLUDE`; missing → `UNKNOWN → REVIEW`. **Duplication rate is never gated** — it is shown as a labelled diagnostic only.

## B4. 384-well mode

A 384-well plate is four interleaved 96-well sub-plates, each sequenced as its own FASTQ pair. Per-sub-plate preprocessing is unchanged; one plate-level AneuFinder + QC/CN review runs over all 384 wells.

The whole 96/384 branch reduces to two strings defined at the top of the `Snakefile`:

```python
SUB_BASE   # str(PLATE_DIR / "{sub}")  @384   |   str(PLATE_DIR)  @96  <- NO wildcard
LOG_TAG    # "{sub}" @384              |   PLATE @96
def wid(sub, well): return f"{sub}_{well}" if IS_384 else well
```

**In 96 mode `SUB_BASE` has no wildcard**, so every 96-mode path renders character-identical to what it was before 384 mode existed — same on-disk outputs, same `.snakemake/metadata` keys, no spurious re-runs. That property is the safety mechanism; preserve it. `wildcard_constraints` for `sub` and `well` are mandatory (Snakemake wildcards default to `.+`, which matches `/`).

`wid()` is the single definition of well identity: staged BAM basename → AneuFinder model id → `.RData` basename → plot basename → HTML well key → `included_wells.tsv` row → `sample_id`. That chain is why `run_aneufinder.R` needs zero edits for 384 mode.

**Regression gate before any `Snakefile` change:** a `MODE=pre_review` dry run against an already-complete 96-well plate must still say "Nothing to be done"; then a 384 dry run must contain no `align` / `deduplicate_bam` / `fastqc` / `multiqc` / `demultiplex` jobs.

## B5. Known defects in the current pipeline

These are established, not suspected. Each is quantified in `02_FINDINGS.md`.

1. **The demux assigns on BC1 alone.** `workflow/scripts/preprocessing/extract_umi_barcode.py:181` reads
   ```python
   if bc1 in barcode_map:
       well_id = barcode_map[bc1]
   ```
   BC2 is parsed, written into the read header as `BC1:xxxxxxxx BC2:yyyyyyyy`, and then never used. **No read has ever been discarded for barcode disagreement.** Enforcing `BC1 == BC2` would drop 4.50–8.18 % of assigned reads per sub-plate and remove 33 % of a failing well's contaminating molecules against 2.3 % of its genuine ones.
2. **There is no 3′ adapter trimming.** With 79 aligned bases and end-to-end alignment, any NlaIII fragment shorter than the read runs into the adapter and cannot align. That is **5.59 M fragments, 40 % of the digest**, structurally invisible.
3. **`bowtie2 -X/--maxins 500` is left at its default**, so any pair longer than 500 bp is never flagged properly paired. The effect is small here (~0.5 % of molecules) but it should be a conscious choice.
4. **`check_gc_rds` is dead in the DAG** — its shell references a `{params.gc_rds}` that `params:` never defines, and nothing consumes its flag. **Leave it dead.**
5. **GC template generation is manual.** `AneuFinder` and `BSgenome.Hsapiens.UCSC.hg38` conflict in the current `renv`, so `ad_hoc_checks/generate_gc_corr.R` is run by hand and the resulting `.rds` must never be declared as a Snakemake-owned output.

## B6. Sequencing submission

Recorded per plate in `PM25186_*.xlsx`. For the plates 21/22 run (`23GK3YLT3`) the eight libraries were sequenced in **one lane** with **combinatorial (non-unique) dual indices** — only 3 distinct i7 and 3 distinct i5 between them:

```
plate21_1 TAAGGCGA+GCGATCTA    plate22_1 TAAGGCGA+AGAGGATA
plate21_2 CGTACTAG+GCGATCTA    plate22_2 CGTACTAG+AGAGGATA
plate21_3 TAAGGCGA+ATAGAGAG    plate22_3 AGGCAGAA+ATAGAGAG
plate21_4 CGTACTAG+ATAGAGAG    plate22_4 AGGCAGAA+AGAGGATA
```

**This is the worst case for index hopping**: a single hop lands squarely on another real library. A4 measured the consequence (`02_FINDINGS.md` §5). Use **unique dual indices** on every future submission.

## B7. Where to find each per-well number

| you want | file | column |
|---|---|---|
| reads assigned per well | `<sub>/demux/demux_stats.json` | pairs per well |
| dimer rate | `<sub>/filtered/adapter_filter_summary.tsv` | `kept_pairs`, `dimer_rate` |
| mapped and unique reads | `<sub>/dedup/dedup_summary.tsv` | `total_reads`, `unique_reads`, `dedup_rate` |
| the gated metric | `<sub>/bam/{well}.stats.txt` | `SN reads mapped:` |
| all of the above, joined | `ad_hoc_checks/plate_qc_summary/output_plate21_22/well_metrics.csv` | one row per well |
| stage-by-stage retention | same directory | `stage_summary.csv` |
| per-well CATG concordance | `ad_hoc_checks/nlaIII_digest_check/tables/well_level_metrics.csv` | `pct_exact`, sites hit, fragments recovered |
| per-well contamination | `ad_hoc_checks/background_floor/tables/per_well_contamination.csv` | `c_measured`, `c_corrected`, genuine molecules, verdict |
| per-well insert size | `ad_hoc_checks/insert_size_vs_spri/tables/insert_density_by_well.csv.gz` | 1 bp resolution, 0–1000 bp |
| CellenONE per-well measurements | `<plate>/cellenone/wells_raw.tsv` | diameter, elongation, circularity, intensity, per-channel fluorescence, **and that well's own criteria** |
| human decisions | `<plate>/qc_decisions.csv` | `sample_id,well,decision,reason,notes` |
