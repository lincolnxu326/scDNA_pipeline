# Reading the QC review: gates, colours, and where every number comes from

This report mixes three sources that are easy to confuse:

1. numbers the **CellenONE instrument measured** and wrote into its run folder,
2. numbers **we derived** from the raw images because the instrument does not record
   them, and
3. **sequencing** numbers, which have nothing to do with CellenONE at all.

Everything below says which is which. If you only read one thing, read
[The two PASS/FAIL words](#the-two-passfail-words-that-are-not-the-same-thing).

---

## 1. The two PASS/FAIL words that are not the same thing

The report shows two independent verdicts per well. They share the words PASS and
FAIL and mean completely different things. **A well is routinely `FAIL` on one and
`PASS` on the other, and that is not a contradiction.**

| | **Gate status** | **Droplet call** (image call) |
|---|---|---|
| answers | "did this cell sequence well enough to use?" | "what was actually in the droplet?" |
| computed from | sequencing read counts | the CellenONE camera image |
| values | `PASS` / `WARN` / `FAIL` | `SINGLE` / `PASS` / `CONTAMINATION` / `FAIL` / `NO_OBJECT` |
| shown as | the plate-map colour, and the big reads figure | the badge on the cell image |
| **decides anything?** | **YES** — it is the only thing that pre-fills the decision | **NO** — display only |

### Gate status — the only thing that gates

Read-count only. `usable_reads` = mapped reads in the deduplicated BAM.

| usable_reads | gate status | decision pre-filled |
|---|---|---|
| ≥ 100,000 | `PASS` | PASS |
| 50,000 – 99,999 | `WARN` | REVIEW |
| < 50,000 | `FAIL` | EXCLUDE |

Cutoffs live in `config.yaml` under `qc_review.usable_reads_*_cutoff`. Duplication
rate is **not** gated — it is shown as a labelled diagnostic only.

### Droplet call — never gates anything

The nozzle is at the **left**. Cells sediment leftwards toward it, so an object far to
the right is still up the capillary and will not be ejected into the well.

| call | meaning |
|---|---|
| `SINGLE` | exactly one object — the ideal case |
| `PASS` | several objects, but the extras are right of the purple line and will not be ejected |
| `CONTAMINATION` | a second object sits between the lines and may be ejected too |
| `FAIL` | every object is already left of the ejection line |
| `NO_OBJECT` | nothing detected |
| `NO_IMAGE` | no CellenONE image for this well |

This call is **display-only by deliberate design.** It never touches the gate status
or the pre-filled decision. Reading "droplet PASS" as "this well passed QC" is the
single easiest mistake to make here.

---

## 2. The plate map

Each well is coloured by its **current decision**:

| colour | decision |
|---|---|
| green | PASS |
| red | EXCLUDE |
| blue | REVIEW |
| amber | REPEAT |

Wells are shown faint when they fall outside the subplate tab you have selected
(`All | SL1 | SL2 | SL3 | SL4`); dimmed wells are still clickable. Nothing else
changes a well's appearance.

A small white square in a well's corner marks progress: hollow = seen, filled = you
have set a decision.

> **Previously** the map also faded wells that were not in a "needs you" queue, which
> combined the read gate with the droplet call. That subset was removed: no metric has
> earned that authority yet, and while this is exploratory the report shows the plate
> as it is rather than steering you toward a subset. If a shading difference is
> visible in an older report, that is what it was.

---

## 3. What CellenONE gives us

The `.Run` folder is configured per run and **values genuinely change between runs —
and sometimes within one run.** Nothing below is hardcoded; it is all read from the
run's own files.

### 3.1 Per-well measurements — `Reordered_*_isolated.xls`

One row per (drop, channel). **Only ever describes the ONE cell the instrument chose
to isolate** — verified: 0 of 384 wells have more than one Transmission row.

| column | meaning | in the report as |
|---|---|---|
| `X`, `Y` | isolated cell's position in the 952×471 frame | used to identify which detected object is the isolated one |
| `Diameter` | its diameter (µm) | second metrics column, "Cell ø" |
| `Elongation` | 1.0 = perfectly round | second metrics column |
| `Circularity` | shape regularity | second metrics column |
| `Intensity` | transmission darkness | second metrics column |
| `Intensity` (Blue/Orange rows) | fluorescence per channel | badge on the fluorescence image |
| `Well` | plate position, e.g. `A-1` | mapped to `(subplate, well)` |
| `DropNo` | drop attempt index | run statistics |

### 3.2 The criteria — also in that file, **per row**

| column | plate21 | plate22 | meaning |
|---|---|---|---|
| `DetDiaMinTrans` | 5 µm | 5 µm | **detection gate** — smallest thing counted as an object |
| `IsoDiaMinTrans` / `IsoDiaMaxTrans` | 12–25 µm | **varies, see below** | **isolation gate** — what qualifies as a dispensable cell |
| `IsoIntMinFlu` / `IsoIntMaxFlu` | 10–255 | 10–255 | fluorescence *positivity* thresholds — **not a display range** |
| `EjBound` | 330 | 330 | ⚠ **not** the ejection line; do not use it |
| `SedBound` | 250 | 250 | sedimentation parameter |

> **⚠ The isolation gate can change mid-run.** plate22 was run with three different
> settings as the operator adjusted them:
>
> | IsoDiaMin | IsoDiaMax | wells |
> |---|---|---|
> | 9 µm | 25 µm | 20 |
> | 9 µm | 30 µm | 72 |
> | 12 µm | 30 µm | 292 |
>
> Each well is therefore judged against **its own row's criteria**, carried in
> `wells_raw.tsv` as `det_dia_min` / `iso_dia_min` / `iso_dia_max`. `run_meta.json`
> keeps a modal summary plus a `criteria_varied` record, and the ingest logs a warning
> when it happens. plate21 did not vary — a single 12–25 µm gate for all 384 wells.

### 3.3 Run-level — the `.par` file and the PDF report

| | |
|---|---|
| `Ejection (pix)` = 396 | the green line |
| `Sedimentation (pix)` = 250 | |
| `FluLEDsOrder` | which LEDs were used, by hardware index |
| Channel table | `Transmission: Positive` · `Blue: Inspect` · `Orange: Inspect` |
| §4.3 Isolation report | Detected 1220 · Fitting criteria 1071 · Isolated 384 · Non-printed 836 |

**"Inspect" means the channel was recorded but did NOT gate cell selection.** On both
plates only Transmission was `Positive`, so fluorescence never influenced which cells
the instrument picked.

> **The folder never says what the colours mean biologically.** It records the
> hardware channel (Blue / Orange) and its threshold — no fluorophore, filter or
> wavelength anywhere. That mapping lives in the experiment design and has to be
> supplied by whoever set the run up.

### 3.4 Other files

`Reordered_*_geoprops.xls` lists every detected object across all drop attempts —
used only for run-level statistics (isolation frequency) and as a cross-check.
`Tscatter.xls` / `Fscatter.xls` back the PDF's scatter plots.

---

## 4. What we derived

Everything here is computed by the pipeline from the raw images, because the
instrument does not record it.

### 4.1 Why detection is necessary at all

CellenONE's tables describe **one** object per well. They cannot answer "how many
objects were in this droplet, and where were they?" — which is the question the
droplet call needs. So we detect objects ourselves:

1. **Median background** over ~60 transmission frames of the run. Better than
   CellenONE's own `_Background.png`, which leaves static artifacts at x ≈ 668 and
   x ≈ 692 in every well.
2. `background − image > dark_threshold` (default 8), then morphological cleanup.
3. Equivalent diameter `2·√(area/π)`, converted to µm.

### 4.2 Pixel scale — derived, anchored to the instrument

`px_per_um` is fitted per run by regressing our detected diameter (px) against the
table's `Diameter` (µm), pairing each image to its table row **by nearest centroid**.

plate21: 1.407 px/µm · plate22: 1.392 px/µm.

> Pairing by "largest blob" would be wrong: the dispensed cell often is not the
> biggest object in the frame. In `plate21_1/W23` the real cell is at x=382
> (area 275) while a larger blob sits at x=718 (area 623).

### 4.3 Object categories — the circle colours

Each detected object is classified using CellenONE's own vocabulary
(cellenREPORT §4.3), so the words match the vendor report:

| colour | category | rule | source |
|---|---|---|---|
| **green** | `isolated` | our detected object nearest the table's `X,Y` (within 20 px) | position **from CellenONE**, matching **derived** |
| **yellow** | `fitting` | our diameter inside that well's isolation gate | threshold **from CellenONE**, measurement **derived** |
| **pink** | `detected` | anything else above the detection gate | threshold **from CellenONE**, measurement **derived** |

Written to `cellenone/objects.tsv` as the `category` column.

**Cross-check:** `isolated` comes to exactly **384** on both plates — one per well —
matching the PDF's "Isolated particles: 384". Our other counts are lower than the
PDF's because CellenONE counts across all drop attempts, including the non-printed
ones; we only image the 384 printed wells.

> **Caveat.** CellenONE's real "fitting criteria" gate uses diameter **and elongation
> and intensity**. We measure equivalent diameter only, so yellow-vs-pink is a
> **diameter-only approximation** of theirs.

### 4.4 Fluorescence display — derived

The images are background-subtracted per run and scaled to a high percentile of the
residual, then painted onto the transmission frame in colour with the measured
intensity in the corner.

> This deliberately does **not** use `IsoIntMinFlu`/`IsoIntMaxFlu` (10–255). Those are
> positivity thresholds, not a display range — real signal sits ~40–80 over a ~20–30
> background, so stretching over 10–255 renders every frame nearly black. Per-frame
> auto-scaling would be worse still: it turns an empty well's sensor noise into a
> convincing bright "cell".

A channel with no signal anywhere in the run is dropped rather than shown as an empty
overlay. On both plates Red was never used (0/384 wells), so it is not rendered even
though 27 stray Red PNGs exist.

### 4.5 Everything else we derive

| value | from |
|---|---|
| `gc_content` | weighted mean of the GC distribution samtools already writes |
| `mapping_rate` | the **raw** BAM — the dedup BAM has had unmapped reads removed, so it reads 100 % for every well |
| `n_objects`, `rightmost_x`, `n_in_iso_window` | our detection |
| droplet call | our detection + the two ejection lines |
| channel offsets | median fluorescence-minus-transmission position across the run |
| frame geometry | modal frame shape; the camera ROI can drift mid-run (plate22: 473×952 for 380 wells, 471 for 3, 474 for 1) |

---

## 5. Quick reference — which file holds what

| file | contents |
|---|---|
| `cellenone/wells_raw.tsv` | one row per well: CellenONE's measurements **and that well's own criteria** |
| `cellenone/objects.tsv` | one row per **detected** object: position, diameter, `in_iso_window`, `category` |
| `cellenone/cellenone_wells.tsv` | per-well summary the report consumes: counts, droplet call, image paths |
| `cellenone/run_meta.json` | lines and their source, criteria (+ `criteria_varied`), channel offsets, pixel scale, run statistics |
| `cellenone/images/` | rendered JPEGs, one per well per channel |
| `qc_review/review.html` | the report |
| `qc_decisions.csv` | your decisions — the reproducible artefact |

---

## 6. The three things worth remembering

1. **Gate status decides; the droplet call does not.** Only read counts pre-fill a
   decision.
2. **Plate-map colour is the decision, nothing else.** Faintness only means "outside
   the selected subplate tab".
3. **Instrument thresholds are read per run and per well, never assumed** — they
   changed three times inside plate22 alone.
