#!/usr/bin/env python3
"""
Detect objects in the CellenONE drop images, apply the ejection-line QC rule, and
render the per-well composite JPEGs the review report shows.

Why detection is needed at all
------------------------------
CellenONE's own tables record exactly ONE Transmission row per printed drop — the
cell it decided to isolate. They therefore cannot answer the question a reviewer
actually has: *how many objects were in this drop, and where were they?* That has to
come from the pixels. The tables remain valuable as ground truth for the one object
they do describe, which is what calibrates the pixel scale below.

Detection
---------
1. **Median background** over evenly-spaced `_Trans_*.png` of the run. This is
   markedly better than CellenONE's own `_Background.png`, which leaves static
   artifacts at x ~ 668 and x ~ 692 in every well and would score them as objects.
2. `dark = median - image > threshold`, then close(5x5) -> fill holes -> open(3x3)
   to knit each cell into one blob and drop single-pixel noise; label; equivalent
   diameter `2*sqrt(area/pi)`.
3. **Pixel scale** by regressing detected equivalent diameter (px) on the table's
   `Diameter` (um), pairing each image to the table row by NEAREST CENTROID.
   Pairing by "largest component" is WRONG and corrupts the scale: in 5 of 30 sampled
   wells the largest blob is not the isolated cell (e.g. M-3, where the table says
   X=382 but the largest component sits at X=718).

The line rule
-------------
The nozzle tip is at the LEFT and cells sediment leftwards toward it, so an object
far to the right is still up the capillary and will not be ejected.

    exactly 1 object                                    -> SINGLE
    >=2 objects, rightmost x > purple                   -> PASS
    rightmost between green and purple                  -> CONTAMINATION
    all objects at x <= green                           -> FAIL
    0 objects                                           -> NO_OBJECT

Any object down to the run's own detection minimum counts.

DISPLAY ONLY: nothing here feeds the automated status or the default decision.
"""

import argparse
import csv
import json
import logging
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

from collections import Counter

import numpy as np
import yaml
from PIL import Image, ImageDraw, ImageFont
from scipy import ndimage as ndi

CHANNELS = ("trans", "blue", "green", "orange", "red")


def setup_logging():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(levelname)s - %(message)s",
                        handlers=[logging.StreamHandler(sys.stdout)])
    return logging.getLogger("cellenone_images")


def parse_args():
    p = argparse.ArgumentParser(description="Render CellenONE cell images + image QC call.")
    p.add_argument("--cellenone-dir", required=True,
                   help="<plate>/cellenone (holds wells_raw.tsv + run_meta.json)")
    p.add_argument("--config", required=True, help="Pipeline config.yaml")
    p.add_argument("--threads", type=int, default=4)
    return p.parse_args()


# ---------------------------------------------------------------------------
# Detection
# ---------------------------------------------------------------------------

def conform(a: np.ndarray, shape) -> np.ndarray:
    """Crop or edge-pad a frame to the run's canonical geometry.

    The camera ROI can drift by a row or two part-way through a run: plate22 has 380
    frames at 473x952, 3 at 471 and 1 at 474. Refusing to render the whole plate over
    a 2-row difference is the wrong trade, and every frame shares the same origin, so
    trimming or padding the bottom edge leaves x (which the ejection-line rule uses)
    and y untouched.
    """
    h, w = shape
    out = a[:min(a.shape[0], h), :min(a.shape[1], w)]
    ph, pw = h - out.shape[0], w - out.shape[1]
    if ph > 0 or pw > 0:
        out = np.pad(out, ((0, max(0, ph)), (0, max(0, pw))), mode="edge")
    return out


def build_median_background(image_dir: Path, names, sample: int, logger):
    """Median over evenly-spaced Transmission frames.

    Every frame has at most a couple of small dark objects, so the per-pixel median
    across a few dozen of them is the empty capillary — including the static
    artifacts CellenONE's own background image fails to remove.
    """
    names = sorted(n for n in names if n)
    if not names:
        raise SystemExit("ERROR: no Transmission images to build a background from")
    if len(names) > sample:
        step = len(names) / sample
        names = [names[int(i * step)] for i in range(sample)]
    stack = []
    for n in names:
        try:
            stack.append(np.asarray(Image.open(image_dir / n).convert("L"), dtype=np.uint8))
        except OSError:
            continue
    if not stack:
        raise SystemExit("ERROR: could not read any Transmission image")
    counts = Counter(a.shape for a in stack)
    shape = counts.most_common(1)[0][0]
    if len(counts) > 1:
        logger.info("Transmission frames vary in size %s; conforming to the modal %s",
                    dict(counts), shape)
    stack = [conform(a, shape) for a in stack]
    logger.info("Median background from %d frames, shape %s", len(stack), shape)
    return np.median(np.stack(stack), axis=0).astype(np.float32), shape


def detect_objects(img: np.ndarray, background: np.ndarray, threshold: float,
                   min_area_px: int = 4):
    """Return [{x, y, area_px, diameter_px}], brightest-background minus image."""
    # Safety net: this is the one place a frame meets the run background, so conform
    # here too rather than trusting every caller to have done it. A mismatched ROI
    # would otherwise raise a broadcast error and abort the whole plate.
    if img.shape != background.shape:
        img = conform(img, background.shape)
    dark = (background - img.astype(np.float32)) > threshold
    dark = ndi.binary_closing(dark, np.ones((5, 5), bool))
    dark = ndi.binary_fill_holes(dark)
    dark = ndi.binary_opening(dark, np.ones((3, 3), bool))
    labels, n = ndi.label(dark)
    if n == 0:
        return []
    out = []
    objects = ndi.find_objects(labels)
    for i, sl in enumerate(objects, start=1):
        if sl is None:
            continue
        area = int((labels[sl] == i).sum())
        if area < min_area_px:
            continue
        cy, cx = ndi.center_of_mass(labels == i)
        out.append({"x": float(cx), "y": float(cy), "area_px": area,
                    "diameter_px": float(2.0 * np.sqrt(area / np.pi))})
    out.sort(key=lambda o: o["x"])
    return out


def image_call(objects, green_x: float, purple_x: float) -> str:
    """The line rule. See the module docstring for the reasoning."""
    if not objects:
        return "NO_OBJECT"
    if len(objects) == 1:
        return "SINGLE"
    rightmost = max(o["x"] for o in objects)
    if rightmost > purple_x:
        return "PASS"
    if rightmost > green_x:
        return "CONTAMINATION"
    return "FAIL"


# ---------------------------------------------------------------------------
# Rendering
# ---------------------------------------------------------------------------

def stretch(img: np.ndarray, lo_pct=1.0, hi_pct=99.5) -> np.ndarray:
    """Percentile contrast stretch to 0..1 (the raw frames are low-contrast)."""
    a = img.astype(np.float32)
    lo, hi = np.percentile(a, lo_pct), np.percentile(a, hi_pct)
    if hi <= lo:
        return np.zeros_like(a)
    return np.clip((a - lo) / (hi - lo), 0, 1)


def channel_signal(img: np.ndarray, background: np.ndarray, scale: float,
                   shift=(0.0, 0.0), gamma: float = 1.0):
    """Background-subtract a fluorescence frame and normalise it to 0..1.

    The window must come from the RUN, not from the isolation criteria and not from
    the frame itself:

    * `IsoIntMinFlu`/`IsoIntMaxFlu` (10/255) are the thresholds CellenONE uses to
      decide whether a cell is positive — not a display range. Real signal here sits
      ~40-80 against a ~20-30 background, so stretching over 10..255 renders every
      fluorescence frame almost black. That was the original bug.
    * Auto-scaling each frame to its own max is worse: a well with no cell would have
      its sensor noise stretched to full brightness and look like a bright cell.

    So `background` is the per-run median frame and `scale` a high percentile of the
    residual across the run: a well with signal is bright, an empty well stays black,
    and two wells are directly comparable by eye.

    Returns the SIGNAL only — colouring happens in `overlay_channel()`, because every
    view now paints the stain onto the transmission frame rather than showing it on
    black.

    `shift` is the per-run (dx, dy) optical offset measured at ingest time; a whole-
    pixel roll is enough, no feature-based registration is needed.
    """
    a = np.clip(img.astype(np.float32) - background, 0.0, None)
    a = np.clip(a / max(float(scale), 1e-6), 0.0, 1.0)
    if gamma and gamma != 1.0:
        # Lifts faint signal without clipping the peak, so a dim cell is still
        # obviously coloured rather than a barely-tinted smudge.
        a = np.power(a, float(gamma))
    dx, dy = int(round(shift[0])), int(round(shift[1]))
    if dx or dy:
        # The fluorescence frame sits at trans + offset, so shift it BACK by -offset.
        a = np.roll(np.roll(a, -dy, axis=0), -dx, axis=1)
    return a


def build_flu_reference(image_dir: Path, rows, sample: int, logger, shape=None):
    """Per-channel (median background, display scale) for the fluorescence channels.

    A channel is only included when the run actually recorded signal in it. On
    plate21 the Red LED was never configured — 0 of 384 wells have a non-zero Red
    intensity and the frames are flat noise — yet stray Red PNGs exist for 27 wells.
    Rendering those produced a uniform red wash that also swamped the merge, so an
    unused channel is dropped outright rather than shown as an empty overlay.
    """
    ref = {}
    for ch in ("blue", "green", "orange", "red"):
        signal = [r for r in rows
                  if (r.get(f"{ch}_intensity") or "").strip() not in ("", "0", "0.0", "0.00")]
        names = sorted({r[f"img_{ch}"] for r in rows if r.get(f"img_{ch}")})
        if not names:
            continue
        if not signal:
            logger.info("Channel %s: 0/%d wells have signal — not rendered "
                        "(unused LED on this run)", ch, len(rows))
            continue
        step = max(1, len(names) // sample)
        picked = names[::step][:sample]
        stack = []
        for n in picked:
            try:
                stack.append(np.asarray(Image.open(image_dir / n).convert("L"), dtype=np.float32))
            except OSError:
                continue
        if not stack:
            continue
        ref_shape = shape or Counter(a.shape for a in stack).most_common(1)[0][0]
        stack = [conform(a, ref_shape) for a in stack]
        bg = np.median(np.stack(stack), axis=0).astype(np.float32)
        # High percentile of the residual across sampled frames: bright enough that a
        # real cell saturates, high enough that background noise stays dark.
        tops = [float(np.percentile(np.clip(a - bg, 0, None), 99.9)) for a in stack]
        scale = max(float(np.percentile(tops, 95)), 8.0)
        ref[ch] = (bg, scale)
        logger.info("Channel %s: %d/%d wells with signal, background mean %.1f, "
                    "display scale %.1f", ch, len(signal), len(rows), bg.mean(), scale)
    return ref


def overlay_channel(base_rgb: np.ndarray, signal: np.ndarray, tint) -> np.ndarray:
    """Paint one fluorescence channel onto the transmission frame in its own colour.

    Alpha-composite toward the full-saturation tint (`out = base*(1-a) + tint*a`)
    rather than a screen blend. Screening over a bright transmission background washes
    the hue out to white exactly where the signal is strongest — the cell — which is
    the one place the colour has to read. With alpha compositing a strong pixel becomes
    the pure tint and a weak one stays transmission, so "is there stain on this cell"
    is answerable at a glance.
    """
    a = np.clip(signal, 0.0, 1.0)[..., None]
    colour = np.asarray(tint, dtype=np.float32)[None, None, :]
    return base_rgb * (1.0 - a) + colour * a


def overlay_channels(base_rgb: np.ndarray, signals: dict, tints: dict) -> np.ndarray:
    """Paint several channels at once, mixing their colours where they coincide.

    Compositing the channels one after another would let the last one painted hide the
    others: a cell positive in both Blue and Orange would render pure Orange and the
    Blue would simply be gone. Instead the per-channel signals set both the total
    opacity and the mix, so a double-positive cell reads as a blend and a
    single-positive one keeps its own hue.
    """
    if not signals:
        return base_rgb
    total = np.zeros(base_rgb.shape[:2], dtype=np.float32)
    weighted = np.zeros_like(base_rgb)
    for ch, sig in signals.items():
        a = np.clip(sig, 0.0, 1.0).astype(np.float32)
        total += a
        weighted += a[..., None] * np.asarray(tints.get(ch, (1.0, 1.0, 1.0)),
                                              dtype=np.float32)[None, None, :]
    mix = weighted / np.maximum(total, 1e-6)[..., None]   # colour blend, not brightness
    alpha = np.clip(total, 0.0, 1.0)[..., None]           # opacity from combined signal
    return base_rgb * (1.0 - alpha) + mix * alpha


_FONT_CACHE = {}


def _font(size: int):
    """A legible bitmap font, falling back through what Pillow can offer."""
    if size in _FONT_CACHE:
        return _FONT_CACHE[size]
    f = None
    for path in ("/usr/share/fonts/dejavu/DejaVuSansMono-Bold.ttf",
                 "/usr/share/fonts/truetype/dejavu/DejaVuSansMono-Bold.ttf"):
        try:
            f = ImageFont.truetype(path, size)
            break
        except OSError:
            continue
    if f is None:
        try:
            f = ImageFont.load_default(size)      # Pillow >= 10.1
        except TypeError:
            f = ImageFont.load_default()
    _FONT_CACHE[size] = f
    return f


def draw_badge(im: Image.Image, text: str):
    """Intensity readout, top-right, on a translucent grey plate.

    Drawn after the resize so it is a fixed on-screen size, and composited through an
    RGBA layer so the plate is genuinely see-through — the capillary wall behind it
    stays visible instead of being blanked out.
    """
    if not text:
        return im
    layer = Image.new("RGBA", im.size, (0, 0, 0, 0))
    d = ImageDraw.Draw(layer)
    font = _font(max(11, int(im.size[1] * 0.075)))
    try:
        x0, y0, x1, y1 = d.textbbox((0, 0), text, font=font)
        tw, th = x1 - x0, y1 - y0
    except AttributeError:
        tw, th = d.textsize(text, font=font)
    pad, margin = 5, 6
    bx1, by0 = im.size[0] - margin, margin
    bx0, by1 = bx1 - (tw + 2 * pad), by0 + (th + 2 * pad)
    d.rectangle([bx0, by0, bx1, by1], fill=(28, 30, 32, 150))
    d.text((bx0 + pad, by0 + pad - 1), text, font=font, fill=(255, 255, 255, 235))
    return Image.alpha_composite(im.convert("RGBA"), layer).convert("RGB")


def draw_lines(draw: ImageDraw.ImageDraw, green_x, purple_x, height, scale):
    for x, colour in ((green_x, (40, 200, 60)), (purple_x, (240, 120, 230))):
        xs = x * scale
        draw.line([(xs, 0), (xs, height)], fill=colour, width=max(1, int(2 * scale)))


# CellenONE's own particle taxonomy (cellenREPORT S4.3 "Isolation report"), reused so
# our outlines read the same way as the vendor's scatter plots:
#   isolated  green   the object actually dispensed into this well
#   fitting   yellow  met the morphology gate but was not the one dispensed
#   detected  pink    detected above the minimum, outside the isolation window
# NB the vendor's gate also uses elongation and intensity; we measure equivalent
# diameter only, so "fitting" here is a diameter-only approximation of theirs.
CATEGORY_COLOUR = {
    "isolated": (60, 220, 90),
    "fitting": (255, 210, 40),
    "detected": (255, 105, 180),
}


def annotate_objects(draw: ImageDraw.ImageDraw, objects, scale, categories=None):
    for i, o in enumerate(objects, start=1):
        cat = (categories or {}).get(i, "detected")
        colour = CATEGORY_COLOUR.get(cat, CATEGORY_COLOUR["detected"])
        r = max(4.0, o["diameter_px"] * 0.75) * scale
        x, y = o["x"] * scale, o["y"] * scale
        draw.ellipse([x - r, y - r, x + r, y + r], outline=colour,
                     width=max(1, int(2 * scale)))
        draw.text((x + r + 2, y - r), str(i), fill=colour)


# Shared read-only worker state. Set once in the parent before the pool is created:
# on Linux the default `fork` start method gives every child a copy-on-write view, so
# the ~1.8 MB background array is never pickled (384 tasks would otherwise ship it
# 384 times).
_CTX = {}


def _worker_init(image_dir, outdir, bg, meta, cfg, flu_ref=None):
    _CTX.update(image_dir=Path(image_dir), outdir=Path(outdir),
                bg=bg, meta=meta, cfg=cfg, flu_ref=flu_ref or {})


def render_well(row):
    """Worker: detect + render one well. Returns (well_row, object_rows)."""
    image_dir, outdir = _CTX["image_dir"], _CTX["outdir"]
    bg, meta, cfg = _CTX["bg"], _CTX["meta"], _CTX["cfg"]
    lines = meta["lines"]
    green_x, purple_x = float(lines["green_x"]), float(lines["purple_x"])
    crit = meta.get("criteria", {}) or {}
    offsets = meta.get("channel_offsets", {}) or {}

    img_cfg = cfg.get("image", {}) or {}
    scale = float(img_cfg.get("scale", 0.5))
    quality = int(img_cfg.get("quality", 72))
    want = [c for c in (img_cfg.get("channels") or ["merge", "trans"])]
    annotate = bool(img_cfg.get("annotate", True))
    tints = cfg.get("tint", {}) or {}
    threshold = float(cfg.get("dark_threshold", 8))
    # <1 lifts faint stain so the colour reads; the peak still saturates.
    flu_gamma = float(cfg.get("flu_gamma", 0.55))

    def well_crit(key, fallback):
        """This well's own gate, falling back to the run summary.

        The criteria columns are per ROW in the run table, and an operator can change
        them mid-run, so a single run-level number is wrong for whichever wells were
        dispensed under the other setting."""
        v = row.get(key)
        try:
            f = float(v)
        except (TypeError, ValueError):
            return fallback
        return f if f > 0 else fallback

    min_dia_um = cfg.get("min_object_diameter_um")
    if min_dia_um is None:
        min_dia_um = well_crit("det_dia_min", crit.get("DetDiaMinTrans", 5.0))
    px_per_um = float(meta.get("px_per_um") or 1.434)

    out = {"id": row["id"], "subplate": row["subplate"], "well": row["well"],
           "pos384": row["pos384"], "n_objects": 0, "rightmost_x": "",
           "rightmost_diameter_um": "", "n_in_iso_window": 0,
           "image_call": "NO_IMAGE", "image_paths": "",
           # CellenONE's own measurements of the cell it chose to isolate, passed
           # straight through from the run table. These describe the isolated cell;
           # the n_objects/rightmost_* fields above describe what OUR detection found
           # in the whole frame. Both are useful and they are not the same thing.
           "diameter_um": row.get("diameter_um", ""),
           "elongation": row.get("elongation", ""),
           "circularity": row.get("circularity", ""),
           "intensity": row.get("intensity", ""),
           "blue_intensity": row.get("blue_intensity", ""),
           "green_intensity": row.get("green_intensity", ""),
           "orange_intensity": row.get("orange_intensity", ""),
           "red_intensity": row.get("red_intensity", "")}
    objs_out = []

    trans_name = row.get("img_trans") or ""
    if not trans_name or not (image_dir / trans_name).exists():
        return out, objs_out

    try:
        trans = np.asarray(Image.open(image_dir / trans_name).convert("L"), dtype=np.uint8)
    except OSError:
        return out, objs_out
    trans = conform(trans, bg.shape)

    objects = detect_objects(trans, bg, threshold)
    # Any object down to the run's own detection minimum counts for the line rule.
    min_dia_px = float(min_dia_um) * px_per_um
    objects = [o for o in objects if o["diameter_px"] >= min_dia_px]

    iso_lo = well_crit("iso_dia_min", crit.get("IsoDiaMinTrans"))
    iso_hi = well_crit("iso_dia_max", crit.get("IsoDiaMaxTrans"))

    # Which object did CellenONE actually dispense? Its (X, Y) for this well is in the
    # run table, so the nearest detected centroid is the isolated cell. Same pairing
    # rule as the pixel-scale calibration, and for the same reason: the dispensed cell
    # is often not the largest blob in the frame.
    iso_index = None
    try:
        tx, ty = float(row.get("x")), float(row.get("y"))
    except (TypeError, ValueError):
        tx = ty = None
    if tx is not None and objects:
        best = min(range(len(objects)),
                   key=lambda k: (objects[k]["x"] - tx) ** 2 + (objects[k]["y"] - ty) ** 2)
        d = ((objects[best]["x"] - tx) ** 2 + (objects[best]["y"] - ty) ** 2) ** 0.5
        if d <= 20:                      # no confident pairing -> claim nothing
            iso_index = best + 1

    n_in_window = 0
    categories = {}
    for i, o in enumerate(objects, start=1):
        dia_um = o["diameter_px"] / px_per_um
        in_window = (iso_lo is not None and iso_hi is not None
                     and iso_lo <= dia_um <= iso_hi)
        n_in_window += int(in_window)
        cat = "isolated" if i == iso_index else ("fitting" if in_window else "detected")
        categories[i] = cat
        objs_out.append({
            "id": row["id"], "object": i, "x": round(o["x"], 1), "y": round(o["y"], 1),
            "area_px": o["area_px"], "diameter_px": round(o["diameter_px"], 2),
            "diameter_um": round(dia_um, 2), "in_iso_window": int(in_window),
            "category": cat,
            "side": ("right_of_purple" if o["x"] > purple_x else
                     "between_lines" if o["x"] > green_x else "left_of_green"),
        })

    call = image_call(objects, green_x, purple_x)

    # ---- render -----------------------------------------------------------
    h, w = trans.shape
    ow, oh = max(1, int(w * scale)), max(1, int(h * scale))
    base = stretch(trans)
    rgb_trans = np.repeat(base[..., None], 3, axis=2)

    flu = {}
    # NB: do NOT unpack into `bg`/`scale` here — those names already hold the
    # transmission background and the image downscale factor, and `scale` is still
    # needed below for draw_lines()/annotate_objects(). Shadowing it silently drew the
    # ejection lines at x*44 instead of x*0.5, i.e. off-canvas.
    for ch, (flu_bg, flu_scale) in (_CTX.get("flu_ref") or {}).items():
        name = row.get(f"img_{ch}") or ""
        if not name or not (image_dir / name).exists():
            continue
        try:
            a = np.asarray(Image.open(image_dir / name).convert("L"), dtype=np.uint8)
        except OSError:
            continue
        a = conform(a, trans.shape)
        if flu_bg.shape != trans.shape:
            continue
        off = offsets.get(ch, {}) or {}
        flu[ch] = channel_signal(a, flu_bg, flu_scale,
                                 shift=(off.get("dx", 0.0), off.get("dy", 0.0)),
                                 gamma=flu_gamma)

    written = []
    outdir.mkdir(parents=True, exist_ok=True)

    def save(kind, arr, badge=""):
        im = Image.fromarray((np.clip(arr, 0, 1) * 255).astype(np.uint8), "RGB")
        if (ow, oh) != (w, h):
            im = im.resize((ow, oh), Image.LANCZOS)
        d = ImageDraw.Draw(im)
        draw_lines(d, green_x, purple_x, oh, scale)
        if annotate and kind in ("merge", "trans"):
            annotate_objects(d, objects, scale, categories)
        if badge:
            im = draw_badge(im, badge)
        name = f"{row['id']}_{kind}.jpg"
        im.save(outdir / name, "JPEG", quality=quality, optimize=True)
        written.append(name)

    def intensity(ch):
        v = row.get(f"{ch}_intensity")
        try:
            f = float(v)
        except (TypeError, ValueError):
            return None
        return f if f > 0 else None      # 0 == nothing recorded, not a measurement

    for kind in want:
        if kind == "trans":
            save("trans", rgb_trans)
        elif kind == "merge":
            # Every channel painted onto the transmission frame, in its own colour.
            merged = overlay_channels(rgb_trans, flu, tints)
            parts = [f"{c[:1].upper()} {intensity(c):.0f}" for c in sorted(flu)
                     if intensity(c) is not None]
            save("merge", merged, badge="  ".join(parts))
        elif kind in flu:
            # A single channel is ALSO shown over transmission, not on black: the
            # reviewer needs to see which object in the capillary is stained, and a
            # bare fluorescence frame gives no anatomical context to place it in.
            single = overlay_channel(rgb_trans, flu[kind], tints.get(kind, (1.0, 1.0, 1.0)))
            iv = intensity(kind)
            save(kind, single,
                 badge=f"{kind.upper()} {iv:.1f}" if iv is not None else f"{kind.upper()} —")

    out.update({
        "n_objects": len(objects),
        "rightmost_x": round(max((o["x"] for o in objects), default=0), 1) if objects else "",
        "rightmost_diameter_um": (
            round(max(objects, key=lambda o: o["x"])["diameter_px"] / px_per_um, 2)
            if objects else ""),
        "n_in_iso_window": n_in_window,
        "image_call": call,
        "image_paths": ";".join(written),
    })
    return out, objs_out


# ---------------------------------------------------------------------------
# Pixel-scale calibration
# ---------------------------------------------------------------------------

def calibrate_px_per_um(rows, image_dir: Path, bg, threshold, logger, n_sample=40):
    """Regress detected equivalent diameter (px) on the table's Diameter (um).

    Pairing is by NEAREST CENTROID to the table's (X, Y), never by largest component:
    the isolated cell is often not the biggest blob in the frame, and calibrating on
    the wrong blob silently corrupts every diameter downstream.
    """
    xs, ys = [], []
    used = 0
    for row in rows:
        if used >= n_sample:
            break
        name = row.get("img_trans") or ""
        dia_um = row.get("diameter_um") or ""
        try:
            dia_um = float(dia_um)
            tx, ty = float(row["x"]), float(row["y"])
        except (TypeError, ValueError):
            continue
        if not name or not (image_dir / name).exists() or dia_um <= 0:
            continue
        try:
            img = np.asarray(Image.open(image_dir / name).convert("L"), dtype=np.uint8)
        except OSError:
            continue
        # Same ROI drift as everywhere else — conform before differencing against the
        # background, or a plate whose camera shifted mid-run cannot be calibrated.
        img = conform(img, bg.shape)
        objects = detect_objects(img, bg, threshold)
        if not objects:
            continue
        nearest = min(objects, key=lambda o: (o["x"] - tx) ** 2 + (o["y"] - ty) ** 2)
        if abs(nearest["x"] - tx) > 20 or abs(nearest["y"] - ty) > 20:
            continue          # no confident pairing; do not let it bias the fit
        xs.append(dia_um)
        ys.append(nearest["diameter_px"])
        used += 1

    if len(xs) < 5:
        logger.warning("Only %d calibration pairs; falling back to px_per_um=1.434", len(xs))
        return 1.434, len(xs)
    xs_a, ys_a = np.asarray(xs), np.asarray(ys)
    # Through the origin: a 0 um object is 0 px.
    scale = float((xs_a * ys_a).sum() / (xs_a * xs_a).sum())
    logger.info("Pixel scale: %.3f px/um (from %d nearest-centroid pairs)", scale, len(xs))
    return scale, len(xs)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    logger = setup_logging()
    args = parse_args()
    cdir = Path(args.cellenone_dir)

    meta = json.loads((cdir / "run_meta.json").read_text())
    with open(cdir / "wells_raw.tsv", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    logger.info("Loaded %d well records", len(rows))

    with open(args.config) as fh:
        cfg = (yaml.safe_load(fh) or {}).get("cellenone", {}) or {}

    image_dir = Path(meta.get("image_dir") or meta["run_dir"])
    threshold = float(cfg.get("dark_threshold", 8))

    bg, frame_shape = build_median_background(
        image_dir, [r.get("img_trans") for r in rows],
        int(cfg.get("median_sample", 60)), logger)

    px_per_um, n_cal = calibrate_px_per_um(rows, image_dir, bg, threshold, logger)
    meta["px_per_um"] = px_per_um
    meta["px_per_um_pairs"] = n_cal

    outdir = cdir / "images"
    outdir.mkdir(parents=True, exist_ok=True)
    flu_ref = build_flu_reference(image_dir, rows, int(cfg.get("median_sample", 60)),
                                  logger, shape=frame_shape)
    meta["flu_channels"] = sorted(flu_ref)
    init_args = (str(image_dir), str(outdir), bg, meta, cfg, flu_ref)

    well_rows, object_rows = [], []
    if args.threads > 1:
        with ProcessPoolExecutor(max_workers=args.threads,
                                 initializer=_worker_init, initargs=init_args) as ex:
            for wr, orows in ex.map(render_well, rows, chunksize=4):
                well_rows.append(wr)
                object_rows.extend(orows)
    else:
        _worker_init(*init_args)
        for r in rows:
            wr, orows = render_well(r)
            well_rows.append(wr)
            object_rows.extend(orows)

    well_fields = ["id", "subplate", "well", "pos384", "n_objects", "rightmost_x",
                   "rightmost_diameter_um", "n_in_iso_window", "image_call", "image_paths",
                   "diameter_um", "elongation", "circularity", "intensity",
                   "blue_intensity", "green_intensity", "orange_intensity",
                   "red_intensity"]
    with open(cdir / "cellenone_wells.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=well_fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(well_rows)

    obj_fields = ["id", "object", "x", "y", "area_px", "diameter_px", "diameter_um",
                  "in_iso_window", "category", "side"]
    with open(cdir / "objects.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=obj_fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(object_rows)

    counts = {}
    for r in well_rows:
        counts[r["image_call"]] = counts.get(r["image_call"], 0) + 1
    meta["image_calls"] = counts
    (cdir / "run_meta.json").write_text(json.dumps(meta, indent=2) + "\n")

    logger.info("Image calls: %s",
                ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    logger.info("Wrote %d well rows, %d object rows, images in %s",
                len(well_rows), len(object_rows), outdir)


if __name__ == "__main__":
    main()
