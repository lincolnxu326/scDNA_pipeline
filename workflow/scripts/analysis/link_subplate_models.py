#!/usr/bin/env python3
"""
Populate a 384 plate's AneuFinder MODELS directory from its subplates' existing models.

Why this is sound
-----------------
AneuFinder is strictly **per-cell**: every BAM is independently binned, GC-corrected
and segmented, and no cross-well information enters a well's model. A model computed
during a 96-well subplate run is therefore bit-for-bit the model a plate-level run
would produce. Only the batch-wide artefacts differ (the cluster PDF and the
genome-wide heatmap), and `render_well_profiles.R` regenerates those at plate level
from the `.RData` files anyway.

So for a 384 plate whose four subplates have already been processed, the plate-level
first pass is a naming exercise, not a compute one:

    plate21/aneufinder/MODELS/method-edivisive/plate21_1_W01.RData
       -> ../../../../plate21_1/aneufinder/MODELS/method-edivisive/W01.RData

Because `render_well_profiles.R` derives each plot's name from the `.RData` basename,
these links alone yield correctly namespaced plate-level plots — with zero AneuFinder
compute and zero R edits.

Links are RELATIVE, so the plate directory stays movable, and only symlinks are ever
created or removed; no model file is copied or modified.
"""

import argparse
import os
import sys
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Symlink per-subplate AneuFinder models into a plate-level MODELS dir.")
    p.add_argument("--plate-dir", required=True, help="384 plate directory")
    p.add_argument("--subplates", required=True, nargs="+", help="Subplate directory names")
    p.add_argument("--methods", required=True,
                   help="Comma-separated AneuFinder methods, e.g. 'edivisive'")
    p.add_argument("--outdir", required=True,
                   help="Plate-level aneufinder output dir (MODELS/ is created inside)")
    p.add_argument("--aneufinder-subdir", default="aneufinder",
                   help="Per-subplate AneuFinder dir name [default: %(default)s]")
    return p.parse_args()


def main():
    args = parse_args()
    plate_dir = Path(args.plate_dir)
    outdir = Path(args.outdir)
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    if not methods:
        sys.exit("ERROR: --methods is empty")

    total_linked = 0
    for method in methods:
        dest = outdir / "MODELS" / f"method-{method}"
        dest.mkdir(parents=True, exist_ok=True)

        # A re-run must not inherit models from a subplate set that has since changed.
        for entry in dest.iterdir():
            if entry.is_symlink():
                entry.unlink()

        linked = 0
        per_sub = []
        for sub in args.subplates:
            src_dir = plate_dir / sub / args.aneufinder_subdir / "MODELS" / f"method-{method}"
            if not src_dir.is_dir():
                sys.exit(
                    f"ERROR: no models for subplate {sub}, method {method}: {src_dir}\n"
                    f"Run that subplate through MODE=aneufinder_first first, or set "
                    f"`aneufinder.plate384_models: rerun` to compute plate-level models."
                )
            n = 0
            for rdata in sorted(src_dir.glob("*.RData")):
                # The staged id is what every downstream name derives from.
                link = dest / f"{sub}_{rdata.name}"
                if link.is_symlink() or link.exists():
                    link.unlink()
                link.symlink_to(os.path.relpath(rdata, dest))
                n += 1
            if n == 0:
                sys.exit(f"ERROR: {src_dir} contains no .RData models")
            per_sub.append(f"{sub}={n}")
            linked += n

        print(f"method-{method}: linked {linked} models ({', '.join(per_sub)}) -> {dest}")
        total_linked += linked

    if total_linked == 0:
        sys.exit("ERROR: nothing linked; the plate-level MODELS dir would be empty.")
    print(f"Reused {total_linked} existing per-subplate model(s); AneuFinder was not run.")


if __name__ == "__main__":
    main()
