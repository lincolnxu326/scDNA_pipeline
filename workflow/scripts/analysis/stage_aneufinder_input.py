#!/usr/bin/env python3
"""
Build the symlink-only input directory that run_aneufinder.R globs.

Both AneuFinder passes need the same thing: a directory of `<id>.bam` (+ `.bam.bai`)
symlinks pointing at the real per-well dedup BAMs. This script is that staging step
for both, and for both plate formats.

Why it exists
-------------
1. **Dead wells abort the batch.** AneuFinder aborts the ENTIRE run if any single
   well has ~no reads (GC-loess yields "NAs found in reads"). Dead wells sit ~1000x
   below the read-count PASS gate (auto-FAIL) and carry no usable CN signal, so
   `--min-reads` drops them and lets the run complete.
2. **384 mode needs namespaced ids.** run_aneufinder.R derives the model id from the
   staged BAM's basename, and render_well_profiles.R derives the plot name from the
   model's `.RData` basename. Staging `plate21_2_W07.bam` is therefore the whole
   mechanism by which plate-level output is correctly namespaced — with zero R edits.

**Symlinks only. No BAM is ever copied, moved or mutated.**

Examples
--------
    # first pass, 96-well plate
    stage_aneufinder_input.py --bam-dir plate21_1=<plate>/bam \
        --id-template '{well}' --out-dir <plate>/aneufinder/input_bams

    # first pass, 384-well plate
    stage_aneufinder_input.py \
        --bam-dir plate21_1=<plate>/plate21_1/bam ... --bam-dir plate21_4=... \
        --id-template '{sub}_{well}' --out-dir <plate>/aneufinder/input_bams

    # second (post-review) pass — the well set comes from included_wells.tsv
    stage_aneufinder_input.py --bam-dir ... --wells-file <plate>/qc_review/included_wells.tsv \
        --min-reads 0 --out-dir <plate>/aneufinder_reviewed/input_bams
"""

import argparse
import os
import re
import sys
from pathlib import Path

FLAGSTAT_MAPPED = re.compile(r"^(\d+)\s.*\smapped \(")


def parse_args():
    p = argparse.ArgumentParser(
        description="Stage a symlink-only AneuFinder input directory.")
    p.add_argument("--bam-dir", action="append", required=True, metavar="SUB=PATH",
                   help="Per-subplate dedup BAM directory, e.g. plate21_1=/…/plate21_1/bam. "
                        "Repeatable; in 96 mode pass it once.")
    p.add_argument("--id-template", default="{well}",
                   help="Staged BAM basename template over {sub} and {well} "
                        "[default: %(default)s]")
    p.add_argument("--wells-file", default="",
                   help="included_wells.tsv restricting the well set (2- or 3-column). "
                        "Default: every W*.bam found in the BAM dirs.")
    p.add_argument("--min-reads", type=int, default=100,
                   help="Skip wells with fewer mapped reads in their flagstat; "
                        "0 disables the filter [default: %(default)s]")
    p.add_argument("--out-dir", required=True, help="Directory of symlinks to create")
    p.add_argument("--report", default="", help="Optional audit TSV listing every well considered")
    p.add_argument("--label", default="AneuFinder input",
                   help="Prefix for the summary line [default: %(default)s]")
    return p.parse_args()


def parse_bam_dirs(specs) -> dict:
    """['sub=path', …] -> {sub: Path}, order preserved."""
    out = {}
    for spec in specs:
        if "=" not in spec:
            sys.exit(f"ERROR: --bam-dir must be SUB=PATH, got {spec!r}")
        sub, path = spec.split("=", 1)
        sub, path = sub.strip(), path.strip()
        if not sub or not path:
            sys.exit(f"ERROR: --bam-dir must be SUB=PATH, got {spec!r}")
        if sub in out:
            sys.exit(f"ERROR: --bam-dir given twice for subplate {sub!r}")
        out[sub] = Path(path)
    return out


def mapped_reads(bam_dir: Path, well: str):
    """Mapped-read count from `<bam_dir>/<well>.flagstat.txt`, or None if unreadable.

    Same `/ mapped \\(/` parse the inline shell loop used, so the kept/excluded set
    is unchanged for existing plates.
    """
    fs = bam_dir / f"{well}.flagstat.txt"
    try:
        with open(fs) as fh:
            for line in fh:
                m = FLAGSTAT_MAPPED.match(line)
                if m:
                    return int(m.group(1))
    except OSError:
        return None
    return None


def read_wells_file(path: Path):
    """Read an included_wells.tsv; return [(subplate_or_None, well), …].

    Accepts both schemas:
        sample_id, well              (96 mode / legacy)
        sample_id, subplate, well    (384 mode)
    Tolerates CRLF (a CSV round-tripped through Excel on Windows).
    """
    rows = []
    with open(path, newline="") as fh:
        header = [h.strip().lstrip("\ufeff") for h in fh.readline().rstrip("\r\n").split("\t")]
        if "well" in header:
            widx = header.index("well")
            sidx = header.index("subplate") if "subplate" in header else None
        else:
            # headerless or unexpected: fall back to positional, and keep this line
            widx, sidx = (2, 1) if len(header) >= 3 else (1, None)
            fh.seek(0)
        for line in fh:
            parts = [c.strip() for c in line.rstrip("\r\n").split("\t")]
            if not parts or len(parts) <= widx or not parts[widx]:
                continue
            if parts[widx] == "well":      # header line reached via the fallback path
                continue
            sub = parts[sidx] if (sidx is not None and len(parts) > sidx) else None
            rows.append((sub or None, parts[widx]))
    return rows


def clear_symlinks(out_dir: Path) -> int:
    """Remove existing symlinks so a re-run cannot inherit a stale well set.

    Only symlinks are removed — a regular file in the staging directory is left
    alone and reported, because deleting one could destroy real data.
    """
    removed = 0
    if not out_dir.exists():
        return 0
    for entry in out_dir.iterdir():
        if entry.is_symlink():
            entry.unlink()
            removed += 1
        elif entry.is_file():
            print(f"WARNING: leaving non-symlink file in staging dir: {entry}", file=sys.stderr)
    return removed


def main():
    args = parse_args()
    bam_dirs = parse_bam_dirs(args.bam_dir)
    out_dir = Path(args.out_dir)

    # ---- work list --------------------------------------------------------
    work = []          # (sub, well)
    if args.wells_file:
        wf = Path(args.wells_file)
        if not wf.exists():
            sys.exit(f"ERROR: wells file not found: {wf}")
        for sub, well in read_wells_file(wf):
            if sub is None:
                if len(bam_dirs) != 1:
                    sys.exit(
                        f"ERROR: {wf} has no `subplate` column but {len(bam_dirs)} BAM "
                        f"directories were given; the well set is ambiguous. Re-derive the "
                        f"included wells from a 6-column decisions CSV."
                    )
                sub = next(iter(bam_dirs))
            if sub not in bam_dirs:
                sys.exit(f"ERROR: {wf} refers to unknown subplate {sub!r}; "
                         f"known: {', '.join(bam_dirs)}")
            work.append((sub, well))
    else:
        for sub, d in bam_dirs.items():
            if not d.is_dir():
                sys.exit(f"ERROR: BAM directory not found: {d}")
            for bam in sorted(d.glob("W*.bam")):
                work.append((sub, bam.stem))

    if not work:
        sys.exit(f"ERROR: no wells to stage (bam dirs: {', '.join(str(d) for d in bam_dirs.values())})")

    # ---- stage ------------------------------------------------------------
    out_dir.mkdir(parents=True, exist_ok=True)
    clear_symlinks(out_dir)

    kept, excluded, missing, audit = 0, [], [], []
    for sub, well in work:
        bam_dir = bam_dirs[sub]
        bam = bam_dir / f"{well}.bam"
        stage_id = args.id_template.format(sub=sub, well=well)
        if not bam.exists():
            missing.append(f"{sub}/{well}")
            audit.append((stage_id, sub, well, str(bam), "", "missing_bam"))
            continue
        n = mapped_reads(bam_dir, well)
        if args.min_reads > 0 and (n is None or n < args.min_reads):
            excluded.append(f"{well}({n if n is not None else 0})" if len(bam_dirs) == 1
                            else f"{sub}/{well}({n if n is not None else 0})")
            audit.append((stage_id, sub, well, str(bam), "" if n is None else n, "excluded_low_reads"))
            continue

        # Absolute but NOT symlink-resolved, matching the `ln -sf "$abs_path"` this
        # replaced: /nemo and /camp are two mount views of the same tree and resolving
        # would silently rewrite one into the other.
        link = out_dir / f"{stage_id}.bam"
        if link.is_symlink() or link.exists():
            link.unlink()
        link.symlink_to(os.path.abspath(bam))
        bai = Path(f"{bam}.bai")
        if bai.exists():
            bai_link = out_dir / f"{stage_id}.bam.bai"
            if bai_link.is_symlink() or bai_link.exists():
                bai_link.unlink()
            bai_link.symlink_to(os.path.abspath(bai))
        kept += 1
        audit.append((stage_id, sub, well, str(bam), "" if n is None else n, "kept"))

    # ---- report -----------------------------------------------------------
    # Wording kept close to the inline shell loop this replaced, so existing logs
    # stay greppable.
    print(f"{args.label}: kept {kept} wells"
          + (f" (>= {args.min_reads} mapped)" if args.min_reads > 0 else "")
          + f"; excluded dead wells:{''.join(' ' + e for e in excluded)}")
    if missing:
        print(f"WARNING: {len(missing)} requested well(s) had no BAM: {', '.join(missing)}",
              file=sys.stderr)

    if args.report:
        rp = Path(args.report)
        rp.parent.mkdir(parents=True, exist_ok=True)
        with open(rp, "w") as fh:
            fh.write("id\tsubplate\twell\tbam\tmapped_reads\tstatus\n")
            for r in audit:
                fh.write("\t".join(str(x) for x in r) + "\n")
        print(f"Staging audit: {rp}")

    if kept == 0:
        print(f"ERROR: nothing staged into {out_dir}; AneuFinder has no input.", file=sys.stderr)
        sys.exit(1)
    print(f"Staged {kept} BAM symlink(s) into {out_dir}")


if __name__ == "__main__":
    main()
