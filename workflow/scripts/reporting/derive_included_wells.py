#!/usr/bin/env python3
"""
Derive the set of wells to include in the post-review (second-pass) AneuFinder run
from a human-edited, already-validated qc_decisions.csv (per-plate, in the data dir).

Inclusion semantics (see config/README.md):
  * PASS    -> included (always)
  * REVIEW  -> included only with --include-review (config qc_review.include_review)
  * EXCLUDE -> excluded
  * REPEAT  -> excluded, but logged ("flagged for rerun") so those wells are visible

Writes a TSV (`sample_id\twell`) of the included wells. This script is the gate that
turns human decisions into the second-pass input set; run validate_qc_decisions.py
first (the Snakemake rule depends on the validated flag).
"""

import argparse
import csv
import sys
from pathlib import Path

EXPECTED_HEADER = ["sample_id", "well", "decision", "reason", "notes"]


def parse_args():
    p = argparse.ArgumentParser(description="Derive included wells for the second AneuFinder pass.")
    p.add_argument("--decisions", required=True, help="Path to qc_decisions.csv")
    p.add_argument("--out", required=True, help="Output TSV of included wells")
    p.add_argument("--include-review", action="store_true",
                   help="Also include wells with decision REVIEW (default: PASS only)")
    return p.parse_args()


def main():
    args = parse_args()
    path = Path(args.decisions)
    if not path.exists():
        print(f"ERROR: decisions file not found: {path}", file=sys.stderr)
        sys.exit(1)

    included = []        # (sample_id, well)
    repeats = []
    counts = {"PASS": 0, "EXCLUDE": 0, "REVIEW": 0, "REPEAT": 0}

    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = [h.strip() for h in next(reader)]
        except StopIteration:
            print("ERROR: decisions file is empty", file=sys.stderr)
            sys.exit(1)
        if header != EXPECTED_HEADER:
            print(f"ERROR: header must be {EXPECTED_HEADER}, got {header}", file=sys.stderr)
            sys.exit(1)

        for row in reader:
            if not row or all(c.strip() == "" for c in row):
                continue
            if len(row) < 4:
                continue
            sample_id, well, decision = row[0].strip(), row[1].strip(), row[2].strip()
            if decision in counts:
                counts[decision] += 1
            if decision == "PASS":
                included.append((sample_id, well))
            elif decision == "REVIEW" and args.include_review:
                included.append((sample_id, well))
            elif decision == "REPEAT":
                repeats.append(well)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample_id", "well"])
        for sample_id, well in included:
            w.writerow([sample_id, well])

    print(f"Decisions: PASS={counts['PASS']} REVIEW={counts['REVIEW']} "
          f"EXCLUDE={counts['EXCLUDE']} REPEAT={counts['REPEAT']}")
    print(f"Included {len(included)} wells "
          f"({'PASS + REVIEW' if args.include_review else 'PASS only'}) -> {out}")
    if repeats:
        print(f"REPEAT (excluded, flagged for rerun): {', '.join(sorted(repeats))}")
    if not included:
        # Not an error here — let the second-pass rule decide. But make it loud.
        print("WARNING: no wells included for the second AneuFinder pass.", file=sys.stderr)


if __name__ == "__main__":
    main()
