#!/usr/bin/env python3
"""
Derive the set of wells to include in the post-review (second-pass) AneuFinder run
from a human-edited, already-validated qc_decisions.csv (per-plate, in the data dir).

Inclusion semantics (see config/README.md):
  * PASS    -> included (always)
  * REVIEW  -> included only with --include-review (config qc_review.include_review)
  * EXCLUDE -> excluded
  * REPEAT  -> excluded, but logged ("flagged for rerun") so those wells are visible

Writes a TSV of the included wells. This script is the gate that turns human
decisions into the second-pass input set; run validate_qc_decisions.py first (the
Snakemake rule depends on the validated flag).

The OUTPUT SCHEMA MIRRORS THE INPUT:
    sample_id,well,…            ->  `sample_id\twell`             (byte-identical to before)
    sample_id,subplate,well,…   ->  `sample_id\tsubplate\twell`   (384: a well id is
                                    only unique within its subplate)
"""

import argparse
import csv
import sys
from pathlib import Path

# Kept identical to validate_qc_decisions.py's copy, deliberately: the two scripts
# must never disagree about what a given decisions file means.
BASIC_HEADER = ["sample_id", "well", "decision", "reason", "notes"]
SUBPLATE_HEADER = ["sample_id", "subplate", "well", "decision", "reason", "notes"]


def detect_schema(header):
    """Identify the decisions-CSV schema; returns (name, {column: index})."""
    header = [h.strip() for h in header]
    if header == BASIC_HEADER:
        return "basic", {c: i for i, c in enumerate(BASIC_HEADER)}
    if header == SUBPLATE_HEADER:
        return "subplate", {c: i for i, c in enumerate(SUBPLATE_HEADER)}
    raise ValueError(
        f"header must be {BASIC_HEADER} or {SUBPLATE_HEADER}, got {header}")


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

    included = []        # (sample_id, subplate|None, well)
    repeats = []
    counts = {"PASS": 0, "EXCLUDE": 0, "REVIEW": 0, "REPEAT": 0}

    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            print("ERROR: decisions file is empty", file=sys.stderr)
            sys.exit(1)
        try:
            schema, idx = detect_schema(header)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)
        n_cols = len(SUBPLATE_HEADER) if schema == "subplate" else len(BASIC_HEADER)

        for row in reader:
            if not row or all(c.strip() == "" for c in row):
                continue
            if len(row) < n_cols - 1:      # notes may be absent/empty on the last column
                continue
            cells = [c.strip() for c in row]
            sample_id = cells[idx["sample_id"]]
            well = cells[idx["well"]]
            decision = cells[idx["decision"]]
            subplate = cells[idx["subplate"]] if schema == "subplate" else None
            if decision in counts:
                counts[decision] += 1
            if decision == "PASS" or (decision == "REVIEW" and args.include_review):
                included.append((sample_id, subplate, well))
            elif decision == "REPEAT":
                repeats.append(f"{subplate}/{well}" if subplate else well)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")  # unix EOL (no trailing \r)
        if schema == "subplate":
            w.writerow(["sample_id", "subplate", "well"])
            for sample_id, subplate, well in included:
                w.writerow([sample_id, subplate, well])
        else:
            w.writerow(["sample_id", "well"])
            for sample_id, _subplate, well in included:
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
