#!/usr/bin/env python3
"""
Validate a human-edited qc_decisions.csv (per-plate, in the data dir) against the well set
and controlled vocabularies.

Run on demand (not part of `rule all`). On success it writes a small flag file;
on any problem it prints a per-row report and exits non-zero so a downstream
rule can depend on a validated decisions file.

CSV contract (see config/README.md):
    sample_id,well,decision,reason,notes
"""

import argparse
import csv
import sys
from pathlib import Path

DECISIONS = {"PASS", "EXCLUDE", "REVIEW", "REPEAT"}
# `reason` may be empty, or a ';'-separated list of these tokens.
REASONS = {
    "low_read_count",
    "noisy_profile",
    "poor_bin_distribution",
    "low_complexity",
    "suspected_doublet_or_mixed_well",
    "sample_swap_suspected",
    "manual_exception",
    "other",
    "missing_qc_metric",
}
EXPECTED_HEADER = ["sample_id", "well", "decision", "reason", "notes"]


def parse_args():
    p = argparse.ArgumentParser(description="Validate a per-plate qc_decisions.csv.")
    p.add_argument("--decisions", required=True, help="Path to qc_decisions.csv")
    p.add_argument("--wells", required=True, nargs="+", help="Allowed well IDs")
    p.add_argument("--out-flag", required=False, help="Flag file written on success")
    p.add_argument("--require-all-wells", action="store_true",
                   help="Error if any pipeline well is missing from the CSV")
    return p.parse_args()


def main():
    args = parse_args()
    allowed_wells = set(args.wells)
    path = Path(args.decisions)
    errors = []

    if not path.exists():
        print(f"ERROR: decisions file not found: {path}", file=sys.stderr)
        sys.exit(1)

    with open(path, newline="") as fh:
        reader = csv.reader(fh)
        try:
            header = next(reader)
        except StopIteration:
            print("ERROR: decisions file is empty", file=sys.stderr)
            sys.exit(1)

        header = [h.strip() for h in header]
        if header != EXPECTED_HEADER:
            print(f"ERROR: header must be {EXPECTED_HEADER}, got {header}", file=sys.stderr)
            sys.exit(1)

        seen = set()
        n_rows = 0
        for lineno, row in enumerate(reader, start=2):
            if not row or all(c.strip() == "" for c in row):
                continue
            n_rows += 1
            if len(row) != 5:
                errors.append(f"line {lineno}: expected 5 columns, got {len(row)}")
                continue
            sample_id, well, decision, reason, _notes = [c.strip() for c in row]

            if well not in allowed_wells:
                errors.append(f"line {lineno}: unknown well '{well}'")
            if well in seen:
                errors.append(f"line {lineno}: duplicate well '{well}'")
            seen.add(well)

            if decision not in DECISIONS:
                errors.append(
                    f"line {lineno}: invalid decision '{decision}' "
                    f"(allowed: {', '.join(sorted(DECISIONS))})"
                )
            # reason may be empty, or a ';'-separated list of allowed tokens
            if reason:
                for tok in reason.split(";"):
                    tok = tok.strip()
                    if tok and tok not in REASONS:
                        errors.append(
                            f"line {lineno}: invalid reason '{tok}' "
                            f"(allowed: {', '.join(sorted(REASONS))})"
                        )

    if args.require_all_wells:
        missing = allowed_wells - seen
        if missing:
            errors.append(f"missing decisions for {len(missing)} wells: "
                          f"{', '.join(sorted(missing))}")

    if errors:
        print(f"VALIDATION FAILED for {path} ({len(errors)} problem(s)):", file=sys.stderr)
        for e in errors:
            print(f"  - {e}", file=sys.stderr)
        sys.exit(1)

    print(f"OK: {path} validated ({n_rows} rows, {len(seen)} wells).")
    if args.out_flag:
        flag = Path(args.out_flag)
        flag.parent.mkdir(parents=True, exist_ok=True)
        flag.write_text(f"validated {n_rows} rows\n")
        print(f"Wrote flag: {flag}")


if __name__ == "__main__":
    main()
