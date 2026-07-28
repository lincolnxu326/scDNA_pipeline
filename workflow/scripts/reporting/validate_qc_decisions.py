#!/usr/bin/env python3
"""
Validate a human-edited qc_decisions.csv (per-plate, in the data dir) against the well set
and controlled vocabularies.

Run on demand (not part of `rule all`). On success it writes a small flag file;
on any problem it prints a per-row report and exits non-zero so a downstream
rule can depend on a validated decisions file.

CSV contract (see config/README.md). Two schemas are accepted:
    sample_id,well,decision,reason,notes                    (96-well / legacy)
    sample_id,subplate,well,decision,reason,notes           (384-well)
Both remain valid indefinitely; a 96-well plate keeps producing the 5-column form.
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
BASIC_HEADER = ["sample_id", "well", "decision", "reason", "notes"]
SUBPLATE_HEADER = ["sample_id", "subplate", "well", "decision", "reason", "notes"]


def detect_schema(header):
    """Identify the decisions-CSV schema.

    Returns ("basic"|"subplate", {column_name: index}). Shared verbatim with
    derive_included_wells.py so the two can never disagree about a file.
    """
    header = [h.strip() for h in header]
    if header == BASIC_HEADER:
        return "basic", {c: i for i, c in enumerate(BASIC_HEADER)}
    if header == SUBPLATE_HEADER:
        return "subplate", {c: i for i, c in enumerate(SUBPLATE_HEADER)}
    raise ValueError(
        f"header must be {BASIC_HEADER} or {SUBPLATE_HEADER}, got {header}")


def parse_args():
    p = argparse.ArgumentParser(description="Validate a per-plate qc_decisions.csv.")
    p.add_argument("--decisions", required=True, help="Path to qc_decisions.csv")
    p.add_argument("--wells", required=True, nargs="+", help="Allowed well IDs")
    p.add_argument("--subplates", nargs="*", default=[],
                   help="Allowed subplate names (384 mode). Empty => not checked.")
    p.add_argument("--out-flag", required=False, help="Flag file written on success")
    p.add_argument("--require-all-wells", action="store_true",
                   help="Error if any pipeline well is missing from the CSV")
    return p.parse_args()


def main():
    args = parse_args()
    allowed_wells = set(args.wells)
    allowed_subplates = set(args.subplates or [])
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

        try:
            schema, idx = detect_schema(header)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            sys.exit(1)
        n_cols = len(SUBPLATE_HEADER) if schema == "subplate" else len(BASIC_HEADER)

        seen = set()
        n_rows = 0
        for lineno, row in enumerate(reader, start=2):
            if not row or all(c.strip() == "" for c in row):
                continue
            n_rows += 1
            if len(row) != n_cols:
                errors.append(f"line {lineno}: expected {n_cols} columns, got {len(row)}")
                continue
            cells = [c.strip() for c in row]
            well = cells[idx["well"]]
            decision = cells[idx["decision"]]
            reason = cells[idx["reason"]]
            subplate = cells[idx["subplate"]] if schema == "subplate" else None

            if well not in allowed_wells:
                errors.append(f"line {lineno}: unknown well '{well}'")
            if subplate is not None and allowed_subplates and subplate not in allowed_subplates:
                errors.append(f"line {lineno}: unknown subplate '{subplate}'")
            # A well id is only unique WITHIN a subplate, so 384 keys on the pair.
            key = (subplate, well) if schema == "subplate" else well
            if key in seen:
                errors.append(f"line {lineno}: duplicate well '{well}'"
                              + (f" in subplate '{subplate}'" if subplate else ""))
            seen.add(key)

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
        if schema == "subplate":
            subs = sorted(allowed_subplates) or sorted({s for s, _ in seen if s})
            expected = {(s, w) for s in subs for w in allowed_wells}
            missing = expected - seen
            labels = sorted(f"{s}/{w}" for s, w in missing)
        else:
            missing = allowed_wells - seen
            labels = sorted(missing)
        if missing:
            errors.append(f"missing decisions for {len(missing)} wells: "
                          f"{', '.join(labels)}")

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
