#!/usr/bin/env python3
"""
384-well plate layout: the mapping between a 384 plate position (A1..P24) and the
(subplate, well) pair the sequencing pipeline actually processes.

A 384-well plate is dispensed as four interleaved 96-well subplates. The pipeline
processes each subplate as its own FASTQ pair (`plate21_1` .. `plate21_4`), so
384 mode has to translate between the two coordinate systems constantly.

The mapping is DETERMINISTIC — no lookup table is needed. For a 384 position with
row index ``r`` (0..15 -> A..P) and column index ``c`` (0..23 -> 1..24):

    subplate_index = (r % 2) + 2 * (c // 12)      # 0 -> SL1 ... 3 -> SL4
    well_number    = (c % 12) * 8 + (r // 2) + 1  # W01..W96

i.e. SL1 = odd rows (A,C,E,...) x cols 1-12, SL2 = even rows x cols 1-12,
SL3 = odd rows x cols 13-24, SL4 = even rows x cols 13-24, and within a subplate
the fill is column-major down its 8 rows. Verified against all 384 rows of
`plate_mapping_384.xlsx` (384/384 match) — see ``--check-xlsx``.

An override TSV (`plate_layout_tsv` in config.yaml) supports future non-standard
layouts; when given it replaces the formula entirely.

IMPORTANT: this module is imported at Snakefile PARSE time as well as by the
reporting scripts, which run in different conda envs. It must therefore stay
**stdlib-only**.
"""

import argparse
import re
import sys
from pathlib import Path

# 384-well geometry
N_ROWS = 16          # A..P
N_COLS = 24          # 1..24
N_SUBPLATES = 4
SUB_ROWS = 8         # rows per subplate
SUB_COLS = 12        # cols per subplate
WELLS_PER_SUB = SUB_ROWS * SUB_COLS   # 96

# "A-1" (CellenONE image/table token) or "A1" / "A01" (plate map spreadsheets)
_POS_RE = re.compile(r"^\s*([A-Pa-p])\s*-?\s*0*([0-9]{1,2})\s*$")
# "W01".."W96"
_WELL_RE = re.compile(r"^\s*[Ww]0*([0-9]{1,3})\s*$")
# "<plate>_<n>" subplate directory names
_SUB_RE = re.compile(r"^(?P<stem>.+)_(?P<idx>[0-9]+)$")


# ---------------------------------------------------------------------------
# Position parsing / formatting
# ---------------------------------------------------------------------------

def parse_pos(token: str):
    """Parse 'A1' / 'A01' / 'A-1' / 'P-24' into a 0-based (row, col).

    Raises ValueError on anything outside the 16x24 grid.
    """
    m = _POS_RE.match(str(token))
    if not m:
        raise ValueError(f"not a 384 well position: {token!r}")
    row = ord(m.group(1).upper()) - ord("A")
    col = int(m.group(2)) - 1
    if not (0 <= row < N_ROWS and 0 <= col < N_COLS):
        raise ValueError(f"position out of range for a 384 plate: {token!r}")
    return row, col


def format_pos(row: int, col: int, sep: str = "") -> str:
    """(0,0) -> 'A1'; sep='-' gives the CellenONE form 'A-1'."""
    return f"{chr(ord('A') + row)}{sep}{col + 1}"


def parse_well(token: str) -> int:
    """'W01' -> 1. Raises ValueError outside 1..96."""
    m = _WELL_RE.match(str(token))
    if not m:
        raise ValueError(f"not a well id: {token!r}")
    n = int(m.group(1))
    if not (1 <= n <= WELLS_PER_SUB):
        raise ValueError(f"well number out of range 1..{WELLS_PER_SUB}: {token!r}")
    return n


def format_well(n: int) -> str:
    """1 -> 'W01'."""
    return f"W{n:02d}"


# ---------------------------------------------------------------------------
# The formula, both directions
# ---------------------------------------------------------------------------

def pos_to_subwell(row: int, col: int):
    """(row, col) on the 384 plate -> (sub_index 0..3, well_id 'W01'..'W96')."""
    if not (0 <= row < N_ROWS and 0 <= col < N_COLS):
        raise ValueError(f"position out of range: ({row}, {col})")
    sub_index = (row % 2) + 2 * (col // SUB_COLS)
    well_number = (col % SUB_COLS) * SUB_ROWS + (row // 2) + 1
    return sub_index, format_well(well_number)


def subwell_to_pos(sub_index: int, well):
    """(sub_index 0..3, 'W01'..'W96') -> (row, col) on the 384 plate."""
    if not (0 <= sub_index < N_SUBPLATES):
        raise ValueError(f"subplate index out of range 0..3: {sub_index}")
    n = parse_well(well) if isinstance(well, str) else int(well)
    if not (1 <= n <= WELLS_PER_SUB):
        raise ValueError(f"well number out of range 1..{WELLS_PER_SUB}: {well!r}")
    n0 = n - 1
    row = 2 * (n0 % SUB_ROWS) + (sub_index % 2)
    col = SUB_COLS * (sub_index // 2) + (n0 // SUB_ROWS)
    return row, col


# ---------------------------------------------------------------------------
# Subplate discovery
# ---------------------------------------------------------------------------

def _fastq_pair_present(d: Path) -> bool:
    """True when <d>/<d.name>_R1.fastq.gz and _R2.fastq.gz both exist.

    The plate naming convention requires the FASTQ basename to equal the
    directory name, so this is an exact check rather than a glob.
    """
    return (d / f"{d.name}_R1.fastq.gz").exists() and (d / f"{d.name}_R2.fastq.gz").exists()


def discover_subplates(plate_dir, expected=None) -> list:
    """Return the ordered subplate directory names under a 384 plate directory.

    Looks for `<plate>_<n>` subdirectories that carry a complete FASTQ pair and
    sorts them by the trailing integer. `expected` (config `subplates:`) short-
    circuits discovery but is still validated to exist.
    """
    plate_dir = Path(plate_dir)
    plate = plate_dir.name

    if expected:
        subs = [str(s) for s in expected]
        missing = [s for s in subs if not (plate_dir / s).is_dir()]
        if missing:
            raise FileNotFoundError(
                f"config `subplates` lists directories that do not exist under {plate_dir}: "
                + ", ".join(missing)
            )
        return subs

    found = []
    for child in sorted(plate_dir.iterdir() if plate_dir.is_dir() else []):
        if not child.is_dir():
            continue
        m = _SUB_RE.match(child.name)
        if not m or m.group("stem") != plate:
            continue
        if not _fastq_pair_present(child):
            continue
        found.append((int(m.group("idx")), child.name))

    if not found:
        raise FileNotFoundError(
            f"No subplate directories found under {plate_dir}.\n"
            f"384 mode expects `{plate}_1` .. `{plate}_N` subdirectories, each holding\n"
            f"`<subplate>_R1.fastq.gz` + `<subplate>_R2.fastq.gz`.\n"
            f"If this is a single 96-well plate, run with plate_format=96 instead."
        )

    found.sort(key=lambda t: t[0])
    return [name for _, name in found]


# ---------------------------------------------------------------------------
# Full layout
# ---------------------------------------------------------------------------

def layout_384(subplates, wells=None, plate=None) -> list:
    """Build the full 384-position layout.

    Returns a list of dicts, one per (subplate, well), each with:
        id         '<subplate>_<well>'  (== sample_id; the AneuFinder model id)
        subplate   'plate21_2'
        sub_index  0-based index into `subplates`
        sub_label  'SL2'
        well       'W07'
        pos384     'A1'
        row, col   0-based grid coordinates on the 384 plate
        sample_id  same as id (kept explicit; the HTML/CSV use this name)

    `wells` defaults to W01..W96 in order. Ordering of the result follows
    (subplate, well), i.e. the order the pipeline processes them in.
    """
    if wells is None:
        wells = [format_well(i) for i in range(1, WELLS_PER_SUB + 1)]
    if len(subplates) > N_SUBPLATES:
        raise ValueError(
            f"a 384 plate has at most {N_SUBPLATES} subplates, got {len(subplates)}: {subplates}"
        )

    out = []
    for si, sub in enumerate(subplates):
        for w in wells:
            row, col = subwell_to_pos(si, w)
            out.append({
                "id": f"{sub}_{w}",
                "subplate": sub,
                "sub_index": si,
                "sub_label": f"SL{si + 1}",
                "well": w,
                "pos384": format_pos(row, col),
                "row": row,
                "col": col,
                "sample_id": f"{sub}_{w}",
            })
    return out


def load_layout_tsv(path, subplates=None) -> list:
    """Load an override layout TSV (columns: subplate, well, pos384).

    Used for future non-standard 384 layouts. Column order does not matter; a
    `sub_index` column is honoured if present, otherwise the index comes from
    `subplates` (or first-seen order).
    """
    path = Path(path)
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        header = [h.strip() for h in header]
        need = ["subplate", "well", "pos384"]
        missing = [c for c in need if c not in header]
        if missing:
            raise ValueError(
                f"{path}: layout TSV needs columns {need}; missing {missing} (got {header})"
            )
        idx = {c: header.index(c) for c in header}
        for lineno, line in enumerate(fh, start=2):
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(need):
                raise ValueError(f"{path}:{lineno}: expected {len(header)} columns, got {len(parts)}")
            sub = parts[idx["subplate"]].strip()
            well = parts[idx["well"]].strip()
            pos = parts[idx["pos384"]].strip()
            row, col = parse_pos(pos)
            if "sub_index" in idx and parts[idx["sub_index"]].strip():
                si = int(parts[idx["sub_index"]])
            elif subplates and sub in subplates:
                si = list(subplates).index(sub)
            else:
                si = None
            rows.append({"subplate": sub, "well": well, "pos384": format_pos(row, col),
                         "row": row, "col": col, "sub_index": si})

    order = list(subplates) if subplates else []
    for r in rows:
        if r["sub_index"] is None:
            if r["subplate"] not in order:
                order.append(r["subplate"])
            r["sub_index"] = order.index(r["subplate"])
        r["id"] = f"{r['subplate']}_{r['well']}"
        r["sample_id"] = r["id"]
        r["sub_label"] = f"SL{r['sub_index'] + 1}"
    return rows


def build_layout(subplates, wells=None, layout_tsv=None) -> list:
    """`layout_384()` unless an override TSV is configured."""
    if layout_tsv:
        return load_layout_tsv(layout_tsv, subplates=subplates)
    return layout_384(subplates, wells=wells)


# ---------------------------------------------------------------------------
# CLI (self-check)
# ---------------------------------------------------------------------------

def _check_roundtrip() -> int:
    bad = 0
    for r in range(N_ROWS):
        for c in range(N_COLS):
            si, w = pos_to_subwell(r, c)
            r2, c2 = subwell_to_pos(si, w)
            if (r, c) != (r2, c2):
                print(f"MISMATCH {format_pos(r, c)}: SL{si+1}/{w} -> {format_pos(r2, c2)}")
                bad += 1
    print(f"round-trip: {N_ROWS * N_COLS - bad}/{N_ROWS * N_COLS} OK, {bad} mismatch(es)")
    return bad


def _read_xlsx_rows(path):
    """Yield each sheet-1 row as a list of cell strings, using only the stdlib.

    An .xlsx is a zip of XML parts; we only need the shared-string table and the
    first worksheet's inline/shared text. Written out here rather than pulled from
    openpyxl so this module stays importable at Snakefile parse time.
    """
    import zipfile
    import xml.etree.ElementTree as ET

    NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

    def cell_text(node):
        return "".join(t.text or "" for t in node.iter(f"{NS}t"))

    with zipfile.ZipFile(path) as z:
        shared = []
        if "xl/sharedStrings.xml" in z.namelist():
            root = ET.fromstring(z.read("xl/sharedStrings.xml"))
            shared = [cell_text(si) for si in root.findall(f"{NS}si")]

        sheets = sorted(n for n in z.namelist()
                        if n.startswith("xl/worksheets/sheet") and n.endswith(".xml"))
        if not sheets:
            return
        root = ET.fromstring(z.read(sheets[0]))
        for row in root.iter(f"{NS}row"):
            values = []
            for c in row.findall(f"{NS}c"):
                if c.get("t") == "s":
                    v = c.find(f"{NS}v")
                    idx = int(v.text) if v is not None and v.text else -1
                    values.append(shared[idx] if 0 <= idx < len(shared) else "")
                elif c.get("t") == "inlineStr":
                    values.append(cell_text(c))
                else:
                    v = c.find(f"{NS}v")
                    values.append(v.text if v is not None and v.text else "")
            yield values


def _check_xlsx(path) -> int:
    """Compare the formula against a `plate_mapping_384.xlsx` (SL1 - W01 -> A1)."""
    sl_re = re.compile(r"SL\s*([1-4])\s*-\s*(W\s*0*[0-9]{1,2})", re.I)

    n_checked = bad = 0
    for row in _read_xlsx_rows(path):
        cells = [str(c).strip() for c in row if c is not None and str(c).strip()]
        sl_hit = pos_hit = None
        for c in cells:
            if sl_hit is None and sl_re.search(c):
                sl_hit = sl_re.search(c)
            elif pos_hit is None:
                try:
                    pos_hit = parse_pos(c)
                except ValueError:
                    pass
        if sl_hit is None or pos_hit is None:
            continue
        si = int(sl_hit.group(1)) - 1
        well = format_well(parse_well(sl_hit.group(2).replace(" ", "")))
        expect = subwell_to_pos(si, well)
        n_checked += 1
        if expect != pos_hit:
            bad += 1
            if bad <= 10:
                print(f"MISMATCH SL{si+1} - {well}: xlsx {format_pos(*pos_hit)} "
                      f"!= formula {format_pos(*expect)}")

    if n_checked == 0:
        print(f"ERROR: no 'SLn - Wnn' + position pairs found in {path}", file=sys.stderr)
        return 1
    print(f"{n_checked - bad}/{n_checked} {'OK' if bad == 0 else 'MATCH'}"
          + (f", {bad} mismatch(es)" if bad else ""))
    return bad


def main():
    p = argparse.ArgumentParser(description="384-well plate layout helper / self-check.")
    p.add_argument("--check-xlsx", metavar="XLSX",
                   help="Validate the formula against a plate_mapping_384.xlsx")
    p.add_argument("--discover", metavar="PLATE_DIR",
                   help="Print the discovered subplates for a 384 plate directory")
    p.add_argument("--pos", metavar="A1", help="Show the (subplate, well) for a 384 position")
    p.add_argument("--subwell", nargs=2, metavar=("SL", "WELL"),
                   help="Show the 384 position for e.g. `1 W08` (SL1/W08)")
    p.add_argument("--dump", metavar="PLATE_DIR",
                   help="Print the full layout TSV for a 384 plate directory")
    args = p.parse_args()

    rc = 0
    if args.check_xlsx:
        rc |= _check_roundtrip()
        rc |= _check_xlsx(args.check_xlsx)
    if args.discover:
        for s in discover_subplates(args.discover):
            print(s)
    if args.pos:
        r, c = parse_pos(args.pos)
        si, w = pos_to_subwell(r, c)
        print(f"{format_pos(r, c)}\tSL{si + 1}\t{w}")
    if args.subwell:
        si = int(args.subwell[0]) - 1
        r, c = subwell_to_pos(si, args.subwell[1])
        print(f"SL{si + 1}\t{args.subwell[1]}\t{format_pos(r, c)}")
    if args.dump:
        subs = discover_subplates(args.dump)
        print("id\tsubplate\tsub_index\tsub_label\twell\tpos384\trow\tcol")
        for e in layout_384(subs):
            print("\t".join(str(e[k]) for k in
                            ("id", "subplate", "sub_index", "sub_label", "well",
                             "pos384", "row", "col")))
    if not any([args.check_xlsx, args.discover, args.pos, args.subwell, args.dump]):
        rc |= _check_roundtrip()
    sys.exit(1 if rc else 0)


if __name__ == "__main__":
    main()
