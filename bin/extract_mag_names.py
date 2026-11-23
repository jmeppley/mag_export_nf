#!/usr/bin/env python3
"""
Extract unique MAG fasta names from one column of a table and check which
exist in a MAG directory.

Writes:
  mags.txt          names found in mag_dir (first-seen order, de-duplicated)
  missing_mags.txt  names listed in the table but absent from mag_dir

Exits non-zero if --on-missing fail and any names are missing.
"""
import argparse
import csv
import os
import sys


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("table")
    p.add_argument("--column", required=True,
                   help="column name (with --header) or 1-based column index")
    p.add_argument("--header", action="store_true",
                   help="first line of the table is a header")
    p.add_argument("--sep", default="\t", help="field separator (default: tab)")
    p.add_argument("--mag-dir", required=True)
    p.add_argument("--mag-ext", default="",
                   help="suffix to append to names that lack it (e.g. .fa)")
    p.add_argument("--on-missing", choices=["warn", "fail"], default="warn")
    p.add_argument("--found-out", default="mags.txt")
    p.add_argument("--missing-out", default="missing_mags.txt")
    return p.parse_args()


def column_index(header_row, column, has_header):
    """Return the 0-based column index for a name or 1-based index."""
    if has_header and column in header_row:
        return header_row.index(column)
    try:
        idx = int(column)
    except ValueError:
        if has_header:
            sys.exit(f"ERROR: column '{column}' not in header: {header_row}")
        sys.exit(f"ERROR: table has no header, so column must be a 1-based "
                 f"index, got '{column}'")
    if idx < 1:
        sys.exit(f"ERROR: column index must be >= 1, got {idx}")
    return idx - 1


def main():
    args = parse_args()
    sep = args.sep.encode().decode("unicode_escape")  # allow "\t" from config

    names = []
    seen = set()
    with open(args.table, newline="") as handle:
        reader = csv.reader(handle, delimiter=sep)
        col = None
        for lineno, row in enumerate(reader, start=1):
            if not row or not any(cell.strip() for cell in row):
                continue
            if col is None:
                col = column_index(row, args.column, args.header)
                if args.header:
                    continue
            if col >= len(row):
                sys.exit(f"ERROR: {args.table} line {lineno} has {len(row)} "
                         f"columns; need column {col + 1}")
            name = row[col].strip()
            if not name:
                continue
            if args.mag_ext and not name.endswith(args.mag_ext):
                name += args.mag_ext
            if name not in seen:
                seen.add(name)
                names.append(name)

    found = [n for n in names if os.path.isfile(os.path.join(args.mag_dir, n))]
    missing = [n for n in names if n not in set(found)]

    with open(args.found_out, "w") as out:
        out.writelines(n + "\n" for n in found)
    with open(args.missing_out, "w") as out:
        out.writelines(n + "\n" for n in missing)

    print(f"{args.table}: {len(names)} unique names, {len(found)} found, "
          f"{len(missing)} missing", file=sys.stderr)
    if missing and args.on_missing == "fail":
        shown = ", ".join(missing[:10]) + (" ..." if len(missing) > 10 else "")
        sys.exit(f"ERROR: {len(missing)} MAG(s) from {args.table} not found in "
                 f"{args.mag_dir}: {shown}")


if __name__ == "__main__":
    main()
