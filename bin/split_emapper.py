#!/usr/bin/env python3
"""
Split a batched eggNOG-mapper .emapper.annotations file into one file per MAG.

Query IDs must look like "<mag_id>|<protein_id>"; the prefix is removed in the
per-MAG output. Every MAG id given on the command line gets a file (header
only if it had no hits). Comment lines (starting with '#') go to every file.
"""
import argparse
import os


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("annotations")
    p.add_argument("mag_ids", nargs="+")
    p.add_argument("--outdir", default="split")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    header, footer, rows = [], [], {m: [] for m in args.mag_ids}
    seen_data = False
    with open(args.annotations) as handle:
        for line in handle:
            if line.startswith("#"):
                (footer if seen_data else header).append(line)
                continue
            seen_data = True
            mag, sep, rest = line.partition("|")
            if not sep or mag not in rows:
                raise SystemExit(f"ERROR: unexpected query id in line: {line[:80]!r}")
            rows[mag].append(rest)

    for mag, lines in rows.items():
        out = os.path.join(args.outdir, f"{mag}.emapper.annotations")
        with open(out, "w") as fh:
            fh.writelines(header)
            fh.writelines(lines)
            fh.writelines(footer)


if __name__ == "__main__":
    main()
