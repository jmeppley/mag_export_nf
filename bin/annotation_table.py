#!/usr/bin/env python3
"""
Combine pyrodigal gene coordinates (GFF) with eggNOG-mapper annotations into
one tab-separated table, one row per predicted gene.

Every gene in the GFF gets a row; genes without an eggNOG hit have empty
annotation columns. eggNOG-mapper's "-" (no value) is written as empty.

Columns:
  gene_id         {contig_id}_{gene_no}  (matches the protein IDs in the .faa)
  start, end, strand
  OG              first entry of eggNOG_OGs (COG/KOG or 5-character eggNOG family)
  domain          Bacteria | Archaea | Eukaryota | Viruses (from eggNOG_OGs)
  max_annot_lvl   taxon name of eggNOG-mapper's max annotation level
  COG_category
  eggNOG_OGs      full column
  Description
  Preferred_name
  PFAMs
"""
import argparse
import csv
import sys

# NCBI taxids eggNOG uses for the top-level groups
DOMAINS = {"2": "Bacteria", "2157": "Archaea", "2759": "Eukaryota", "10239": "Viruses"}

OUT_COLUMNS = ["gene_id", "start", "end", "strand", "OG", "domain",
               "max_annot_lvl", "COG_category", "eggNOG_OGs", "Description",
               "Preferred_name", "PFAMs"]

EMAPPER_COLUMNS = ["eggNOG_OGs", "max_annot_lvl", "COG_category",
                   "Description", "Preferred_name", "PFAMs"]


def clean(value):
    value = (value or "").strip()
    return "" if value == "-" else value


def read_gff(path):
    """Yield (gene_id, start, end, strand) for CDS features, in file order."""
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "CDS":
                continue
            attrs = dict(a.split("=", 1) for a in fields[8].split(";") if "=" in a)
            # pyrodigal IDs are "{contig_index}_{gene_no}"
            gene_no = attrs["ID"].rsplit("_", 1)[-1]
            yield f"{fields[0]}_{gene_no}", fields[3], fields[4], fields[6]


def read_emapper(path):
    """Return {query: {column: value}} from an .emapper.annotations file."""
    hits, columns = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#query"):
                columns = line[1:].split("\t")
                continue
            if not line or line.startswith("#"):
                continue
            if columns is None:
                sys.exit(f"ERROR: no '#query' header line in {path}")
            row = dict(zip(columns, line.split("\t")))
            hits[row["query"]] = row
    if columns is not None:
        # match column names case-insensitively (e.g. Description vs description)
        lookup = {c.lower(): c for c in columns}
        missing = [c for c in EMAPPER_COLUMNS if c.lower() not in lookup]
        if missing:
            sys.exit(f"ERROR: {path} lacks column(s): {', '.join(missing)}")
        hits = {q: {c: r.get(lookup[c.lower()], "") for c in EMAPPER_COLUMNS}
                for q, r in hits.items()}
    return hits


def summarize(hit):
    ogs = clean(hit.get("eggNOG_OGs"))
    entries = [e for e in ogs.split(",") if e] if ogs else []
    first_og = entries[0].split("@", 1)[0] if entries else ""
    domain = ""
    for entry in entries:
        taxid = entry.split("@", 1)[-1].split("|", 1)[0]
        if taxid in DOMAINS:
            domain = DOMAINS[taxid]
            break
    max_lvl = clean(hit.get("max_annot_lvl"))
    return {
        "OG": first_og,
        "domain": domain,
        "max_annot_lvl": max_lvl.split("|", 1)[-1],
        "COG_category": clean(hit.get("COG_category")),
        "eggNOG_OGs": ogs,
        "Description": clean(hit.get("Description")),
        "Preferred_name": clean(hit.get("Preferred_name")),
        "PFAMs": clean(hit.get("PFAMs")),
    }


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("gff")
    p.add_argument("emapper_annotations")
    p.add_argument("-o", "--output", default="-")
    args = p.parse_args()

    hits = read_emapper(args.emapper_annotations)
    out = sys.stdout if args.output == "-" else open(args.output, "w", newline="")
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(OUT_COLUMNS)
    n_genes = n_hits = 0
    gene_ids = set()
    for gene_id, start, end, strand in read_gff(args.gff):
        n_genes += 1
        gene_ids.add(gene_id)
        hit = hits.get(gene_id)
        n_hits += hit is not None
        ann = summarize(hit) if hit else {}
        writer.writerow([gene_id, start, end, strand] +
                        [ann.get(c, "") for c in OUT_COLUMNS[4:]])
    if out is not sys.stdout:
        out.close()

    unmatched = set(hits) - gene_ids
    if unmatched:
        sys.exit(f"ERROR: {len(unmatched)} emapper queries not in the GFF, "
                 f"e.g. {sorted(unmatched)[:3]}")
    print(f"{args.gff}: {n_genes} genes, {n_hits} with eggNOG hits", file=sys.stderr)


if __name__ == "__main__":
    main()
