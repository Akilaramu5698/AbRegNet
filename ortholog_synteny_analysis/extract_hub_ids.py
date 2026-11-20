#!/usr/bin/env python3
"""
Extract locus_tag and NCBI protein accession for a list of gene symbols
from a genome annotation table.

Usage
-----
    python extract_hub_ids.py annotation.csv hub_genes.txt -o hub_ids.csv



import argparse
import csv
import sys
from collections import Counter

# Column name synonyms -> canonical name. Extend if your file differs.
COLUMN_ALIASES = {
    "symbol":       ["gene_name", "gene", "symbol", "gene_symbol", "name"],
    "locus_tag":    ["locus_tag", "locustag", "locus", "old_locus_tag"],
    "protein":      ["ncbi_protein", "protein_id", "protein_accession",
                     "ncbi_protein_accession", "refseq_protein"],
}
# Extra columns copied through to the output when present.
PASSTHROUGH = ["start", "end", "strand", "gene_product", "product",
               "COG", "COG_Family_Name"]


def sniff_delimiter(path, encoding="utf-8"):
    """Guess the field separator from the header line."""
    with open(path, newline="", encoding=encoding) as fh:
        head = fh.readline()
    counts = {d: head.count(d) for d in ["\t", ",", ";"]}
    delim, n = max(counts.items(), key=lambda kv: kv[1])
    if n == 0:
        raise SystemExit(f"Could not detect a delimiter in {path!r} "
                         f"(header has no tab, comma or semicolon).")
    return delim


def resolve_columns(fieldnames):
    """Map canonical names onto this file's actual header names."""
    lower = {f.strip().lower(): f for f in fieldnames}
    resolved = {}
    for canon, aliases in COLUMN_ALIASES.items():
        for a in aliases:
            if a.lower() in lower:
                resolved[canon] = lower[a.lower()]
                break
        else:
            raise SystemExit(
                f"Required column for {canon!r} not found. Looked for "
                f"{aliases}; file has {list(fieldnames)}")
    return resolved


def read_symbols(path):
    """One symbol per line; blank lines and '#' comments ignored."""
    out = []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                out.append(s)
    seen, uniq = set(), []
    for s in out:                       # keep input order, drop repeats
        if s.lower() not in seen:
            seen.add(s.lower())
            uniq.append(s)
    return uniq


def extract(annotation, symbols_file, out_csv, encoding="utf-8"):
    symbols = read_symbols(symbols_file)
    delim = sniff_delimiter(annotation, encoding)

    with open(annotation, newline="", encoding=encoding) as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        cols = resolve_columns(reader.fieldnames)
        extras = [c for c in PASSTHROUGH if c in reader.fieldnames]
        rows = [{k: (v or "").strip() for k, v in r.items() if k is not None}
                for r in reader]

    want = {s.lower(): s for s in symbols}
    hits = {s: [] for s in symbols}
    for r in rows:
        key = r[cols["symbol"]].lower()
        if key in want:
            hits[want[key]].append(r)

    out_fields = (["query_symbol", "matched_symbol", "locus_tag",
                   "ncbi_protein"] + extras + ["n_copies", "status"])
    with open(out_csv, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=out_fields)
        w.writeheader()
        for s in symbols:
            rs = hits[s]
            if not rs:
                w.writerow({"query_symbol": s, "matched_symbol": "",
                            "locus_tag": "", "ncbi_protein": "",
                            "n_copies": 0, "status": "not_found"})
                continue
            for r in rs:
                rec = {"query_symbol": s,
                       "matched_symbol": r[cols["symbol"]],
                       "locus_tag": r[cols["locus_tag"]],
                       "ncbi_protein": r[cols["protein"]],
                       "n_copies": len(rs),
                       "status": "ok" if len(rs) == 1 else "multi_copy"}
                rec.update({c: r[c] for c in extras})
                w.writerow(rec)

    # ---- report -------------------------------------------------------
    n_found = sum(1 for s in symbols if hits[s])
    total = sum(len(v) for v in hits.values())
    print(f"annotation : {annotation}  (delimiter {delim!r}, {len(rows)} rows)")
    print(f"symbols    : {len(symbols)} queried, {n_found} found, "
          f"{len(symbols) - n_found} not found")
    print(f"output     : {out_csv}  ({total} rows)\n")
    for s in symbols:
        rs = hits[s]
        if not rs:
            print(f"  {s:<8} NOT FOUND")
        else:
            ids = "; ".join(f"{r[cols['locus_tag']]}/{r[cols['protein']]}"
                            for r in rs)
            flag = "" if len(rs) == 1 else "  <-- multi-copy"
            print(f"  {rs[0][cols['symbol']]:<8} n={len(rs):<3}{ids}{flag}")

    missing = [s for s in symbols if not hits[s]]
    multi = [s for s in symbols if len(hits[s]) > 1]
    blank = total - sum(1 for s in symbols for r in hits[s]
                        if r[cols["protein"]])
    if multi:
        print(f"\nWARNING: {len(multi)} symbol(s) map to more than one locus "
              f"({', '.join(multi)}).\n"
              "         The symbol alone does not identify a single protein; "
              "disambiguate by\n         locus tag before any protein-level "
              "work.")
    if blank:
        print(f"\nWARNING: {blank} matched row(s) have an empty protein "
              "accession (non-coding feature?).")
    if missing:
        print(f"\nWARNING: not found: {', '.join(missing)}")
    return total


def main():
    ap = argparse.ArgumentParser(
        description="Pull locus_tag + NCBI protein accession for a list of "
                    "gene symbols out of an annotation table.")
    ap.add_argument("annotation", help="annotation table (.csv or .tsv)")
    ap.add_argument("symbols", help="text file, one gene symbol per line")
    ap.add_argument("-o", "--out", default="hub_ids.csv",
                    help="output CSV (default: hub_ids.csv)")
    ap.add_argument("--encoding", default="utf-8")
    a = ap.parse_args()
    extract(a.annotation, a.symbols, a.out, a.encoding)


if __name__ == "__main__":
    main()
