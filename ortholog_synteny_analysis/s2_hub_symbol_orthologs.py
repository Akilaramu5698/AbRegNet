"""
s4_hub_symbol_orthologs.py -- RBH of the hub PROTEIN SEQUENCES against all six
strains, and propagation of the hub gene symbol to every ortholog that clears
the cutoff.

1. Reads hub_ids.csv (gene, locus_tag, ncbi_protein) and hub_proteins.faa.
 
2. For every strain (the reference included, as a self-check):
     forward  hub protein  ->  that strain's full proteome
     reverse  the hit      ->  the reference full proteome
   Both directions use BLOSUM62, gap open 11, gap extend 1, and the same
   cutoff: E <= 1e-10, query coverage >= 85%, identity >= 85%.
   A hit is an ortholog only if it clears the cutoff in BOTH directions AND
   the reverse best hit returns to the originating reference locus (RBH).
3. Every confirmed ortholog inherits the hub's gene symbol. 

Outputs (results_ortho/)
  hub_symbol_orthologs_long.csv   one row per hub locus x strain, with
                                  identity, coverage, E, bitscore, rbh_confirmed
  hub_symbol_orthologs_matrix.csv hub locus x strain matrix of ortholog loci
  symbol_assignments.csv          strain locus -> assigned gene symbol
                                  (the propagated annotation)
  symbol_assignment_conflicts.csv target loci claimed by more than one hub locus
  hub_symbol_summary.csv          per gene symbol: strains covered, mean identity
  s4_report.md

Run:  python s4_hub_symbol_orthologs.py
"""
from __future__ import annotations

import os
import shutil
import subprocess

import pandas as pd

IN = "data_strains"
OUT = "results_ortho"
BLASTDIR = os.path.join(OUT, "blast_raw")
REFERENCE = "ATCC_19606"

HUB_IDS = "hub_ids.csv"
HUB_FASTA = "hub_proteins.faa"
CUTOFF = {"evalue": 1e-10, "qcovs": 85.0, "pident": 85.0}
MATRIX, GAPOPEN, GAPEXTEND = "BLOSUM62", 11, 1
THREADS = max(1, (os.cpu_count() or 4) - 1)
COLS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovs"]


def load_hub_ids(path: str = HUB_IDS) -> pd.DataFrame:
    sep = "\t" if path.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(path, sep=sep)
    c = {k.lower().strip(): k for k in df.columns}
    g = c.get("gene") or c.get("hub") or c.get("query") or df.columns[0]
    l = c.get("locus_tag") or c.get("locus") or df.columns[1]
    p = c.get("ncbi_protein") or c.get("protein_id") or c.get("protein")
    out = pd.DataFrame({"gene": df[g].astype(str).str.strip().str.lower(),
                        "locus_tag": df[l].astype(str).str.strip(),
                        "ncbi_protein": (df[p].astype(str).str.strip()
                                         if p else "")})
    return out.drop_duplicates().reset_index(drop=True)


def read_fasta(path: str) -> list:
    recs, h, seq = [], None, []
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            if h is not None:
                recs.append((h, "".join(seq)))
            h, seq = line[1:].strip(), []
        elif h is not None:
            seq.append(line.strip())
    if h is not None:
        recs.append((h, "".join(seq)))
    return recs


def build_query(hub_ids: pd.DataFrame, fasta: str, qpath: str):
    """Write a query FASTA headed '<gene>|<locus_tag>'; report unmatched."""
    by_locus = {r.locus_tag: r for r in hub_ids.itertuples()}
    by_prot = {r.ncbi_protein: r for r in hub_ids.itertuples() if r.ncbi_protein}
    by_prot_nover = {k.split(".")[0]: v for k, v in by_prot.items()}
    used, unmatched, written = set(), [], 0
    with open(qpath, "w", encoding="utf-8") as fh:
        for header, seq in read_fasta(fasta):
            hit = None
            for tok in [header] + header.replace("|", " ").replace(",", " ").split():
                tok = tok.strip()
                hit = (by_locus.get(tok) or by_prot.get(tok)
                       or by_prot_nover.get(tok.split(".")[0]))
                if hit is not None:
                    break
            if hit is None or not seq:
                unmatched.append(header)
                continue
            fh.write(f">{hit.gene}|{hit.locus_tag}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")
            used.add(hit.locus_tag)
            written += 1
    missing = [r.locus_tag for r in hub_ids.itertuples() if r.locus_tag not in used]
    return written, unmatched, missing


def makedb(strain: str) -> str:
    db = os.path.join(BLASTDIR, strain)
    if not (os.path.exists(db + ".pin") or os.path.exists(db + ".pdb")):
        subprocess.run(["makeblastdb", "-in", os.path.join(IN, f"{strain}.faa"),
                        "-dbtype", "prot", "-out", db],
                       check=True, capture_output=True, text=True)
    return db


def blastp(query: str, db: str, out: str) -> pd.DataFrame:
    subprocess.run(["blastp", "-query", query, "-db", db,
                    "-outfmt", "6 " + " ".join(COLS),
                    "-matrix", MATRIX, "-gapopen", str(GAPOPEN),
                    "-gapextend", str(GAPEXTEND),
                    "-evalue", str(CUTOFF["evalue"]), "-max_target_seqs", "5",
                    "-num_threads", str(THREADS), "-out", out],
                   check=True, capture_output=True, text=True)
    return pd.read_csv(out, sep="\t", names=COLS)


def passes(df: pd.DataFrame) -> pd.DataFrame:
    return df[(df.evalue <= CUTOFF["evalue"]) & (df.qcovs >= CUTOFF["qcovs"]) &
              (df.pident >= CUTOFF["pident"])]


def best(df: pd.DataFrame) -> pd.DataFrame:
    df = df.sort_values(["qseqid", "bitscore", "evalue"],
                        ascending=[True, False, True])
    return df.groupby("qseqid", as_index=False).first()


def write_subset(strain: str, ids, path: str):
    keep, out, cur = set(ids), [], None
    for header, seq in read_fasta(os.path.join(IN, f"{strain}.faa")):
        tok = header.split()[0]
        if tok in keep:
            out.append((tok, seq))
    with open(path, "w", encoding="utf-8") as fh:
        for t, s in out:
            fh.write(f">{t}\n")
            for i in range(0, len(s), 60):
                fh.write(s[i:i + 60] + "\n")
    return len(out)


def main():
    if shutil.which("blastp") is None:
        raise SystemExit("BLAST+ not found on PATH")
    os.makedirs(BLASTDIR, exist_ok=True)
    hub_ids = load_hub_ids()
    qpath = os.path.join(BLASTDIR, "hub_query.faa")
    n, unmatched, missing = build_query(hub_ids, HUB_FASTA, qpath)
    print(f"  {HUB_IDS}: {len(hub_ids)} hub loci over {hub_ids.gene.nunique()} "
          f"symbols; query built from {n} sequences")
    if unmatched:
        print(f"  FASTA records not matched to any hub id ({len(unmatched)}): "
              + ", ".join(unmatched[:5]) + (" ..." if len(unmatched) > 5 else ""))
    if missing:
        print(f"  hub loci with no sequence in {HUB_FASTA} ({len(missing)}): "
              + ", ".join(missing[:5]) + (" ..." if len(missing) > 5 else ""))

    strains = sorted(f[:-4] for f in os.listdir(IN) if f.endswith(".faa"))
    rows = []
    for s in strains:
        fwd_all = blastp(qpath, makedb(s),
                         os.path.join(BLASTDIR, f"hubq__vs__{s}.tsv"))
        fwd = best(passes(fwd_all))
        rpath = os.path.join(BLASTDIR, f"hubq_rev_{s}.faa")
        write_subset(s, set(fwd["sseqid"]), rpath)
        rev = best(passes(blastp(rpath, makedb(REFERENCE),
                                 os.path.join(BLASTDIR, f"{s}__vs__{REFERENCE}_hubrev.tsv"))))
        rev_i = rev.set_index("qseqid")
        for _, r in fwd.iterrows():
            gene, ref_locus = r["qseqid"].split("|", 1)
            back = rev_i.at[r["sseqid"], "sseqid"] if r["sseqid"] in rev_i.index else None
            rows.append({
                "gene": gene, "reference_locus": ref_locus, "strain": s,
                "ortholog_locus": r["sseqid"],
                "pident": round(float(r["pident"]), 2),
                "qcovs": float(r["qcovs"]), "evalue": float(r["evalue"]),
                "bitscore": float(r["bitscore"]),
                "reverse_best_hit": back,
                "rbh_confirmed": back == ref_locus,
                "passes_cutoff": True,
                "assigned_symbol": gene if back == ref_locus else ""})
        print(f"  {s}: {len(fwd)} hub queries with a passing hit, "
              f"{sum(1 for x in rows if x['strain'] == s and x['rbh_confirmed'])} RBH-confirmed")

    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(OUT, "hub_symbol_orthologs_long.csv"), index=False)

    ok = df[df["rbh_confirmed"]]
    mat = ok.pivot_table(index=["gene", "reference_locus"], columns="strain",
                         values="ortholog_locus", aggfunc="first").reset_index()
    mat.to_csv(os.path.join(OUT, "hub_symbol_orthologs_matrix.csv"), index=False)

    assign = (ok[["strain", "ortholog_locus", "assigned_symbol",
                  "reference_locus", "pident", "qcovs", "evalue"]]
              .rename(columns={"ortholog_locus": "locus_tag"})
              .sort_values(["strain", "assigned_symbol", "locus_tag"]))
    assign.to_csv(os.path.join(OUT, "symbol_assignments.csv"), index=False)

    conf = (assign.groupby(["strain", "locus_tag"])
            .agg(n_hub_loci_claiming=("reference_locus", "nunique"),
                 symbols=("assigned_symbol", lambda s: ";".join(sorted(set(s)))),
                 claiming_reference_loci=("reference_locus",
                                          lambda s: ";".join(sorted(set(s)))))
            .reset_index())
    conf = conf[conf["n_hub_loci_claiming"] > 1]
    conf.to_csv(os.path.join(OUT, "symbol_assignment_conflicts.csv"), index=False)

    summ = (ok.groupby("gene")
            .agg(n_reference_loci=("reference_locus", "nunique"),
                 n_strains_covered=("strain", "nunique"),
                 n_orthologs=("ortholog_locus", "nunique"),
                 mean_pident=("pident", "mean"), min_pident=("pident", "min"))
            .round(2).reset_index())
    summ.to_csv(os.path.join(OUT, "hub_symbol_summary.csv"), index=False)

    L = ["# Hub symbol propagation by reciprocal best hit", "",
         f"Cutoff: E <= {CUTOFF['evalue']}, query coverage >= {CUTOFF['qcovs']}%, "
         f"identity >= {CUTOFF['pident']}% in both directions, plus the reverse "
         "best hit returning to the originating reference locus.", "",
         "| gene | reference loci | strains covered | orthologs assigned | "
         "mean identity | min identity |", "|---|---|---|---|---|---|"]
    for _, r in summ.iterrows():
        L.append(f"| {r['gene']} | {r['n_reference_loci']} | "
                 f"{r['n_strains_covered']} | {r['n_orthologs']} | "
                 f"{r['mean_pident']} | {r['min_pident']} |")
    if len(conf):
        L += ["", f"## Conflicting assignments ({len(conf)})", "",
              "These target loci are the RBH of more than one hub locus, so the "
              "symbol they inherit is ambiguous:", "",
              "| strain | locus | symbols | claiming reference loci |",
              "|---|---|---|---|"]
        for _, r in conf.iterrows():
            L.append(f"| {r['strain']} | {r['locus_tag']} | {r['symbols']} | "
                     f"{r['claiming_reference_loci']} |")
    with open(os.path.join(OUT, "s4_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print(f"  wrote {OUT}/symbol_assignments.csv ({len(assign)} assignments), "
          f"{len(conf)} conflicts")


if __name__ == "__main__":
    main()
