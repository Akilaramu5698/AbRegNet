"""
s3_microsynteny.py -- conserved chromosomal microsynteny around the hub loci.


Outputs (results_ortho/)
  hub_locus_mapping.csv       hub -> reference locus, and how it was mapped
  hub_locus_candidates.csv    ambiguous product-keyword candidates to review
  hub_ortholog_table.csv      hub -> ortholog locus in each strain + identity
  microsynteny_scores.csv     hub x strain: S_syn and its components
  microsynteny_matrix.csv     hub x strain matrix of S_syn
  microsynteny_windows.csv    the actual gene windows, per hub per strain
  s3_report.md

Run:  python s3_microsynteny.py
"""
from __future__ import annotations

import os

import numpy as np
import pandas as pd

IN = "data_strains"
OUT = "results_ortho"
REFERENCE = "ATCC_19606"
WINDOW = 2                      # +/- this many CDS
OVERRIDES = "hub_locus_overrides.csv"
SCORE_ALL_PARALOGS = True   # score EVERY paralog listed in hub_ids.csv
HUB_THRESH = {"evalue": 1e-10, "qcovs": 85.0, "pident": 85.0}
MATRIX, GAPOPEN, GAPEXTEND = "BLOSUM62", 11, 1

HUBS = ["fadr", "rcsc", "beti", "dmlr", "pgrr", "yybr",
        "arsr2", "paax", "lrp", "rhas", "rpoh"]



# --------------------------------------------------------------- hub input
HUB_IDS = "hub_ids.csv"          # columns: gene, locus_tag, ncbi_protein
HUB_FASTA = "hub_proteins.faa"   # one record per hub locus (see below)


def load_hub_ids(path: str = HUB_IDS) -> pd.DataFrame:
    """Read the hub id table: one row per (hub symbol, reference locus).

    Accepts gene|hub|query for the symbol column, locus_tag|locus for the
    locus, and ncbi_protein|protein_id for the accession. Multi-copy symbols
    are kept as separate rows - 9 of the 11 hub symbols are multi-copy in this
    genome, so a symbol does not identify a locus and every paralog is carried
    through to the output.
    """
    if not os.path.exists(path):
        return pd.DataFrame(columns=["gene", "locus_tag", "ncbi_protein"])
    sep = "\t" if path.endswith((".tsv", ".tab")) else ","
    df = pd.read_csv(path, sep=sep)
    c = {k.lower().strip(): k for k in df.columns}
    gcol = c.get("gene") or c.get("hub") or c.get("query") or df.columns[0]
    lcol = c.get("locus_tag") or c.get("locus") or df.columns[1]
    pcol = c.get("ncbi_protein") or c.get("protein_id") or c.get("protein")
    out = pd.DataFrame({
        "gene": df[gcol].astype(str).str.strip().str.lower(),
        "locus_tag": df[lcol].astype(str).str.strip(),
        "ncbi_protein": (df[pcol].astype(str).str.strip() if pcol else "")})
    return out.drop_duplicates().reset_index(drop=True)


def read_fasta(path: str) -> dict:
    """FASTA -> {header_token: sequence} for every whitespace/pipe token."""
    recs, cur = {}, None
    for line in open(path, encoding="utf-8"):
        if line.startswith(">"):
            cur = line[1:].strip()
            recs[cur] = []
        elif cur is not None:
            recs[cur].append(line.strip())
    out = {}
    for h, parts in recs.items():
        seq = "".join(parts)
        for tok in [h] + h.replace("|", " ").split():
            out.setdefault(tok.strip(), seq)
    return out


def verify_hub_ids(hub_ids: pd.DataFrame, ref_cds: pd.DataFrame,
                   fasta_path: str = HUB_FASTA) -> pd.DataFrame:
    """Check each locus exists in the reference annotation, that the protein
    accession agrees, and that a sequence is available for it."""
    by_locus = ref_cds.set_index("locus_tag")
    fa = read_fasta(fasta_path) if os.path.exists(fasta_path) else {}
    rows = []
    for _, r in hub_ids.iterrows():
        present = r["locus_tag"] in by_locus.index
        acc = by_locus.at[r["locus_tag"], "protein_id"] if present else ""
        rows.append({
            "gene": r["gene"], "locus_tag": r["locus_tag"],
            "ncbi_protein": r["ncbi_protein"],
            "locus_in_reference_annotation": present,
            "annotation_protein_id": acc,
            "protein_id_matches": bool(present and str(acc).split(".")[0] ==
                                       str(r["ncbi_protein"]).split(".")[0]),
            "sequence_available": bool(fa.get(r["locus_tag"])
                                       or fa.get(r["ncbi_protein"])),
        })
    return pd.DataFrame(rows)


def blast_hub_fasta(ref_strain: str) -> pd.DataFrame:
    """BLASTP the hub protein FASTA against the reference proteome."""
    import shutil
    import subprocess
    if not os.path.exists(HUB_FASTA):
        return pd.DataFrame()
    if shutil.which("blastp") is None:
        print(f"  {HUB_FASTA} found but BLAST+ is not on PATH - skipping")
        return pd.DataFrame()
    db = os.path.join(OUT, "blast_raw", ref_strain)
    if not (os.path.exists(db + ".pin") or os.path.exists(db + ".pdb")):
        os.makedirs(os.path.join(OUT, "blast_raw"), exist_ok=True)
        subprocess.run(["makeblastdb", "-in", os.path.join(IN, f"{ref_strain}.faa"),
                        "-dbtype", "prot", "-out", db], check=True,
                       capture_output=True, text=True)
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovs"]
    out = os.path.join(OUT, "blast_raw", "hubs__vs__" + ref_strain + ".tsv")
    subprocess.run(["blastp", "-query", HUB_FASTA, "-db", db,
                    "-outfmt", "6 " + " ".join(cols),
                    "-matrix", MATRIX, "-gapopen", str(GAPOPEN),
                    "-gapextend", str(GAPEXTEND),
                    "-evalue", str(HUB_THRESH["evalue"]),
                    "-max_target_seqs", "5",
                    "-num_threads", str(max(1, (os.cpu_count() or 4) - 1)),
                    "-out", out], check=True, capture_output=True, text=True)
    df = pd.read_csv(out, sep="\t", names=cols)
    df = df[(df.evalue <= HUB_THRESH["evalue"]) &
            (df.qcovs >= HUB_THRESH["qcovs"]) &
            (df.pident >= HUB_THRESH["pident"])]
    df = df.sort_values(["qseqid", "bitscore", "evalue"],
                        ascending=[True, False, True])
    return df.groupby("qseqid", as_index=False).first()


def load_cds() -> dict:
    out = {}
    for f in sorted(os.listdir(IN)):
        if f.endswith("_cds.csv"):
            s = f[:-8]
            df = pd.read_csv(f, sep=",") if False else pd.read_csv(os.path.join(IN, f))
            out[s] = df
    return out


def ortholog_maps(rbh: pd.DataFrame, strains: list) -> dict:
    """{(strainA, strainB): {locusA: locusB}} for all ordered pairs present."""
    maps = {}
    for _, r in rbh.iterrows():
        maps.setdefault((r.strain_A, r.strain_B), {})[r.locus_A] = r.locus_B
        maps.setdefault((r.strain_B, r.strain_A), {})[r.locus_B] = r.locus_A
    return maps


def map_hubs(ref_cds: pd.DataFrame, ref_strain: str = REFERENCE):
    """Return (mapping DataFrame, candidates DataFrame)."""
    gene = ref_cds["gene"].fillna("").str.lower()
    prod = ref_cds["product"].fillna("").str.lower()
    rows, cands = [], []
    # Hub loci come straight from hub_ids.csv (gene, locus_tag, ncbi_protein);
    # no symbol or sequence matching is needed for the reference strain. Every
    # paralog listed is scored separately.
    cand_rows = []
    hub_ids = load_hub_ids()
    if SCORE_ALL_PARALOGS and len(hub_ids):
        ver = verify_hub_ids(hub_ids, ref_cds)
        ver.to_csv(os.path.join(OUT, "hub_id_verification.csv"), index=False)
        bad = ver[~ver["locus_in_reference_annotation"]]
        mism = ver[ver["locus_in_reference_annotation"] & ~ver["protein_id_matches"]]
        print(f"  {HUB_IDS}: {len(ver)} hub loci over {ver['gene'].nunique()} "
              f"symbols; {len(bad)} not in the reference annotation, "
              f"{len(mism)} with a protein-id mismatch")
        counts = hub_ids.groupby("gene")["locus_tag"].nunique().to_dict()
        for _, r in hub_ids.iterrows():
            if r["locus_tag"] in set(ref_cds["locus_tag"]):
                cand_rows.append({"hub": r["gene"], "reference_locus": r["locus_tag"],
                                  "mapped_by": "hub_ids_csv",
                                  "n_candidates": counts[r["gene"]],
                                  "ncbi_protein": r["ncbi_protein"]})
    blast_map = {}
    bh = blast_hub_fasta(ref_strain)
    if len(bh):
        blast_map = {str(q).strip().lower(): r for q, r in
                     zip(bh["qseqid"], bh.to_dict("records"))}
        print(f"  {HUB_FASTA}: {len(blast_map)} hubs mapped by BLASTP")
    overrides = {}
    if os.path.exists(OVERRIDES):
        ov = pd.read_csv(OVERRIDES)
        overrides = {str(h).strip().lower(): str(l).strip()
                     for h, l in zip(ov["hub"], ov["locus_tag"])}
    by_hub_cands = {}
    for r in cand_rows:
        by_hub_cands.setdefault(r["hub"], []).append(r)
    for h in HUBS:
        if h not in overrides and h not in blast_map and h in by_hub_cands:
            rows += by_hub_cands[h]
            continue
        if h in overrides:
            rows.append({"hub": h, "reference_locus": overrides[h],
                         "mapped_by": "user_override", "n_candidates": 1})
            continue
        if h in blast_map:
            b = blast_map[h]
            rows.append({"hub": h, "reference_locus": b["sseqid"],
                         "mapped_by": "blastp_hub_fasta", "n_candidates": 1,
                         "pident": round(float(b["pident"]), 2),
                         "qcovs": float(b["qcovs"]), "evalue": float(b["evalue"]),
                         "bitscore": float(b["bitscore"])})
            continue
        exact = ref_cds[gene == h]
        if len(exact) == 1:
            rows.append({"hub": h, "reference_locus": exact.iloc[0]["locus_tag"],
                         "mapped_by": "exact_gene_symbol", "n_candidates": 1})
            continue
        hit = ref_cds[prod.str.contains(h, regex=False)]
        for _, c in hit.iterrows():
            cands.append({"hub": h, "candidate_locus": c["locus_tag"],
                          "gene": c["gene"], "product": c["product"]})
        rows.append({"hub": h, "reference_locus": np.nan,
                     "mapped_by": ("ambiguous_product_match" if len(hit)
                                   else "unmapped"),
                     "n_candidates": len(hit)})
    return pd.DataFrame(rows), pd.DataFrame(cands)


def window_of(cds: pd.DataFrame, locus: str, w: int = None):
    """Rows of the +/-w CDS window around `locus`, plus the hub row.

    w defaults to the module-level WINDOW read AT CALL TIME, so a caller that
    rebinds WINDOW (the window-sensitivity check in s5) actually changes the
    window - a default argument would have frozen it at import.
    """
    w = WINDOW if w is None else w
    hit = cds.index[cds["locus_tag"] == locus]
    if len(hit) == 0:
        return None, None
    i = hit[0]
    hub = cds.loc[i]
    same = cds[cds["replicon"] == hub["replicon"]]
    idx = hub["cds_index"]
    win = same[(same["cds_index"] >= idx - w) & (same["cds_index"] <= idx + w)].copy()
    win["offset"] = win["cds_index"] - idx
    return hub, win


def score_hub(hub_name, ref_locus, ref_cds, tgt_strain, tgt_cds, omap):
    """S_syn and components for one hub in one target strain."""
    ref_hub, ref_win = window_of(ref_cds, ref_locus)
    if ref_hub is None:
        return None, []
    tgt_locus = omap.get(ref_locus)
    base = {"hub": hub_name, "reference_locus": ref_locus,
            "strain": tgt_strain, "ortholog_locus": tgt_locus,
            "hub_has_ortholog": tgt_locus is not None}
    if tgt_locus is None:
        base.update({"S_syn": 0.0, "S_orthology": 0.0, "S_in_window": 0.0,
                     "S_order": 0.0, "S_strand": 0.0, "inverted": None,
                     "n_flankers": len(ref_win) - 1})
        return base, []
    tgt_hub, tgt_win = window_of(tgt_cds, tgt_locus)
    inverted = (ref_hub["strand"] != tgt_hub["strand"])
    sign = -1 if inverted else 1
    tgt_by_locus = {r["locus_tag"]: r for _, r in tgt_win.iterrows()}
    tgt_all = set(tgt_cds["locus_tag"])

    n = ortho = inwin = order = strand_ok = 0
    detail = []
    for _, f in ref_win.iterrows():
        if f["offset"] == 0:
            continue
        n += 1
        o = omap.get(f["locus_tag"])
        has_o = o is not None and o in tgt_all
        in_w = has_o and o in tgt_by_locus
        exp_off = sign * f["offset"]
        at_off = in_w and tgt_by_locus[o]["offset"] == exp_off
        exp_strand = (f["strand"] if not inverted
                      else ("+" if f["strand"] == "-" else "-"))
        s_ok = in_w and tgt_by_locus[o]["strand"] == exp_strand
        ortho += has_o; inwin += in_w; order += at_off; strand_ok += s_ok
        detail.append({"hub": hub_name, "strain": tgt_strain,
                       "ref_locus": f["locus_tag"], "ref_offset": f["offset"],
                       "ref_strand": f["strand"], "ref_gene": f["gene"],
                       "ortholog_locus": o, "has_ortholog": has_o,
                       "in_window": in_w,
                       "target_offset": tgt_by_locus[o]["offset"] if in_w else None,
                       "expected_offset": exp_off,
                       "order_conserved": at_off,
                       "target_strand": tgt_by_locus[o]["strand"] if in_w else None,
                       "expected_strand": exp_strand,
                       "strand_conserved": s_ok,
                       "counts_towards_S_syn": bool(at_off and s_ok),
                       "neighbourhood_inverted": inverted})
    conserved = sum(1 for d in detail if d["counts_towards_S_syn"])
    base.update({"n_flankers": n, "inverted": inverted,
                 "S_syn": round(conserved / n, 4) if n else np.nan,
                 "S_orthology": round(ortho / n, 4) if n else np.nan,
                 "S_in_window": round(inwin / n, 4) if n else np.nan,
                 "S_order": round(order / n, 4) if n else np.nan,
                 "S_strand": round(strand_ok / n, 4) if n else np.nan})
    return base, detail


def main():
    cds = load_cds()
    if REFERENCE not in cds:
        raise SystemExit(f"reference {REFERENCE} not found in {IN}")
    rbh = pd.read_csv(os.path.join(OUT, "rbh_pairs.csv"))
    maps = ortholog_maps(rbh, list(cds))
    ref_cds = cds[REFERENCE]

    mapping, candidates = map_hubs(ref_cds)
    mapping.to_csv(os.path.join(OUT, "hub_locus_mapping.csv"), index=False)
    candidates.to_csv(os.path.join(OUT, "hub_locus_candidates.csv"), index=False)
    mapped = mapping.dropna(subset=["reference_locus"])
    print(f"  hubs mapped to a reference locus: {len(mapped)}/{len(mapping)}")
    unmapped = mapping[mapping["reference_locus"].isna()]
    if len(unmapped):
        print("  NOT mapped (supply hub_locus_overrides.csv): "
              + ", ".join(f"{r.hub}({r.mapped_by},{r.n_candidates} candidates)"
                          for _, r in unmapped.iterrows()))

    scores, details, ortho_rows = [], [], []
    others = [s for s in cds if s != REFERENCE]
    for _, m in mapped.iterrows():
        row = {"hub": m["hub"], f"{REFERENCE}_locus": m["reference_locus"]}
        for s in others:
            omap = maps.get((REFERENCE, s), {})
            sc, det = score_hub(m["hub"], m["reference_locus"], ref_cds, s,
                                cds[s], omap)
            if sc:
                scores.append(sc)
                details += det
                row[f"{s}_locus"] = sc["ortholog_locus"]
        ortho_rows.append(row)

    sc = pd.DataFrame(scores)
    sc.to_csv(os.path.join(OUT, "microsynteny_scores.csv"), index=False)
    pd.DataFrame(details).to_csv(os.path.join(OUT, "microsynteny_windows.csv"),
                                 index=False)
    pd.DataFrame(ortho_rows).to_csv(os.path.join(OUT, "hub_ortholog_table.csv"),
                                    index=False)
    if len(sc):
        mat = sc.pivot_table(index="hub", columns="strain", values="S_syn")
        mat["mean_S_syn"] = mat.mean(axis=1).round(4)
        mat.reset_index().to_csv(os.path.join(OUT, "microsynteny_matrix.csv"),
                                 index=False)

        L = ["# Microsynteny around the hub loci", "",
             f"Reference: {REFERENCE}. Window: +/-{WINDOW} CDS "
             f"({2*WINDOW} flanking genes). S_syn = fraction of flankers "
             "conserving orthology, relative order and strand.", "",
             "| hub | strain | ortholog | inverted | S_syn | orthology | "
             "in window | order | strand |", "|---|---|---|---|---|---|---|---|---|"]
        for _, r in sc.iterrows():
            L.append(f"| {r['hub']} | {r['strain']} | {r['ortholog_locus']} | "
                     f"{r['inverted']} | {r['S_syn']} | {r['S_orthology']} | "
                     f"{r['S_in_window']} | {r['S_order']} | {r['S_strand']} |")
        if len(unmapped):
            L += ["", "## Hubs without a reference locus", "",
                  "These were not scored; map them in hub_locus_overrides.csv "
                  "(columns hub,locus_tag) and re-run:", ""]
            for _, r in unmapped.iterrows():
                L.append(f"- {r['hub']} ({r['mapped_by']}, "
                         f"{r['n_candidates']} product-keyword candidates)")
        with open(os.path.join(OUT, "s3_report.md"), "w", encoding="utf-8") as fh:
            fh.write("\n".join(L))
        print(f"  wrote {OUT}/microsynteny_scores.csv ({len(sc)} hub x strain rows)")


if __name__ == "__main__":
    main()
