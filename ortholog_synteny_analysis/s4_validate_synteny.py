"""
s5_validate_synteny.py -- validation of the microsynteny analysis.

Four checks, in order of how much they matter.

1. GENOME BACKGROUND (the essential one). S_syn is recomputed for a random
   sample of non-hub reference loci. If the whole genome scores ~0.9, a hub
   scoring 0.9 says nothing. Each hub locus gets an empirical two-sided
   position against that background: p_high = fraction of background loci
   with S_syn >= the hub's, p_low = fraction with S_syn <= it.

2. POSITIVE AND NEGATIVE CONTROLS. Ribosomal protein genes sit in long,
   strongly conserved operons and should score near 1. Transposases,
   integrases and phage genes sit in mobile elements and should score low.
   If these two classes are not separated, the metric is not measuring
   synteny and nothing else in the table can be trusted.

3. WINDOW-SIZE SENSITIVITY. S_syn is recomputed at +/-1, 2, 3 and 5 CDS.
   A conclusion that only holds at one window size is a window artefact.

4. SELF-CONSISTENCY. The reference scored against itself must give S_syn = 1
   for every locus; anything else is a bug in the scoring or the orthology.



Outputs (results_ortho/validation/)
  background_scores.csv        S_syn for every sampled background locus
  hub_vs_background.csv        per hub locus: S_syn, background mean, p_high, p_low
  control_scores.csv           ribosomal / mobile-element control loci
  window_sensitivity.csv       hub S_syn at +/-1, 2, 3, 5
  self_consistency.csv         reference-vs-reference check
  validation_report.md

Run:  python s5_validate_synteny.py
"""
from __future__ import annotations

import os
import random
import shutil
import subprocess

import sys

import numpy as np
import pandas as pd

# make the sibling scripts importable even when Python's safe-path mode keeps
# the script directory off sys.path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import s3_microsynteny as s3

IN = "data_strains"
OUT = os.path.join("results_ortho", "validation")
BLASTDIR = os.path.join("results_ortho", "blast_raw")
REFERENCE = "ATCC_19606"

N_BACKGROUND = 1000
WINDOWS = [1, 2, 3, 5]
SEED = 0
CUTOFF = {"evalue": 1e-10, "qcovs": 85.0, "pident": 85.0}
MATRIX, GAPOPEN, GAPEXTEND = "BLOSUM62", 11, 1
THREADS = max(1, (os.cpu_count() or 4) - 1)
COLS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "qcovs"]

POSITIVE_KEYWORDS = [
    # Ribosomal
    "ribosomal protein",
    "rRNA",
    "ribosomal RNA",
    "50S ribosomal",
    "30S ribosomal",

    # Core transcription/translation machinery
    "RNA polymerase subunit",
    "translation elongation factor",
    "translation initiation factor",
    "tRNA synthetase",
    "aminoacyl-tRNA",

    # DNA replication core
    "DNA gyrase",
    "DNA polymerase III",
    "replication protein DnaA",
    "replication initiator",

    # Essential cell division
    "cell division protein FtsZ",
    "cell division protein FtsA",

    # Core metabolism (housekeeping)
    "ATP synthase",
    "chaperone protein GroEL",
    "chaperone protein DnaK",
    "chaperone protein GroES",
]

NEGATIVE_KEYWORDS = [
    # Transposable elements
    "transposase",
    "insertion sequence",
    "IS element",
    "IS protein",
    "transposon",
    "Tn element",

    # Integrases and recombinases
    "integrase",
    "site-specific recombinase",
    "recombinase",
    "tyrosine recombinase",
    "serine recombinase",
    "resolvase",

    # Phage related
    "phage",
    "bacteriophage",
    "prophage",
    "phage integrase",
    "phage tail",
    "phage capsid",
    "phage coat",
    "phage structural",
    "terminase",
    "endolysin",

    # Plasmid related
    "plasmid replication",
    "plasmid mobilization",
    "mobilization protein",
    "conjugal transfer",
    "conjugative",
    "relaxase",

    # Resistance islands (mobile)
    "resistance island",
    "genomic island",
    "pathogenicity island",

    # Toxin-antitoxin (often mobile)
    "toxin-antitoxin",
    "addiction module",

    # Hypothetical/unknown (often mobile)
    "hypothetical protein",
    "unknown function",
    "uncharacterized protein",
]

def makedb(strain):
    db = os.path.join(BLASTDIR, strain)
    if not (os.path.exists(db + ".pin") or os.path.exists(db + ".pdb")):
        subprocess.run(["makeblastdb", "-in", os.path.join(IN, f"{strain}.faa"),
                        "-dbtype", "prot", "-out", db], check=True,
                       capture_output=True, text=True)
    return db


def blast_best(query_faa, db, out):
    if not os.path.exists(out):
        subprocess.run(["blastp", "-query", query_faa, "-db", db,
                        "-outfmt", "6 " + " ".join(COLS), "-matrix", MATRIX,
                        "-gapopen", str(GAPOPEN), "-gapextend", str(GAPEXTEND),
                        "-evalue", str(CUTOFF["evalue"]), "-max_target_seqs", "5",
                        "-num_threads", str(THREADS), "-out", out],
                       check=True, capture_output=True, text=True)
    df = pd.read_csv(out, sep="\t", names=COLS)
    df = df[(df.evalue <= CUTOFF["evalue"]) & (df.qcovs >= CUTOFF["qcovs"]) &
            (df.pident >= CUTOFF["pident"])]
    df = df.sort_values(["qseqid", "bitscore", "evalue"],
                        ascending=[True, False, True])
    return df.groupby("qseqid", as_index=False).first()


def genome_wide_rbh(strains) -> dict:
    """{strain: {ref_locus: other_locus}} by RBH, reference vs each strain."""
    maps = {}
    for s in strains:
        if s == REFERENCE:
            continue
        fwd = blast_best(os.path.join(IN, f"{REFERENCE}.faa"), makedb(s),
                         os.path.join(BLASTDIR, f"val_{REFERENCE}__vs__{s}.tsv"))
        rev = blast_best(os.path.join(IN, f"{s}.faa"), makedb(REFERENCE),
                         os.path.join(BLASTDIR, f"val_{s}__vs__{REFERENCE}.tsv"))
        rev_i = rev.set_index("qseqid")["sseqid"].to_dict()
        maps[s] = {q: t for q, t in zip(fwd.qseqid, fwd.sseqid)
                   if rev_i.get(t) == q}
        print(f"  genome-wide RBH {REFERENCE} <-> {s}: {len(maps[s])} pairs")
    return maps


def score_set(loci, label, ref_cds, cds, maps, window):
    old = s3.WINDOW
    s3.WINDOW = window
    rows = []
    for locus in loci:
        for strain, omap in maps.items():
            sc, _ = s3.score_hub(label, locus, ref_cds, strain, cds[strain], omap)
            if sc:
                sc["class"] = label
                sc["window"] = window
                rows.append(sc)
    s3.WINDOW = old
    return pd.DataFrame(rows)


def main():
    if shutil.which("blastp") is None:
        raise SystemExit("BLAST+ not found on PATH")
    os.makedirs(OUT, exist_ok=True)
    cds = s3.load_cds()
    ref_cds = cds[REFERENCE]
    strains = [s for s in cds if s != REFERENCE]
    maps = genome_wide_rbh(list(cds))

    hub_ids = s3.load_hub_ids() if hasattr(s3, "load_hub_ids") else pd.DataFrame()
    hub_loci = sorted(set(hub_ids["locus_tag"])) if len(hub_ids) else []
    prod = ref_cds["product"].fillna("").str.lower()
    pos = sorted(ref_cds[prod.str.contains("|".join(POSITIVE_KEYWORDS))]["locus_tag"])
    neg = sorted(ref_cds[prod.str.contains("|".join(
        [k.lower() for k in NEGATIVE_KEYWORDS]))]["locus_tag"])
    pool = [l for l in ref_cds["locus_tag"]
            if l not in set(hub_loci) | set(pos) | set(neg)]
    rng = random.Random(SEED)
    background = rng.sample(pool, min(N_BACKGROUND, len(pool)))
    print(f"  hubs {len(hub_loci)} | background {len(background)} | "
          f"ribosomal {len(pos)} | mobile {len(neg)}")

    # --- 1-2. background and controls at the analysis window
    bg = score_set(background, "background", ref_cds, cds, maps, s3.WINDOW)
    ctl_pos = score_set(pos, "ribosomal_positive_control", ref_cds, cds, maps, s3.WINDOW)
    ctl_neg = score_set(neg, "mobile_negative_control", ref_cds, cds, maps, s3.WINDOW)
    hubs = score_set(hub_loci, "hub", ref_cds, cds, maps, s3.WINDOW)
    bg.to_csv(os.path.join(OUT, "background_scores.csv"), index=False)
    pd.concat([ctl_pos, ctl_neg]).to_csv(os.path.join(OUT, "control_scores.csv"),
                                         index=False)

    bg_mean = bg.groupby("reference_locus")["S_syn"].mean()
    hub_mean = hubs.groupby("reference_locus")["S_syn"].mean()
    bgv = bg_mean.values
    comp = pd.DataFrame({
        "reference_locus": hub_mean.index,
        "hub_mean_S_syn": hub_mean.values.round(4),
        "background_mean_S_syn": round(float(bgv.mean()), 4),
        "background_median": round(float(np.median(bgv)), 4),
        "background_sd": round(float(bgv.std(ddof=1)), 4),
        "p_high": [round(float((bgv >= v).mean()), 4) for v in hub_mean.values],
        "p_low": [round(float((bgv <= v).mean()), 4) for v in hub_mean.values],
    })
    if len(hub_ids):
        comp = comp.merge(hub_ids[["gene", "locus_tag"]]
                          .rename(columns={"locus_tag": "reference_locus"}),
                          on="reference_locus", how="left")
    comp.sort_values("hub_mean_S_syn").to_csv(
        os.path.join(OUT, "hub_vs_background.csv"), index=False)

    # --- 3. window sensitivity
    ws = pd.concat([score_set(hub_loci, "hub", ref_cds, cds, maps, w)
                    for w in WINDOWS])
    ws.to_csv(os.path.join(OUT, "window_sensitivity.csv"), index=False)
    wpiv = ws.groupby(["reference_locus", "window"])["S_syn"].mean().unstack()

    # --- 4. self-consistency
    self_map = {l: l for l in ref_cds["locus_tag"]}
    self_rows = []
    for locus in (hub_loci or list(ref_cds["locus_tag"])[:50]):
        sc, _ = s3.score_hub("self", locus, ref_cds, REFERENCE, ref_cds, self_map)
        if sc:
            self_rows.append(sc)
    selfdf = pd.DataFrame(self_rows)
    selfdf.to_csv(os.path.join(OUT, "self_consistency.csv"), index=False)
    self_ok = bool((selfdf["S_syn"] == 1.0).all()) if len(selfdf) else None

    L = ["# Microsynteny validation", "",
         f"Reference {REFERENCE}; {len(strains)} comparison strains; "
         f"window +/-{s3.WINDOW} CDS for the main table.", "",
         "## 1. Hubs against the genome background", "",
         "| class | n loci | mean S_syn | median | sd |", "|---|---|---|---|---|"]
    for lbl, d in [("hub", hubs), ("background", bg),
                   ("ribosomal (positive control)", ctl_pos),
                   ("mobile element (negative control)", ctl_neg)]:
        if len(d):
            m = d.groupby("reference_locus")["S_syn"].mean()
            L.append(f"| {lbl} | {len(m)} | {m.mean():.3f} | "
                     f"{m.median():.3f} | {m.std(ddof=1):.3f} |")
    L += ["", "Interpretation: the hub row is only informative relative to the "
          "background row. If the two are equal, hub neighbourhoods are no more "
          "conserved than an average gene, and S_syn should be reported as a "
          "descriptive statistic, not as evidence of selection.", "",
          "## 2. Controls", "",
          "Ribosomal loci should score near 1 and mobile-element loci well "
          "below the background. If they do not separate, the metric is not "
          "measuring what it claims.", "",
          "## 3. Window sensitivity (mean S_syn per window)", "",
          "| window | mean hub S_syn |", "|---|---|"]
    for w in WINDOWS:
        if w in wpiv.columns:
            L.append(f"| +/-{w} | {wpiv[w].mean():.3f} |")
    L += ["", "## 4. Self-consistency", "",
          f"Reference scored against itself: all S_syn == 1.0 -> {self_ok}", ""]
    with open(os.path.join(OUT, "validation_report.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print(f"  wrote {OUT}/validation_report.md")


if __name__ == "__main__":
    main()
