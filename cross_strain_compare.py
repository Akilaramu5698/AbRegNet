#!/usr/bin/env python3
"""
Cross-strain motif comparison script

Can be run from THREE locations — auto-detects which:
  1. MOTIF root:       cd MOTIF && python3 cross_strain_compare.py
  2. Antibiotic folder: cd MOTIF/Amikacin && python3 cross_strain_compare.py
  3. Strain folder:     cd MOTIF/Amikacin/gse56218 && python3 cross_strain_compare.py

Output: cross_strain_results/ folder with all comparison tables.
"""

import os, sys, glob
import pandas as pd

OUT = "cross_strain_results"
os.makedirs(OUT, exist_ok=True)

print("=" * 60)
print(" Cross-Strain Motif Comparison")
print("=" * 60)
print(f"\n  Folder: {os.getcwd()}")


# AUTO-DETECT run location and build file list accordingly

def find_files(pattern_list):
    for pattern in pattern_list:
        files = sorted(glob.glob(pattern))
        if files:
            return files, pattern
    return [], None

# Try 3 levels deep → 2 levels → 1 level
enrich_files, e_pat = find_files([
    "*/*/results/motif_enrichment.tsv",   # MOTIF root
    "*/results/motif_enrichment.tsv",     # inside antibiotic folder
    "results/motif_enrichment.tsv",       # inside strain folder
])
edges_files, d_pat = find_files([
    "*/*/results/final_edges.tsv",
    "*/results/final_edges.tsv",
    "results/final_edges.tsv",
])

if not enrich_files and not edges_files:
    print("\n  ERROR: No result files found.")
    print("  Run this from your MOTIF root, antibiotic, or strain folder.")
    sys.exit(1)

# Determine the path parser based on detected level
def parse_path(fpath, pattern):
    """Extract (antibiotic, strain) from file path based on glob pattern depth."""
    parts = fpath.replace("\\", "/").split("/")
    if pattern.startswith("*/*/"):
        # MOTIF root: parts = [antibiotic, strain, results, file]
        return parts[0], parts[1]
    elif pattern.startswith("*/"):
        # Antibiotic folder: parts = [strain, results, file]
        antibiotic = os.path.basename(os.getcwd())
        return antibiotic, parts[0]
    else:
        # Strain folder: parts = [results, file]
        strain    = os.path.basename(os.getcwd())
        antibiotic = os.path.basename(os.path.dirname(os.getcwd()))
        return antibiotic, strain

mode = "MOTIF root" if e_pat and e_pat.startswith("*/*/") else \
       "single antibiotic" if e_pat and e_pat.startswith("*/") else \
       "single strain"
print(f"  Mode detected: {mode}")
print(f"  Found: {len(enrich_files)} enrichment files, {len(edges_files)} edge files")


# LOAD enrichment files

print("\n[1] Loading motif_enrichment.csv files...")
enrich_frames = []
for fpath in enrich_files:
    ab, strain = parse_path(fpath, e_pat)
    try:
        df = pd.read_csv(fpath, sep="\t")
        df.insert(0, "Antibiotic", ab)
        df.insert(1, "Strain",     strain)
        df.insert(2, "Condition",  f"{ab}/{strain}")
        enrich_frames.append(df)
        print(f"  ✓ {ab}/{strain}: {len(df)} rows")
    except Exception as e:
        print(f"  ✗ {ab}/{strain}: {e}")

enrich_all = pd.concat(enrich_frames, ignore_index=True) if enrich_frames else pd.DataFrame()


# LOAD edges files

print("\n[2] Loading final_edges.csv files...")
edges_frames = []
for fpath in edges_files:
    ab, strain = parse_path(fpath, d_pat)
    try:
        df = pd.read_csv(fpath, sep="\t")
        df.insert(0, "Antibiotic", ab)
        df.insert(1, "Strain",     strain)
        df.insert(2, "Condition",  f"{ab}/{strain}")
        edges_frames.append(df)
        print(f"  ✓ {ab}/{strain}: {len(df)} edges")
    except Exception as e:
        print(f"  ✗ {ab}/{strain}: {e}")

edges_all = pd.concat(edges_frames, ignore_index=True) if edges_frames else pd.DataFrame()


# TABLE 1 — Regulator enrichment per strain

print("\n[3] Building Table 1: Regulator enrichment per strain...")

if not enrich_all.empty:
    records = []
    for condition in sorted(enrich_all["Condition"].unique()):
        ab, strain = condition.split("/", 1)
        cdf = enrich_all[enrich_all["Condition"] == condition]
        for reg in sorted(cdf["Regulator"].unique()):
            rdf  = cdf[cdf["Regulator"] == reg]
            sig  = rdf[rdf["Significant"] == True]
            best = rdf.loc[rdf["FDR"].idxmin()] if not rdf.empty else None
            records.append({
                "Antibiotic":            ab,
                "Strain":                strain,
                "Regulator":             reg,
                "Total_motifs_tested":   len(rdf),
                "Sig_motifs_FDR05":      len(sig),
                "Best_FDR":              round(best["FDR"], 6) if best is not None else None,
                "Best_fold_enrichment":  round(best["Fold_enrichment"], 2) if best is not None else None,
                "Best_motif_pattern":    best["Motif"] if best is not None else None,
                "Has_sig_motif":         len(sig) > 0,
            })
    t1 = pd.DataFrame(records)
    t1.to_csv(f"{OUT}/Table1_regulator_enrichment_per_strain.csv", sep=",", index=False)
    print(f"  ✓ Saved Table1_regulator_enrichment_per_strain.csv  ({len(t1)} rows, "
          f"{t1['Antibiotic'].nunique()} antibiotics)")

    # Pivot: antibiotic+regulator rows × strain columns
    try:
        t1b = t1.pivot_table(
            index=["Antibiotic", "Regulator"],
            columns="Strain",
            values="Sig_motifs_FDR05",
            aggfunc="sum"
        ).fillna(0).astype(int).reset_index()
        t1b.to_csv(f"{OUT}/Table1b_sig_motifs_pivot.csv", sep=",", index=False)
        print(f"  ✓ Saved Table1b_sig_motifs_pivot.csv")
    except Exception as e:
        print(f"  ⚠ Could not build pivot: {e}")


# TABLE 2 — Edge motif consistency across strains PER ANTIBIOTIC

print("\n[4] Building Table 2: Edge motif consistency per antibiotic...")

if not edges_all.empty:
    for antibiotic in sorted(edges_all["Antibiotic"].unique()):
        adf     = edges_all[edges_all["Antibiotic"] == antibiotic]
        strains = sorted(adf["Strain"].unique())
        n       = len(strains)
        records = []

        for (reg, tgt), grp in adf.groupby(["Regulator", "Target"]):
            row   = {"Regulator": reg, "Target": tgt}
            n_hit = 0
            for strain in strains:
                sdf     = grp[grp["Strain"] == strain]
                has_hit = (not sdf.empty) and (sdf["Motif_hit"] == "Yes").any()
                best_p  = sdf[sdf["Motif_hit"] == "Yes"]["Motif_pval"].iloc[0] \
                          if has_hit else "—"
                row[f"{strain}_hit"]  = "Yes" if has_hit else ("No" if not sdf.empty else "absent")
                row[f"{strain}_pval"] = best_p
                if has_hit:
                    n_hit += 1

            row["Strains_with_hit"] = n_hit
            row["Total_strains"]    = n
            row["Consistency_pct"]  = f"{round(n_hit/n*100)}%"
            row["Conserved_all"]    = n_hit == n
            records.append(row)

        t2 = pd.DataFrame(records).sort_values(
            ["Regulator", "Strains_with_hit"], ascending=[True, False]
        )
        ab_safe = antibiotic.replace(" ", "_").replace("/", "_")
        t2.to_csv(f"{OUT}/Table2_{ab_safe}_edge_consistency.csv", sep=",", index=False)
        n_cons = t2["Conserved_all"].sum()
        print(f"  ✓ {antibiotic}: {n_cons}/{len(t2)} edges conserved across {n} strain(s)")


# TABLE 3 — Conserved edges only (motif hit in ALL strains)

print("\n[5] Building Table 3: Conserved edges (motif in ALL strains)...")

if not edges_all.empty:
    conserved = []
    for antibiotic in sorted(edges_all["Antibiotic"].unique()):
        adf     = edges_all[edges_all["Antibiotic"] == antibiotic]
        strains = sorted(adf["Strain"].unique())
        n       = len(strains)
        for (reg, tgt), grp in adf.groupby(["Regulator", "Target"]):
            hits = {s: ((not grp[grp["Strain"]==s].empty) and
                        (grp[grp["Strain"]==s]["Motif_hit"] == "Yes").any())
                    for s in strains}
            if all(hits.values()):
                conserved.append({
                    "Antibiotic":    antibiotic,
                    "Regulator":     reg,
                    "Target":        tgt,
                    "Strains_total": n,
                    "Strains_hit":   n,
                    "Conserved":     "YES — all strains",
                })
    if conserved:
        t3 = pd.DataFrame(conserved).sort_values(["Antibiotic","Regulator","Target"])
        t3.to_csv(f"{OUT}/Table3_conserved_edges_all_strains.csv", sep=",", index=False)
        print(f"  ✓ Saved Table3_conserved_edges_all_strains.csv")
        for ab, grp in t3.groupby("Antibiotic"):
            print(f"    {ab}: {len(grp)} conserved edges ({grp['Strains_total'].iloc[0]} strains)")
    else:
        print("  No fully conserved edges found")


# TABLE 4 — Regulator consistency score

print("\n[6] Building Table 4: Regulator consistency score...")

if not enrich_all.empty:
    rec = []
    for (ab, reg), grp in enrich_all.groupby(["Antibiotic", "Regulator"]):
        strains  = grp["Strain"].unique()
        n_total  = len(strains)
        n_sig    = sum(grp[grp["Strain"]==s]["Significant"].any() for s in strains)
        sig_df   = grp[grp["Significant"] == True]
        rec.append({
            "Antibiotic":             ab,
            "Regulator":              reg,
            "Strains_total":          n_total,
            "Strains_sig_FDR05":      n_sig,
            "Consistency_pct":        f"{round(n_sig/n_total*100)}%",
            "Consistently_sig_all":   n_sig == n_total,
            "Best_FDR_overall":       round(sig_df["FDR"].min(), 6) if not sig_df.empty else "—",
            "Avg_fold_enrichment_sig":round(sig_df["Fold_enrichment"].mean(), 2) if not sig_df.empty else "—",
        })
    t4 = pd.DataFrame(rec).sort_values(
        ["Antibiotic", "Strains_sig_FDR05"], ascending=[True, False]
    )
    t4.to_csv(f"{OUT}/Table4_regulator_consistency_score.csv", sep=",", index=False)
    print(f"  ✓ Saved Table4_regulator_consistency_score.csv  ({len(t4)} rows)")
    print(f"\n  Regulator consistency summary:")
    print(f"  {'Antibiotic':<22} {'Regulator':<10} {'Sig strains':<14} {'Consistent?'}")
    print(f"  {'-'*60}")
    for _, row in t4.iterrows():
        flag = "★ ALL strains" if row["Consistently_sig_all"] else ""
        print(f"  {row['Antibiotic']:<22} {row['Regulator']:<10} "
              f"{row['Strains_sig_FDR05']}/{row['Strains_total']} ({row['Consistency_pct']:<8}) {flag}")


# TABLE 5 — All significant motifs catalogue

print("\n[7] Building Table 5: All significant motifs catalogue...")

if not enrich_all.empty:
    t5 = enrich_all[enrich_all["Significant"] == True][
        ["Antibiotic","Strain","Regulator","Motif",
         "Target_with_motif","Target_total",
         "Fold_enrichment","p_value","FDR"]
    ].sort_values(["Antibiotic","Regulator","FDR"])
    t5.to_csv(f"{OUT}/Table5_all_significant_motifs.csv", sep=",", index=False)
    print(f"  ✓ Saved Table5_all_significant_motifs.csv  ({len(t5)} significant motifs)")
    # Summary by antibiotic
    for ab, grp in t5.groupby("Antibiotic"):
        regs = grp.groupby("Regulator").size()
        print(f"    {ab}: {len(grp)} sig motifs — {dict(regs)}")


# FINAL SUMMARY

print("\n" + "=" * 60)
print(" SUMMARY")
print("=" * 60)
if not edges_all.empty:
    n_ab    = edges_all["Antibiotic"].nunique()
    n_cond  = edges_all["Condition"].nunique()
    n_edges = len(edges_all)
    n_hit   = (edges_all["Motif_hit"] == "Yes").sum()
    print(f"\n  Antibiotics:          {n_ab} ({', '.join(sorted(edges_all['Antibiotic'].unique()))})")
    print(f"  Strain conditions:    {n_cond}")
    print(f"  Total edges:          {n_edges}")
    print(f"  Edges with motif hit: {n_hit} ({round(n_hit/n_edges*100,1)}%)")

print(f"\n  Output: {OUT}/")
print("    Table1  — regulator enrichment per strain (FDR, fold, best motif)")
print("    Table1b — sig motifs count pivot (regulator × strain)")
print("    Table2  — per-edge motif hit Yes/No per strain (one file per antibiotic)")
print("    Table3  — edges conserved across ALL strains of same antibiotic")
print("    Table4  — regulator consistency score across strains")
print("    Table5  — all significant motifs (FDR<0.05) with fold enrichment")
print("\n  Done.")
