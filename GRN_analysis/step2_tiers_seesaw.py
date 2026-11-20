"""
step2_tiers_seesaw.py -- Step 2: three-tier hub stratification and
mutual-exclusivity ("seesaw") analysis across the seven GRNs.


Outputs (results/step2_tiers_seesaw/)
  tier_assignments_all_regulators.csv  every regulator, its tier and Z_out
  tier_counts_per_network.csv          Tier 1/2/3 counts per condition
  hub_tier_table.csv                   the condition's own hubs, tiered
  hub_tier_matrix.csv                  hub x condition tier matrix
  hub_kout_matrix.csv / hub_kin_matrix.csv
  dominance_per_network.csv            who leads each network and by how much
  seesaw_pairs.csv                     all hub pairs, seesaw flag + Spearman
  step2_report.md                      tables + verdict

Run:  python step2_tiers_seesaw.py
"""
from __future__ import annotations

import os
from itertools import combinations

import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr

from grn_common import (load_all, regulators, degree_frame, outdir, write,
                        HUBS_BY_CONDITION, hubs_for)

SUB = "step2_tiers_seesaw"
NEAR_APEX_KIN = 1     # k_in <= this, with k_in > 0, counts as near-apex
Z_CUT = 2.0           # out-degree Z-score for an outlier hub
RANK_CUT = 10         # "top-ranked regulator" for the seesaw rule
LOW_FRAC = 0.25       # "suppressed" = <= this fraction of the hub's own max

# editable annotation: drug class / mechanism of action
MECHANISM = {
    "amikacin":      "aminoglycoside - 30S protein synthesis",
    "ciprofloxacin": "fluoroquinolone - DNA gyrase/topoisomerase",
    "doripenem":     "carbapenem - cell wall (PBP)",
    "meropenem":     "carbapenem - cell wall (PBP)",
    "polymyxin":     "polymyxin - outer membrane",
    "tigecycline":   "glycylcycline - 30S protein synthesis",
    "colistin_res":  "polymyxin - outer membrane (resistance)",
}

ALL_HUBS = sorted({h for hs in HUBS_BY_CONDITION.values() for h in hs})


# --------------------------------------------------------------- 2a. tiers
def tier_of(k_in: int, k_out_tf: int) -> str:
    if k_in == 0:
        return "Tier1_apex"
    if k_out_tf > 0:
        return "Tier2_broadcaster"
    return "Tier3_convergence_sink"


def tier_frame(G: nx.DiGraph) -> pd.DataFrame:
    """Tier + Z_out for every active regulator of one network."""
    cond = G.graph.get("name", "")
    df = degree_frame(G)
    tf = df[df["is_regulator"]].copy()
    mu, sd = tf["k_out"].mean(), tf["k_out"].std(ddof=1)
    tf["tier"] = [tier_of(ki, kt) for ki, kt in zip(tf["k_in"], tf["k_out_to_tf"])]
    tf["near_apex"] = (tf["k_in"] > 0) & (tf["k_in"] <= NEAR_APEX_KIN)
    tf["z_out"] = (tf["k_out"] - mu) / sd if sd and sd == sd else 0.0
    tf["outlier_hub_z_ge_cut"] = tf["z_out"] >= Z_CUT
    tf["kout_rank"] = tf["k_out"].rank(ascending=False, method="min").astype(int)
    tf["condition"] = cond
    tf["mechanism"] = MECHANISM.get(cond, "")
    tf["is_declared_hub"] = tf["node"].isin(hubs_for(cond))
    return tf[["condition", "mechanism", "node", "tier", "near_apex",
               "k_in", "k_out", "k_out_to_tf", "k_out_to_leaf", "z_out",
               "outlier_hub_z_ge_cut", "kout_rank", "is_declared_hub"]]


# ------------------------------------------------------- 2b. dominance
def dominance(G: nx.DiGraph) -> dict:
    """Who leads this network, and how concentrated is its out-degree?"""
    cond = G.graph.get("name", "")
    tfs = regulators(G)
    kout = sorted(((G.out_degree(t), t) for t in tfs), reverse=True)
    total = sum(k for k, _ in kout)
    top = kout[0] if kout else (0, "")
    second = kout[1] if len(kout) > 1 else (0, "")
    return {
        "condition": cond,
        "mechanism": MECHANISM.get(cond, ""),
        "declared_hubs": ";".join(hubs_for(cond)),
        "top_regulator": top[1], "top_k_out": top[0],
        "second_regulator": second[1], "second_k_out": second[0],
        "lead_margin": top[0] - second[0],
        "top_share_of_all_out_edges": round(top[0] / total, 4) if total else 0.0,
        "n_active_regulators": len(tfs),
        "hubs_absent_from_this_network": ";".join(
            h for h in ALL_HUBS if h not in G),
    }


# --------------------------------------------------------- 2c. seesaw
def hub_matrices(graphs: dict) -> tuple:
    conds = list(graphs)
    kout = pd.DataFrame(0, index=ALL_HUBS, columns=conds, dtype=int)
    kin = pd.DataFrame(0, index=ALL_HUBS, columns=conds, dtype=int)
    rank = pd.DataFrame(np.nan, index=ALL_HUBS, columns=conds)
    tier = pd.DataFrame("absent", index=ALL_HUBS, columns=conds, dtype=object)
    for cond, G in graphs.items():
        tf = tier_frame(G).set_index("node")
        for h in ALL_HUBS:
            if h in G:
                kout.loc[h, cond] = G.out_degree(h)
                kin.loc[h, cond] = G.in_degree(h)
            if h in tf.index:
                rank.loc[h, cond] = tf.at[h, "kout_rank"]
                tier.loc[h, cond] = tf.at[h, "tier"]
            elif h in G:
                tier.loc[h, cond] = "pure_target"
    return kout, kin, rank, tier


def seesaw_table(kout: pd.DataFrame, rank: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for a, b in combinations(ALL_HUBS, 2):
        ka, kb = kout.loc[a], kout.loc[b]
        ra, rb = rank.loc[a], rank.loc[b]
        if ka.std() > 0 and kb.std() > 0:
            rho, p = spearmanr(ka.values, kb.values)
        else:
            rho, p = np.nan, np.nan
        a_dom = [c for c in kout.columns
                 if ra[c] <= RANK_CUT and kb[c] <= LOW_FRAC * max(kb.max(), 1)]
        b_dom = [c for c in kout.columns
                 if rb[c] <= RANK_CUT and ka[c] <= LOW_FRAC * max(ka.max(), 1)]
        rows.append({
            "hub_A": a, "hub_B": b,
            "spearman_rho_kout": round(float(rho), 3) if rho == rho else np.nan,
            "p_value": round(float(p), 4) if p == p else np.nan,
            "n_cond_A_dominant_B_low": len(a_dom),
            "n_cond_B_dominant_A_low": len(b_dom),
            "seesaw": bool(a_dom) and bool(b_dom),
            "A_dominant_in": ";".join(a_dom),
            "B_dominant_in": ";".join(b_dom),
            "A_mechanisms": ";".join(sorted({MECHANISM.get(c, "") for c in a_dom})),
            "B_mechanisms": ";".join(sorted({MECHANISM.get(c, "") for c in b_dom})),
        })
    df = pd.DataFrame(rows)
    df["seesaw_score"] = df["n_cond_A_dominant_B_low"] + df["n_cond_B_dominant_A_low"]
    return df.sort_values(["seesaw", "seesaw_score"], ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------- report
def write_report(tiers, counts, hubtab, dom, see, tiermat, path):
    L = ["# Step 2 - hub tiers and seesaw dominance", "",
         "## Tier composition per network", "",
         "| condition | mechanism | Tier 1 apex | Tier 2 broadcaster | "
         "Tier 3 sink | outlier hubs (Z>=%.1f) |" % Z_CUT,
         "|---|---|---|---|---|---|"]
    for _, r in counts.iterrows():
        L.append(f"| {r['condition']} | {MECHANISM.get(r['condition'],'')} | "
                 f"{r.get('Tier1_apex',0)} | {r.get('Tier2_broadcaster',0)} | "
                 f"{r.get('Tier3_convergence_sink',0)} | {r['n_outlier_hubs']} |")

    L += ["", "## Declared hubs, tiered", "",
          "| condition | hub | tier | k_in | k_out | k_out to TFs | Z_out | "
          "rank by k_out |", "|---|---|---|---|---|---|---|---|"]
    for _, r in hubtab.iterrows():
        L.append(f"| {r['condition']} | {r['node']} | {r['tier']} | {r['k_in']} | "
                 f"{r['k_out']} | {r['k_out_to_tf']} | {r['z_out']:.2f} | "
                 f"{r['kout_rank']} |")

    L += ["", "## Dominance per network", "",
          "| condition | mechanism | leading regulator | k_out | runner-up | "
          "margin | share of out-edges |", "|---|---|---|---|---|---|---|"]
    for _, r in dom.iterrows():
        L.append(f"| {r['condition']} | {r['mechanism']} | {r['top_regulator']} | "
                 f"{r['top_k_out']} | {r['second_regulator']} ({r['second_k_out']}) | "
                 f"{r['lead_margin']} | {100*r['top_share_of_all_out_edges']:.1f}% |")

    ss = see[see["seesaw"]]
    L += ["", "## Seesaw pairs", "",
          f"Hub pairs tested: {len(see)}; pairs meeting the seesaw rule: {len(ss)}.", ""]
    if len(ss):
        L += ["| hub A | hub B | rho | A dominant in | B dominant in |",
              "|---|---|---|---|---|"]
        for _, r in ss.head(20).iterrows():
            L.append(f"| {r['hub_A']} | {r['hub_B']} | {r['spearman_rho_kout']} | "
                     f"{r['A_dominant_in']} | {r['B_dominant_in']} |")
    cols = list(tiermat.columns)
    L += ["", "## Tier matrix (hub x condition)", "",
          "| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    for _, r in tiermat.iterrows():
        L.append("| " + " | ".join(str(r[c]) for c in cols) + " |")
    L += ["",
          "Note: MECHANISM is manual annotation; every other column is derived "
          "from the edge files. A hub absent from a network is reported as "
          "'absent', which is not the same as a hub present with k_out = 0 "
          "('pure_target').", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print(f"  wrote {path}")


def main():
    graphs = load_all()
    tiers = pd.concat([tier_frame(G) for G in graphs.values()], ignore_index=True)
    dom = pd.DataFrame([dominance(G) for G in graphs.values()])

    counts = (tiers.groupby(["condition", "tier"]).size().unstack(fill_value=0)
              .reset_index())
    counts["n_outlier_hubs"] = (tiers[tiers["outlier_hub_z_ge_cut"]]
                                .groupby("condition").size()
                                .reindex(counts["condition"]).fillna(0).astype(int).values)

    hubtab = tiers[tiers["is_declared_hub"]].copy()
    kout, kin, rank, tier = hub_matrices(graphs)
    see = seesaw_table(kout, rank)

    write(tiers, "tier_assignments_all_regulators.csv", SUB)
    write(counts, "tier_counts_per_network.csv", SUB)
    write(hubtab, "hub_tier_table.csv", SUB)
    write(kout.reset_index().rename(columns={"index": "hub"}), "hub_kout_matrix.csv", SUB)
    write(kin.reset_index().rename(columns={"index": "hub"}), "hub_kin_matrix.csv", SUB)
    tiermat = tier.reset_index().rename(columns={"index": "hub"})
    write(tiermat, "hub_tier_matrix.csv", SUB)
    write(dom, "dominance_per_network.csv", SUB)
    write(see, "seesaw_pairs.csv", SUB)
    write_report(tiers, counts, hubtab, dom, see, tiermat,
                 os.path.join(outdir(SUB), "step2_report.md"))


if __name__ == "__main__":
    main()
