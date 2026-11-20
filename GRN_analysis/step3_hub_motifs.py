"""
step3_hub_motifs.py -- Step 3: feed-forward loops and bi-fan motifs driven by
the identified hubs.

Outputs (results/step3_hub_motifs/)
  ffl_hub_driven.csv           every FFL involving the hub scope, with roles
  ffl_hub_pair_ranking.csv     (X -> Y) hub pairs ranked by FFLs seeded
  bifan_hub_pairs.csv          co-regulating hub pairs, shared targets, C(s,2)
  hub_cooperation_summary.csv  per hub pair: hierarchical / parallel / both
  motif_counts_per_network.csv totals per condition
  ffl_conservation_by_condition_subset.csv
                               per (X -> Y) seed pair: condition sets sharing
                               at least one downstream Z, with is_maximal
  bifan_conservation_by_condition_subset.csv
                               every set of conditions (size >= 2) in which a
                               pair keeps >= 2 shared targets, with is_maximal
  step3_report.md

Run:  python step3_hub_motifs.py
"""
from __future__ import annotations

import os
from collections import defaultdict
from itertools import combinations
from math import comb

import pandas as pd
import networkx as nx

from grn_common import (load_all, regulators, targets_of, outdir, write,
                        HUBS_BY_CONDITION, hubs_for)

SUB = "step3_hub_motifs"
HUB_SCOPE = "condition"        # "condition" | "all"
REQUIRE_HUB_PAIR = False       # True = motif must be driven by two hubs
MIN_SHARED_TARGETS = 2         # minimum co-regulated targets for a bi-fan
ANNOTATION_FILE = "gene_annotations.csv"   # optional: columns gene,function

ALL_HUBS = sorted({h for hs in HUBS_BY_CONDITION.values() for h in hs})


def hub_set(cond: str) -> set:
    return set(ALL_HUBS) if HUB_SCOPE == "all" else set(hubs_for(cond))


def load_annotations() -> dict:
    if not os.path.exists(ANNOTATION_FILE):
        return {}
    ann = pd.read_csv(ANNOTATION_FILE)
    cols = {c.lower().strip(): c for c in ann.columns}
    g, f = cols.get("gene"), cols.get("function")
    if not (g and f):
        return {}
    return {str(k).strip().lower(): str(v) for k, v in zip(ann[g], ann[f])}


def fun_of(ann: dict, genes) -> str:
    """Functional categories carried by a set of genes, deduplicated."""
    vals = sorted({ann[g] for g in genes if g in ann})
    return ";".join(vals)


# ------------------------------------------------------------------- FFLs
def hub_ffls(G: nx.DiGraph, hubs: set, ann: dict) -> pd.DataFrame:
    """FFL triads X -> Y -> Z (+ X -> Z) involving the hub scope."""
    cond = G.graph.get("name", "")
    tfs = regulators(G)
    succ = {n: targets_of(G, n) for n in G}
    rows = []
    for x, y in G.edges():
        if x == y:
            continue
        for z in succ[y] & succ[x]:
            if z in (x, y):
                continue
            hits = [h for h in (x, y, z) if h in hubs]
            if not hits or (REQUIRE_HUB_PAIR and len({x, y} & hubs) < 2):
                continue
            rows.append({
                "condition": cond, "X": x, "Y": y, "Z": z,
                "X_is_hub": x in hubs, "Y_is_hub": y in hubs,
                "Z_is_hub": z in hubs,
                "n_hubs_in_triad": len(hits),
                "ffl_class": "TF_only" if z in tfs else "TF_target",
                "k_out_X": G.out_degree(x), "k_out_Y": G.out_degree(y),
                "triad": f"{x}->{y}->{z}",
                "function_Z": ann.get(z, ""),
                "functions_in_triad": fun_of(ann, (x, y, z)),
            })
    return pd.DataFrame(rows)


def rank_ffl_pairs(ffl: pd.DataFrame) -> pd.DataFrame:
    """Which (X -> Y) hub pairs seed the most FFLs: the hierarchical couplings."""
    if ffl.empty:
        return pd.DataFrame()
    g = (ffl.groupby(["condition", "X", "Y"])
           .agg(n_ffl=("Z", "size"),
                n_tf_only=("ffl_class", lambda s: int((s == "TF_only").sum())),
                n_tf_target=("ffl_class", lambda s: int((s == "TF_target").sum())),
                both_are_hubs=("X_is_hub", "first"),
                shared_targets=("Z", lambda s: ";".join(sorted(set(s))[:40])),
                functions_covered=("function_Z",
                                   lambda s: ";".join(sorted({v for v in s if v}))))
           .reset_index())
    return g.sort_values(["condition", "n_ffl"], ascending=[True, False]).reset_index(drop=True)


# ---------------------------------------------------------------- bi-fans
def hub_bifans(G: nx.DiGraph, hubs: set, ann: dict) -> pd.DataFrame:
    """Hub pairs co-regulating >= MIN_SHARED_TARGETS targets in parallel."""
    cond = G.graph.get("name", "")
    regs = regulators(G)
    succ = {t: targets_of(G, t) for t in regs}
    by_target = defaultdict(set)
    for t, tg in succ.items():
        for z in tg:
            by_target[z].add(t)

    pairs = defaultdict(set)
    for z, rs in by_target.items():
        if len(rs) > 1:
            for a, b in combinations(sorted(rs), 2):
                pairs[(a, b)].add(z)

    rows = []
    for (a, b), shared in pairs.items():
        hits = [h for h in (a, b) if h in hubs]
        if not hits or (REQUIRE_HUB_PAIR and len(hits) < 2):
            continue
        if len(shared) < MIN_SHARED_TARGETS:
            continue
        union = len(succ[a] | succ[b])
        rows.append({
            "condition": cond, "X1": a, "X2": b,
            "X1_is_hub": a in hubs, "X2_is_hub": b in hubs,
            "n_shared_targets": len(shared),
            "n_bifan_motifs": comb(len(shared), 2),
            "jaccard_targets": round(len(shared) / union, 4) if union else 0.0,
            "k_out_X1": G.out_degree(a), "k_out_X2": G.out_degree(b),
            "edge_between_hubs": G.has_edge(a, b) or G.has_edge(b, a),
            "shared_targets": ";".join(sorted(shared)),
            "functions_coupled": fun_of(ann, shared),
        })
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values(["condition", "n_shared_targets"],
                          ascending=[True, False]).reset_index(drop=True)


# ------------------------------------------------- cooperation + conservation
def cooperation(ffl: pd.DataFrame, bif: pd.DataFrame) -> pd.DataFrame:
    """Per hub pair: hierarchical (FFL), parallel (bi-fan), or both."""
    rows = {}
    if not ffl.empty:
        for _, r in ffl.iterrows():
            key = (r["condition"], *sorted((r["X"], r["Y"])))
            d = rows.setdefault(key, {"n_ffl": 0, "n_shared": 0, "n_bifan": 0,
                                      "edge": True, "functions": set()})
            d["n_ffl"] += 1
            if r["functions_in_triad"]:
                d["functions"] |= set(r["functions_in_triad"].split(";"))
    if not bif.empty:
        for _, r in bif.iterrows():
            key = (r["condition"], *sorted((r["X1"], r["X2"])))
            d = rows.setdefault(key, {"n_ffl": 0, "n_shared": 0, "n_bifan": 0,
                                      "edge": r["edge_between_hubs"],
                                      "functions": set()})
            d["n_shared"] = r["n_shared_targets"]
            d["n_bifan"] = r["n_bifan_motifs"]
            d["edge"] = r["edge_between_hubs"]
            if r["functions_coupled"]:
                d["functions"] |= set(r["functions_coupled"].split(";"))
    out = []
    for (cond, a, b), d in rows.items():
        mode = ("both" if d["n_ffl"] and d["n_bifan"] else
                "hierarchical_FFL" if d["n_ffl"] else "parallel_bifan")
        out.append({"condition": cond, "hub_A": a, "hub_B": b,
                    "cooperation_mode": mode,
                    "direct_edge_between_them": d["edge"],
                    "n_ffls_together": d["n_ffl"],
                    "n_shared_targets": d["n_shared"],
                    "n_bifan_motifs": d["n_bifan"],
                    "functions_coupled": ";".join(sorted(x for x in d["functions"] if x))})
    return (pd.DataFrame(out)
            .sort_values(["condition", "n_bifan_motifs", "n_ffls_together"],
                         ascending=[True, False, False]).reset_index(drop=True))


def ffl_conservation_by_condition_subset(ffl: pd.DataFrame,
                                         min_core: int = 1) -> pd.DataFrame:
    """For each (X -> Y) seed pair: the condition sets sharing >= min_core
    downstream targets Z, i.e. where a complete FFL survives.

    The bi-fan version needs two core genes for a square; an FFL needs one, so
    min_core = 1 here. `is_maximal` marks subsets not contained in a larger
    qualifying one, which answers "across how many antibiotics does this FFL
    hold?".
    """
    if ffl.empty:
        return pd.DataFrame()
    by = (ffl.groupby(["X", "Y", "condition"])["Z"]
            .apply(lambda s: set(s)).reset_index())
    rows = []
    for (x, y), g in by.groupby(["X", "Y"]):
        sets = dict(zip(g["condition"], g["Z"]))
        conds = sorted(sets)
        if len(conds) < 2:
            continue
        qualifying = []
        for k in range(2, len(conds) + 1):
            for sub in combinations(conds, k):
                core = set.intersection(*(sets[c] for c in sub))
                if len(core) >= min_core:
                    qualifying.append((set(sub), core))
        for sub, core in qualifying:
            rows.append({
                "seed_pair": f"{x} -> {y}", "X": x, "Y": y,
                "n_conditions_in_subset": len(sub),
                "conditions": ";".join(sorted(sub)),
                "n_conditions_pair_occurs_in": len(conds),
                "n_core_targets_Z": len(core),
                "core_targets_Z": ";".join(sorted(core)),
                "is_maximal": not any(sub < o for o, _ in qualifying),
                "per_condition_Z_counts": ";".join(
                    f"{c}:{len(sets[c])}" for c in sorted(sub))})
    out = pd.DataFrame(rows)
    return (out.sort_values(["n_conditions_in_subset", "n_core_targets_Z",
                             "seed_pair"], ascending=[False, False, True])
            .reset_index(drop=True) if len(out) else out)


def conservation_by_condition_subset(bif: pd.DataFrame,
                                     min_core: int = 2) -> pd.DataFrame:
    """Every SET of conditions (size >= 2) in which a pair keeps >= min_core
    shared targets, so a complete four-node bi-fan survives across that set.

    Condition pairs are only the k = 2 case. With at most 7 conditions a pair
    can occur in, all 2^k subsets are enumerated exactly. `is_maximal` marks
    the subsets that are not contained in a larger qualifying subset - those
    are the answer to "across how many antibiotics does this bi-fan hold?".
    """
    if bif.empty:
        return pd.DataFrame()
    df = bif.copy()
    df["_t"] = df["shared_targets"].fillna("").apply(
        lambda s: {x for x in str(s).split(";") if x})
    df["_pair"] = [tuple(sorted((a, b))) for a, b in zip(df["X1"], df["X2"])]
    rows = []
    for (a, b), g in df.groupby("_pair"):
        sets = dict(zip(g["condition"], g["_t"]))
        conds = sorted(sets)
        if len(conds) < 2:
            continue
        qualifying = []
        for k in range(2, len(conds) + 1):
            for sub in combinations(conds, k):
                core = set.intersection(*(sets[c] for c in sub))
                if len(core) >= min_core:
                    qualifying.append((set(sub), core))
        for sub, core in qualifying:
            maximal = not any(sub < other for other, _ in qualifying)
            rows.append({
                "pair": f"{a} + {b}", "X1": a, "X2": b,
                "n_conditions_in_subset": len(sub),
                "conditions": ";".join(sorted(sub)),
                "n_conditions_pair_occurs_in": len(conds),
                "n_core_genes": len(core),
                "core_genes": ";".join(sorted(core)),
                "n_conserved_bifan_motifs": comb(len(core), 2),
                "is_maximal": maximal,
                "per_condition_shared_counts": ";".join(
                    f"{c}:{len(sets[c])}" for c in sorted(sub))})
    out = pd.DataFrame(rows)
    return (out.sort_values(["n_conditions_in_subset", "n_core_genes", "pair"],
                            ascending=[False, False, True]).reset_index(drop=True)
            if len(out) else out)


def write_report(counts, coop, bif, path):
    L = ["# Step 3 - hub-driven feed-forward loops and bi-fans", "",
         f"Hub scope: {HUB_SCOPE}; motifs must involve at least "
         f"{'two hubs' if REQUIRE_HUB_PAIR else 'one hub'}.", "",
         "## Motif counts per network", "",
         "| condition | hub FFLs | TF-only | TF-target | hub bi-fan pairs | "
         "bi-fan motifs |", "|---|---|---|---|---|---|"]
    for _, r in counts.iterrows():
        L.append(f"| {r['condition']} | {r['n_hub_ffl']} | {r['n_ffl_tf_only']} | "
                 f"{r['n_ffl_tf_target']} | {r['n_bifan_pairs']} | "
                 f"{r['n_bifan_motifs']} |")
    L += ["", "## How hub pairs cooperate", "",
          "| condition | hub A | hub B | mode | direct edge | FFLs | "
          "shared targets | bi-fans |", "|---|---|---|---|---|---|---|---|"]
    for _, r in coop.head(40).iterrows():
        L.append(f"| {r['condition']} | {r['hub_A']} | {r['hub_B']} | "
                 f"{r['cooperation_mode']} | {r['direct_edge_between_them']} | "
                 f"{r['n_ffls_together']} | {r['n_shared_targets']} | "
                 f"{r['n_bifan_motifs']} |")
    L += ["", "Reading: 'hierarchical_FFL' means one hub gates the other's "
          "output (signal filter); 'parallel_bifan' means the two hubs hit the "
          "same targets with no edge between them (combinatorial switch); "
          "'both' means the pair does both at once. Functional categories are "
          "present only if gene_annotations.csv was supplied.", ""]
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print(f"  wrote {path}")


def main():
    graphs = load_all()
    ann = load_annotations()
    print(f"  annotations loaded: {len(ann)} genes"
          if ann else "  no gene_annotations.csv - function columns left blank")

    ffls, bifs, counts = [], [], []
    for cond, G in graphs.items():
        hubs = hub_set(cond)
        f = hub_ffls(G, hubs, ann)
        b = hub_bifans(G, hubs, ann)
        ffls.append(f)
        bifs.append(b)
        counts.append({
            "condition": cond,
            "hub_scope": ";".join(sorted(hubs)),
            "n_hub_ffl": len(f),
            "n_ffl_tf_only": int((f["ffl_class"] == "TF_only").sum()) if len(f) else 0,
            "n_ffl_tf_target": int((f["ffl_class"] == "TF_target").sum()) if len(f) else 0,
            "n_bifan_pairs": len(b),
            "n_bifan_motifs": int(b["n_bifan_motifs"].sum()) if len(b) else 0,
            "max_shared_targets": int(b["n_shared_targets"].max()) if len(b) else 0,
        })

    cat = lambda fs: (pd.concat([x for x in fs if len(x)], ignore_index=True)
                      if any(len(x) for x in fs) else pd.DataFrame())
    ffl_all, bif_all = cat(ffls), cat(bifs)
    counts = pd.DataFrame(counts)
    coop = cooperation(ffl_all, bif_all)

    write(counts, "motif_counts_per_network.csv", SUB)
    write(ffl_all, "ffl_hub_driven.csv", SUB)
    write(rank_ffl_pairs(ffl_all), "ffl_hub_pair_ranking.csv", SUB)
    write(bif_all, "bifan_hub_pairs.csv", SUB)
    write(coop, "hub_cooperation_summary.csv", SUB)
    fsub = ffl_conservation_by_condition_subset(ffl_all)
    write(fsub, "ffl_conservation_by_condition_subset.csv", SUB)
    if len(fsub):
        fm = fsub[fsub["is_maximal"]]
        print(f"  FFL seed pairs with a shared target across conditions: "
              f"{fsub['seed_pair'].nunique()}; {len(fm)} maximal subsets; "
              f"largest = {int(fsub['n_conditions_in_subset'].max())} conditions")
    else:
        print("  no FFL survives across any two conditions")
    sub = conservation_by_condition_subset(bif_all)
    write(sub, "bifan_conservation_by_condition_subset.csv", SUB)
    if len(sub):
        mx = sub[sub["is_maximal"]]
        print(f"  condition subsets with >= 2 core genes: {len(sub)} "
              f"({len(mx)} maximal) over {sub['pair'].nunique()} regulator "
              f"pairs; largest subset = "
              f"{int(sub['n_conditions_in_subset'].max())} conditions")
    else:
        print("  no condition subset keeps >= 2 core genes - no four-node "
              "bi-fan survives between any two antibiotics")
    write_report(counts, coop, bif_all, os.path.join(outdir(SUB), "step3_report.md"))


if __name__ == "__main__":
    main()
