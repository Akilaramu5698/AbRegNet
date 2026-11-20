"""
coverage_leaves.py -- leaf coverage of the significant hub motifs.

  FFL leaves          unique leaves over all Step 1 FFLs
  FFL hub leaves      unique leaves over the Step 1 FFLs that also hold a hub
  Bi-fan leaves       unique leaves over all Step 1 pure bi-fans
  Bi-fan hub leaves   unique leaves over the Step 1 pure bi-fans with a hub
  n_leaves            | FFL leaves UNION bi-fan leaves |
  n_hub_leaves        | FFL hub leaves UNION bi-fan hub leaves |
  coverage            n_hub_leaves / n_leaves


Significance follows the pipeline's two-stage flags: "significant" = Step 1
(q < 0.05 and Z >= 2.0); "hub" = Step 2 (Step 1 plus at least one Master Hub).

Run:  python coverage_leaves.py [results/step4_motif_enrichment]
      (run it from the directory holding grn_common.py and the *.grn.csv files)
"""
from __future__ import annotations

import os
import sys

import pandas as pd

from grn_common import load_all, regulators, outdir, write, hubs_for

# Which hubs count as "hub" for the coverage numerator.
#   "master_all"    the global list of 11, exactly as Step 4 flagged them, so a
#                   gene counts in every condition (DEFAULT - matches the
#                   pipeline's own is_hub columns)
#   "per_condition" only that condition's declared panel from hubs_for(cond),
#                   e.g. ciprofloxacin = pgrr, rcsc, fadr
# Override from the command line: python coverage_leaves.py <results_dir> per_condition
HUB_SCOPE = "master_all"

DISPLAY = {"amikacin": "Amikacin", "ciprofloxacin": "Ciprofloxacin",
           "doripenem": "Doripenem", "meropenem": "Meropenem",
           "polymyxin": "Polymyxin", "tigecycline": "Tigecycline",
           "colistin_res": "Colistin resistance",
           "colistin resistance": "Colistin resistance"}


def norm(s: str) -> str:
    k = "".join(ch for ch in str(s).lower() if ch.isalpha())
    return "colistinres" if k.startswith("colis") else k


def leaf_sets(d: pd.DataFrame, cols: list, non_reg: set, panel=None) -> tuple:
    """(all-significant leaves, hub-motif leaves, leaf -> hub genes).

    The hub names come from the table's own *_is_hub columns rather than a
    second copy of the Master Hub list, so this can never drift out of step
    with the pipeline. A leaf can be reached by more than one hub motif, so
    the mapping holds a SET of hub genes per leaf.
    """
    if len(d) == 0:
        return set(), set(), {}
    s1 = d[d["passes_step1_stats"].astype(bool)]
    s2 = d[d["passes_significance"].astype(bool)]
    pick = lambda t: {g for c in cols for g in t[c] if g in non_reg}
    node_cols = [c[:-len("_is_hub")] for c in d.columns if c.endswith("_is_hub")]
    hubs_of = {}
    rows = []
    for _, r in s2.iterrows():
        if panel is None:
            row_hubs = {r[c] for c in node_cols if bool(r[f"{c}_is_hub"])}
        else:
            # Re-derive from the declared panel. Step 1 is untouched - only the
            # hub test narrows - so a Step 1 motif whose only hub is foreign to
            # this condition drops out of the numerator.
            row_hubs = {r[c] for c in node_cols if r[c] in panel}
        if not row_hubs:
            continue
        rows.append(r)
        for c in cols:
            if r[c] in non_reg:
                hubs_of.setdefault(r[c], set()).update(row_hubs)
    hub_leaves = {g for r in rows for c in cols if r[c] in non_reg
                  for g in (r[c],)}
    return pick(s1), hub_leaves, hubs_of


def run_scope(scope: str, src: str, graphs: dict, ffl, bif):
    """Build the coverage table for one hub scope."""
    rows, members = [], []
    for cond, G in graphs.items():
        key = norm(cond)
        # A leaf is a gene that regulates nothing in THIS condition's network,
        # judged on the whole GRN - not merely absent from the motif tables.
        non_reg = set(G) - regulators(G)
        panel = (None if scope == "master_all"
                 else {h for h in hubs_for(cond) if h in G})
        f_all, f_hub, f_hubs_of = leaf_sets(
            ffl[ffl["_c"] == key], ["Z"], non_reg, panel)
        b_all, b_hub, b_hubs_of = leaf_sets(
            bif[bif["_c"] == key], ["Z1", "Z2"], non_reg, panel)
        n_leaves, n_hub_leaves = len(f_all | b_all), len(f_hub | b_hub)
        rows.append({
            "condition": DISPLAY.get(cond.lower(), cond),
            "hub_scope": scope,
            "hubs_counted": ";".join(sorted(panel)) if panel else "all_11_master",
            "n_genes": G.number_of_nodes(),
            "n_regulators": len(regulators(G)),
            "n_non_regulators": len(non_reg),
            "ffl_leaves": len(f_all), "ffl_hub_leaves": len(f_hub),
            "bifan_leaves": len(b_all), "bifan_hub_leaves": len(b_hub),
            "n_leaves_sum": len(f_all) + len(b_all),
            "n_hub_leaves_sum": len(f_hub) + len(b_hub),
            "n_leaves": n_leaves, "n_hub_leaves": n_hub_leaves,
            "coverage": round(n_hub_leaves / n_leaves, 4) if n_leaves else None,
            "coverage_pct": (f"{100 * n_hub_leaves / n_leaves:.1f}% "
                             f"({n_hub_leaves}/{n_leaves})") if n_leaves
                            else "n/a"})
        for g in sorted(f_all | b_all):
            hf = sorted(f_hubs_of.get(g, ()))
            hb = sorted(b_hubs_of.get(g, ()))
            members.append({"condition": DISPLAY.get(cond.lower(), cond),
                            "hub_scope": scope, "leaf_gene": g,
                            "in_ffl": g in f_all, "in_bifan": g in b_all,
                            "in_hub_ffl": g in f_hub,
                            "in_hub_bifan": g in b_hub,
                            "reached_by_hub_motif": g in (f_hub | b_hub),
                            "hub_genes_ffl": ";".join(hf),
                            "hub_genes_bifan": ";".join(hb),
                            "hub_genes_all": ";".join(sorted(set(hf) | set(hb))),
                            "n_hub_genes": len(set(hf) | set(hb))})

    t = pd.DataFrame(rows)
    num = [c for c in t.columns if c.startswith(("ffl_", "bifan_", "n_"))]
    total = {c: int(t[c].sum()) for c in num}
    total["condition"] = "Total (instances)"
    total["hub_scope"] = scope
    total["hubs_counted"] = ""
    total["coverage"] = round(total["n_hub_leaves"] / total["n_leaves"], 4)
    total["coverage_pct"] = (f"{100 * total['n_hub_leaves'] / total['n_leaves']:.1f}% "
                             f"({total['n_hub_leaves']}/{total['n_leaves']})")
    t = pd.concat([t, pd.DataFrame([total])[t.columns]], ignore_index=True)

    suffix = "" if scope == "master_all" else "_per_condition"
    write(t, f"coverage_by_condition{suffix}.csv", "coverage")
    write(pd.DataFrame(members), f"coverage_leaf_members{suffix}.csv", "coverage")
    print(f"\n  --- hub scope: {scope} ---")
    print(t[["condition", "ffl_leaves", "ffl_hub_leaves", "bifan_leaves",
             "bifan_hub_leaves", "n_leaves", "n_hub_leaves",
             "coverage_pct"]].to_string(index=False))
    return t


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
        "results", "step4_motif_enrichment")
    # No scope argument: do BOTH, so the global-hub and declared-panel numbers
    # are always produced from the same inputs in one run and cannot be mixed
    # up between sessions. Pass one explicitly to build only that table.
    if len(sys.argv) > 2:
        scopes = [sys.argv[2].lower()]
        if scopes[0] not in ("master_all", "per_condition"):
            raise SystemExit("hub scope must be 'master_all' or "
                             f"'per_condition', got {scopes[0]!r}")
    else:
        scopes = ["master_all", "per_condition"]

    ffl = pd.read_csv(os.path.join(src, "ffl_triad_significance.csv"))
    bif = pd.read_csv(os.path.join(src, "bifan_pair_significance.csv"))
    ffl["_c"] = ffl["condition"].map(norm)
    bif["_c"] = bif["condition"].map(norm)
    graphs = load_all()

    tables = {sc: run_scope(sc, src, graphs, ffl, bif) for sc in scopes}

    if len(tables) == 2:
        a, b = tables["master_all"], tables["per_condition"]
        cmp = a[["condition", "n_leaves", "n_hub_leaves", "coverage"]].merge(
            b[["condition", "n_hub_leaves", "coverage"]], on="condition",
            suffixes=("_master_all", "_per_condition"))
        cmp["leaves_lost_by_narrowing"] = (cmp["n_hub_leaves_master_all"]
                                          - cmp["n_hub_leaves_per_condition"])
        write(cmp, "coverage_scope_comparison.csv", "coverage")
        print("\n  --- both scopes side by side ---")
        print(cmp.to_string(index=False))
    print(f"\n  wrote results to {outdir('coverage')}")


if __name__ == "__main__":
    main()
