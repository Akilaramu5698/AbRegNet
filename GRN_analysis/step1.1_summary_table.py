"""
summary_table.py -- the per-antibiotic summary table.

An FFL is a triad X -> Y, Y -> Z, X -> Z: the regulator X reaches Z both
directly and through Y. 

Bi-fans: X1, X2 -> Z1, Z2 (two regulators that both regulate the same two
targets). Two kinds are counted:
  n_bifan_total  every X1,X2 -> Z1,Z2 instance (no restriction)
  n_bifan_pure   PURE bi-fans: X1 and X2 do not regulate each other, Z1 and Z2
                 do not regulate each other, and no Z regulates an X - the four
                 nodes carry exactly the four X -> Z edges and nothing else.
"Leaf" versions require both Z genes to be leaves (out-degree 0). 
Run:  python summary_table.py      (edge files + grn_common.py beside it)
Out:  results/summary_table.csv          one row per antibiotic (the headline
                                          table: nodes, edges, self-loops, DAG,
                                          cycles, depth, FFL counts, hub chains)
      results/motif_counts.csv           FFL and pure bi-fan counts per condition
                                          (instances and unique leaf genes),
                                          with a Total (instances) row
      results/hub_ffl_details.csv        every FFL involving at least one hub of
                                          that condition: X, Y, Z, which
                                          positions are hubs, FFL class
      results/hub_ffl_gene_summary.csv   per condition x hub: how many FFLs it
                                          takes part in and the gene list
"""
import csv
import os
from itertools import islice

import networkx as nx
import pandas as pd

from grn_common import GRN_FILES, HUBS_BY_CONDITION, hubs_for

ALL_HUBS = sorted({h for hs in HUBS_BY_CONDITION.values() for h in hs})
OUT = os.path.join("results", "summary_table.csv")

# An FFL is "hub-rooted" when this position holds one of the condition's hubs:
#   "X"   the top regulator X          "XY"  X or Y          "any"  X, Y or Z
FFL_HUB_POS = "X"


def load(path):
    with open(path, encoding="utf-8") as f:
        raw = [(r["source"].strip().lower(), r["target"].strip().lower())
               for r in csv.DictReader(f)]
    G = nx.DiGraph()
    G.add_edges_from((u, v) for u, v in raw if u != v)
    return G, sum(1 for u, v in raw if u == v)


def depth_chain(G):
    """Max depth in edges, one example chain, and how many chains tie for it."""
    L, C = {}, {}
    for n in nx.topological_sort(G):
        preds = list(G.predecessors(n))
        if not preds:
            L[n], C[n] = 0, 1
        else:
            L[n] = 1 + max(L[q] for q in preds)
            C[n] = sum(C[q] for q in preds if L[q] == L[n] - 1)
    mx = max(L.values())
    ties = sum(C[n] for n in L if L[n] == mx)
    ends = sorted(n for n in L if L[n] == mx)
    return mx, " -> ".join(nx.dag_longest_path(G)), ties, ends


def all_ffls(G):
    """Every FFL as a (X, Y, Z) tuple, with X -> Y -> Z and the shortcut X -> Z."""
    succ = {n: set(G.successors(n)) for n in G}
    out = []
    for x, y in G.edges():
        for z in succ[y] & succ[x]:
            if z not in (x, y):
                out.append((x, y, z))
    return out


def count_ffls(G):
    """All feed-forward loops X -> Y -> Z with the shortcut X -> Z.

    For each edge X -> Y, any common successor of X and Y closes an FFL, so
    the scan is one set intersection per edge rather than a triple loop.
    """
    tfs = {n for n in G if G.out_degree(n) > 0}
    succ = {n: set(G.successors(n)) for n in G}
    total = tf_only = 0
    for x, y in G.edges():
        for z in succ[y] & succ[x]:
            if z in (x, y):
                continue
            total += 1
            if z in tfs:
                tf_only += 1
    return total, tf_only, total - tf_only


def ffl_leaf_sets(G, hubs):
    """Unique leaf genes (Z with out-degree 0) in all FFLs, and in hub-rooted FFLs."""
    hubset = set(hubs)
    leaves, hub_leaves, n_hub_rooted = set(), set(), 0
    for x, y, z in all_ffls(G):
        rooted = {"X": x in hubset,
                  "XY": x in hubset or y in hubset,
                  "any": bool({x, y, z} & hubset)}[FFL_HUB_POS]
        n_hub_rooted += rooted
        if G.out_degree(z) == 0:
            leaves.add(z)
            if rooted:
                hub_leaves.add(z)
    return n_hub_rooted, leaves, hub_leaves


def bifan_counts(G, hubs):
    """Bi-fan counts X1,X2 -> Z1,Z2 (see module docstring for definitions).

    For every unordered pair of regulators the common targets S are found with
    one set intersection; any 2 of them close a bi-fan, so C(|S|, 2) instances.
    For PURE bi-fans X1/X2 must not touch, targets that regulate X1 or X2 are
    dropped, and Z pairs joined by an edge (either direction) are subtracted.
    """
    hubset = set(hubs)
    succ = {n: set(G.successors(n)) for n in G}
    pred = {n: set(G.predecessors(n)) for n in G}
    regs = sorted(n for n in G if succ[n])
    c = dict(total=0, pure=0, pure_leaf=0, pure_hub_rooted=0, pure_both_hubs=0)
    leaves, hub_leaves = set(), set()
    for i, x1 in enumerate(regs):
        for x2 in regs[i + 1:]:
            S = (succ[x1] & succ[x2]) - {x1, x2}
            if len(S) < 2:
                continue
            c["total"] += len(S) * (len(S) - 1) // 2
            if x2 in succ[x1] or x1 in succ[x2]:
                continue                                  # X1 and X2 touch -> not pure
            S = {z for z in S if z not in pred[x1] and z not in pred[x2]}
            if len(S) < 2:
                continue
            adj = set()                                   # Z pairs joined by an edge
            for z in S:
                for w in succ[z] & S:
                    adj.add(frozenset((z, w)))
            n_pure = len(S) * (len(S) - 1) // 2 - len(adj)
            if n_pure == 0:
                continue
            L = [z for z in S if not succ[z]]             # leaves are never adjacent
            n_leaf = len(L) * (len(L) - 1) // 2
            rooted = x1 in hubset or x2 in hubset
            c["pure"] += n_pure
            c["pure_leaf"] += n_leaf
            if rooted:
                c["pure_hub_rooted"] += n_pure
            if x1 in hubset and x2 in hubset:
                c["pure_both_hubs"] += n_pure
            if n_leaf:
                leaves |= set(L)
                if rooted:
                    hub_leaves |= set(L)
    return c, leaves, hub_leaves


def hub_ffls(G, hubs):
    """FFL triads whose three nodes are all hubs (hub -> hub edges only)."""
    H = G.subgraph([h for h in hubs if h in G])
    out = []
    for x, y in H.edges():
        for z in set(H.successors(y)) & set(H.successors(x)):
            if z not in (x, y):
                out.append(f"{x} -> {y} -> {z} (+ {x} -> {z})")
    return len(out), " / ".join(sorted(out))


def longest_hub_chains(G, hubs):
    """Longest directed chains using hub -> hub edges only.

    The induced subgraph on the hubs is tiny, so every simple path is
    enumerated exactly; all chains of maximum length are returned.
    """
    H = G.subgraph([h for h in hubs if h in G])
    if H.number_of_edges() == 0:
        present = [h for h in hubs if h in G]
        return 0, f"(no hub-to-hub edges; {', '.join(present)} act in parallel)"
    best, chains = 0, []
    for s in H:
        for t in H:
            if s == t:
                continue
            for p in nx.all_simple_paths(H, s, t):
                if len(p) - 1 > best:
                    best, chains = len(p) - 1, [p]
                elif len(p) - 1 == best:
                    chains.append(p)
    uniq, seen = [], set()
    for p in chains:                       # drop sub-paths of longer chains
        if tuple(p) not in seen:
            seen.add(tuple(p))
            uniq.append(" -> ".join(p))
    return best, " / ".join(sorted(uniq))


rows, ffl_rows, hub_rows, motif_rows = [], [], [], []
for cond, fn in GRN_FILES.items():
    if not os.path.exists(fn):
        print(f"  skipped (not found): {fn}")
        continue
    G, loops = load(fn)
    is_dag = nx.is_directed_acyclic_graph(G)
    n_cyc = 0 if is_dag else len(list(islice(nx.simple_cycles(G), 100000)))
    n_own, chain_own = longest_hub_chains(G, hubs_for(cond))
    n_all, chain_all = longest_hub_chains(G, ALL_HUBS)
    ffl_total, ffl_tf, ffl_leaf = count_ffls(G)
    if is_dag:
        mx, chain, ties, ends = depth_chain(G)
    else:
        mx, chain, ties, ends = None, "", None, []
    n_ffl_hub_rooted, ffl_leaves, ffl_hub_leaves = ffl_leaf_sets(G, hubs_for(cond))
    bf, bf_leaves, bf_hub_leaves = bifan_counts(G, hubs_for(cond))
    all_leaves = ffl_leaves | bf_leaves
    all_hub_leaves = ffl_hub_leaves | bf_hub_leaves
    motif_rows.append({
        "antibiotic_condition": cond,
        "n_ffl": ffl_total, "n_ffl_hub_rooted": n_ffl_hub_rooted,
        "n_ffl_leaves": len(ffl_leaves), "n_ffl_hub_leaves": len(ffl_hub_leaves),
        "n_bifan_total": bf["total"], "n_bifan_pure": bf["pure"],
        "n_bifan_pure_leaf": bf["pure_leaf"],
        "n_bifan_pure_hub_rooted": bf["pure_hub_rooted"],
        "n_bifan_pure_both_hubs": bf["pure_both_hubs"],
        "n_bifan_leaves": len(bf_leaves), "n_bifan_hub_leaves": len(bf_hub_leaves),
        "n_leaves_ffl_or_bifan": len(all_leaves),
        "n_hub_leaves_ffl_or_bifan": len(all_hub_leaves),
        "hub_leaf_coverage": (f"{100.0 * len(all_hub_leaves) / len(all_leaves):.1f}% "
                              f"({len(all_hub_leaves)}/{len(all_leaves)})"
                              if all_leaves else "n/a"),
    })
    n_ffl_own, ffl_own = hub_ffls(G, hubs_for(cond))
    n_ffl_all, ffl_all = hub_ffls(G, ALL_HUBS)
    # ---- FFLs involving this condition's hubs, gene by gene
    hubs_own = [h for h in hubs_for(cond) if h in G]
    triads = all_ffls(G)
    tfs = {n for n in G if G.out_degree(n) > 0}
    involved = {h: set() for h in hubs_own}
    n_ffl_hub = {h: 0 for h in hubs_own}
    for x, y, z in triads:
        hit = [h for h in hubs_own if h in (x, y, z)]
        if not hit:
            continue
        for h in hit:
            n_ffl_hub[h] += 1
            involved[h] |= {x, y, z} - {h}
        ffl_rows.append({
            "antibiotic_condition": cond, "X": x, "Y": y, "Z": z,
            "hubs_in_triad": ";".join(hit),
            "X_is_hub": x in hubs_own, "Y_is_hub": y in hubs_own,
            "Z_is_hub": z in hubs_own,
            "all_three_are_hubs": len(hit) == 3,
            "ffl_class": "TF_only" if z in tfs else "TF_target",
            "genes": f"{x};{y};{z}",
        })
    for h in hubs_for(cond):
        present = h in G
        hub_rows.append({
            "antibiotic_condition": cond, "hub": h, "present_in_network": present,
            "k_in": G.in_degree(h) if present else 0,
            "k_out": G.out_degree(h) if present else 0,
            "n_ffls_involving_hub": n_ffl_hub.get(h, 0),
            "n_partner_genes_in_ffls": len(involved.get(h, ())),
            "partner_genes_in_ffls": ";".join(sorted(involved.get(h, ()))),
            "network_is_DAG": is_dag,
            "network_nodes": G.number_of_nodes(),
            "network_directed_edges": G.number_of_edges(),
        })

    rows.append({
        "antibiotic_condition": cond,
        "file": fn,
        "nodes": G.number_of_nodes(),
        "directed_edges": G.number_of_edges(),
        "self_loops": loops,
        "is_DAG": is_dag,
        "n_simple_cycles": n_cyc,
        "max_dag_depth_edges": mx,
        "max_dag_depth_genes": mx + 1 if is_dag else None,
        "max_depth_chain_example": chain,
        "n_distinct_max_depth_chains": ties,
        "max_depth_endpoints": ";".join(ends),
        "condition_hub_panel": ";".join(hubs_for(cond)),
        "hub_chain_len_condition_hubs": n_own,
        "chain_condition_hubs": chain_own,
        "hub_chain_len_all_hubs": n_all,
        "chain_all_hubs": chain_all,
        "n_ffl_total": ffl_total,
        "n_ffl_tf_only": ffl_tf,
        "n_ffl_tf_target": ffl_leaf,
        "n_ffl_hub_rooted": n_ffl_hub_rooted,
        "n_bifan_total": bf["total"],
        "n_bifan_pure": bf["pure"],
        "n_bifan_pure_leaf": bf["pure_leaf"],
        "n_bifan_pure_hub_rooted": bf["pure_hub_rooted"],
        "n_bifan_pure_both_hubs": bf["pure_both_hubs"],
        "n_ffl_among_condition_hubs": n_ffl_own,
        "ffl_condition_hubs": ffl_own,
        "n_ffl_among_all_hubs": n_ffl_all,
        "ffl_all_hubs": ffl_all,
    })

df = pd.DataFrame(rows)
os.makedirs("results", exist_ok=True)
df.to_csv(OUT, index=False)
pd.DataFrame(ffl_rows).to_csv(os.path.join("results", "hub_ffl_details.csv"),
                              index=False)
motif_df = pd.DataFrame(motif_rows)
if len(motif_df):
    total = {"antibiotic_condition": "Total (instances)"}
    for col in motif_df.columns[1:-1]:
        total[col] = int(motif_df[col].sum())
    total["hub_leaf_coverage"] = (
        f"{100.0 * total['n_hub_leaves_ffl_or_bifan'] / total['n_leaves_ffl_or_bifan']:.1f}% "
        f"({total['n_hub_leaves_ffl_or_bifan']}/{total['n_leaves_ffl_or_bifan']})"
        if total["n_leaves_ffl_or_bifan"] else "n/a")
    motif_df = pd.concat([motif_df, pd.DataFrame([total])], ignore_index=True)
motif_df.to_csv(os.path.join("results", "motif_counts.csv"), index=False)
hub_df = pd.DataFrame(hub_rows)
hub_df.to_csv(os.path.join("results", "hub_ffl_gene_summary.csv"), index=False)
print(hub_df[["antibiotic_condition", "hub", "present_in_network", "k_out",
              "n_ffls_involving_hub", "n_partner_genes_in_ffls",
              "network_is_DAG", "network_nodes",
              "network_directed_edges"]].to_string(index=False))
print()
print(df[["antibiotic_condition", "nodes", "directed_edges", "self_loops",
          "is_DAG", "n_simple_cycles", "max_dag_depth_edges",
          "n_ffl_total", "n_ffl_tf_only", "n_bifan_total", "n_bifan_pure",
          "n_ffl_among_all_hubs",
          "max_depth_chain_example", "n_distinct_max_depth_chains"]]
      .to_string(index=False))
print("\nMotif counts (structural, whole network):")
print(motif_df.to_string(index=False))
print(f"\nwrote {OUT}, results/motif_counts.csv, results/hub_ffl_details.csv, "
      f"results/hub_ffl_gene_summary.csv")
