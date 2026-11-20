"""
step4_motif_enrichment.py -- Step 4: statistical validation of FFL and bi-fan
abundance against a degree-preserving null model.

Null models (NULL_MODEL at the top of the file)
-----------------------------------------------
"swap" (default)  Double-edge swap: (a -> b, c -> d) becomes (a -> d, c -> b),
                  rejecting self-loops and duplicates. Every node keeps its
                  EXACT k_in and k_out and the graph stays simple, so a hub
                  with k_out = 136 still has 136 out-edges in every replicate.

"configuration"   nx.directed_configuration_model, then collapsed to a simple
                  DiGraph with self-loops removed. Faster, but the collapse
                  deletes parallel edges, which arise preferentially between
                  high-degree nodes - so hubs lose out-edges and the null
                  becomes weaker than the real network. Measured on one of
                  these GRNs this cost 2-4% of edges, concentrated on the
                  hubs, which INFLATES apparent enrichment. If you use it,
                  the realised degree error is reported in
                  null_model_quality.csv - quote it alongside the p-values.

Statistics
----------
Network level, per motif class:
    observed, null mean +/- sd, fold = observed / null mean,
    z = (observed - mean) / sd,
    p = (#{null >= observed} + 1) / (N + 1)        [one-sided upper tail]

Per individual motif:
    ffl_triad_significance.csv   one row per observed X -> Y -> Z triad: in how
                                 many replicates the same triad reappears, with
                                 p, a Bernoulli Z_score and a single
                                 q_value_BH set by Q_SCOPE, plus the two-stage
                                 flags passes_step1_stats (q < 0.05 and
                                 Z >= 2.0), contains_master_hub and
                                 passes_significance (both)
    bifan_pair_significance.csv  the same columns for every observed pure
                                 (mfinder ID204) bi-fan square X1,X2 -> Z1,Z2.
                                 One row per SQUARE, not per regulator pair:
                                 a pair with s shared targets gives C(s,2) rows
    ffl_pair_significance.csv    each observed (X -> Y) seed pair: how many
                                 FFLs it seeds, against the null distribution
                                 of that same pair's FFL count. This table
                                 keeps the older q_value_BH_global /
                                 q_value_BH_within_condition / significant_q05
                                 / n_tests_in_group columns

Two-stage filter (flagged AND written out as separate files)
------------------------------------------------------------
Step 1  BH q < 0.05 (within condition; FFLs and pure bi-fans are separate
        test families) AND Z >= 2.0            -> passes_step1_stats
Step 2  Step 1 survivors that contain >= 1 of the 11 MASTER_HUBS
                                               -> passes_significance
Every motif row is kept in the main CSV with these flags, and the survivors
are ALSO written to *_step1_pass.csv and *_step2_pass.csv, plus
step1_step2_summary.csv with the counts.


Outputs (results/step4_motif_enrichment/)
  motif_enrichment_per_network.csv, hub_motif_enrichment.csv,
  ffl_triad_significance.csv, ffl_pair_significance.csv,
  bifan_pair_significance.csv, All_BiFans_Total.csv,
  null_distributions.csv, null_model_quality.csv, step4_report.md,
  ffl_triad_significance_step1_pass.csv / _step2_pass.csv,
  bifan_pair_significance_step1_pass.csv / _step2_pass.csv,
  step1_step2_summary.csv

Run:  python step4_motif_enrichment.py
"""
from __future__ import annotations

import os
import random
from collections import defaultdict, Counter
from itertools import combinations
from math import comb, sqrt

import numpy as np
import pandas as pd
import networkx as nx

from grn_common import load_all, regulators, outdir, write, hubs_for

# Which null model to use. Default "swap"; override from the command line:
#     python step4_motif_enrichment.py swap
#     python step4_motif_enrichment.py configuration
# Each model writes to its own results folder, so both runs can coexist.
import sys as _sys
NULL_MODEL = (_sys.argv[1].lower() if len(_sys.argv) > 1 else
              os.environ.get("GRN_NULL_MODEL", "swap"))
if NULL_MODEL not in ("swap", "configuration"):
    raise SystemExit(f"NULL_MODEL must be 'swap' or 'configuration', got {NULL_MODEL!r}")
SUB = "step4_motif_enrichment" + ("" if NULL_MODEL == "swap" else "_configuration")
N_RANDOMISATIONS = 1000
SWAP_FACTOR = 10             # SUCCESSFUL swaps per replicate = SWAP_FACTOR * |E|
MAX_ATTEMPT_FACTOR = 1000    # safety cap: attempts <= this x the swap target
SEED = 42
PER_MOTIF_MIN_OBS = 1
MASTER_HUBS = frozenset(["arsr2", "beti", "dmlr", "fadr", "lrp", "paax",
                         "pgrr", "rcsc", "rhas", "rpoh", "yybr"])
# The 11 Master Hubs are a GLOBAL list: a gene is flagged in every condition,
# not only in the conditions where it was nominated. The per-condition panels
# (hubs_for) still drive hub_motif_enrichment.csv, which asks a different
# question - whether a given hub is enriched in its own condition.
Q_SCOPE = "within_condition"   # "within_condition" | "global" - the test family
                               # the workbook's single q_value_BH is corrected over


def adjacency(G: nx.DiGraph) -> dict:
    return {n: set(G.successors(n)) - {n} for n in G}


def ffl_triads(succ: dict) -> list:
    """Every FFL as (X, Y, Z): X -> Y, Y -> Z, X -> Z."""
    out = []
    for x, tx in succ.items():
        for y in tx:
            for z in succ[y] & tx:
                if z != x:
                    out.append((x, y, z))
    return out


def bifan_pairs(succ: dict) -> Counter:
    """{(X1, X2): n_shared_targets} for every co-regulating pair."""
    return Counter({k: len(v) for k, v in bifan_shared_sets(succ).items()})


def bifan_shared_sets(succ: dict) -> dict:
    """{(X1, X2): set of co-regulated targets} - needed for per-square work."""
    by_target = defaultdict(list)
    for t, tg in succ.items():
        for z in tg:
            by_target[z].append(t)
    shared = defaultdict(set)
    for z, regs in by_target.items():
        if len(regs) > 1:
            for a, b in combinations(sorted(regs), 2):
                shared[(a, b)].add(z)
    return shared


def is_pure_bifan(G: nx.DiGraph, x1, x2, z1, z2) -> bool:
    """True when the induced subgraph on the four nodes is EXACTLY the bi-fan.

    The canonical 4-node bi-fan (mfinder motif ID 204) has precisely the four
    edges X1->Z1, X1->Z2, X2->Z1, X2->Z2. If the induced subgraph carries any
    further edge - X1->X2, Z1->Z2, a back-edge, a reciprocal - the instance is
    a different 4-node pattern that merely contains a bi-fan. "Pure" marks the
    induced-subgraph definition; without it the count is the non-induced one.
    """
    want = {(x1, z1), (x1, z2), (x2, z1), (x2, z2)}
    nodes = [x1, x2, z1, z2]
    got = {(u, v) for u in nodes for v in nodes if u != v and G.has_edge(u, v)}
    return got == want


def is_pure_bifan_succ(succ: dict, x1, x2, z1, z2) -> bool:
    """is_pure_bifan() on an adjacency dict (used on every null replicate).

    Pure = the induced subgraph on the four nodes is exactly the four
    X->Z edges: no X1<->X2, no Z1<->Z2, no Z->X edge. The four X->Z edges are
    assumed present (callers pass squares taken from shared-target sets).
    """
    sx1, sx2 = succ.get(x1, ()), succ.get(x2, ())
    sz1, sz2 = succ.get(z1, ()), succ.get(z2, ())
    return not (x2 in sx1 or x1 in sx2 or z2 in sz1 or z1 in sz2
                or x1 in sz1 or x2 in sz1 or x1 in sz2 or x2 in sz2)


def count_pure_bifans(succ: dict) -> int:
    """Number of pure (induced, ID 204) bi-fan squares in a network.

    Per co-regulating pair (X1, X2) with no edge between them, keep the shared
    targets that do not point back at X1/X2; a pair of those is a pure bi-fan
    unless the two targets are linked. So the count is
    C(k, 2) - (number of linked target pairs), with no square enumeration.
    """
    total = 0
    for (x1, x2), zs in bifan_shared_sets(succ).items():
        if len(zs) < 2 or x2 in succ[x1] or x1 in succ[x2]:
            continue
        ok = {z for z in zs if x1 not in succ[z] and x2 not in succ[z]}
        k = len(ok)
        if k < 2:
            continue
        linked = set()
        for z in ok:
            for w in succ[z] & ok:
                linked.add((z, w) if z < w else (w, z))
        total += comb(k, 2) - len(linked)
    return total


def summarise(succ: dict, tfs: set):
    tri = ffl_triads(succ)
    pairs = bifan_pairs(succ)   # counts; bifan_shared_sets() gives the members
    n_tf_only = sum(1 for _, _, z in tri if z in tfs)
    motifs = sum(comb(s, 2) for s in pairs.values() if s >= 2)
    n_pairs = sum(1 for s in pairs.values() if s >= 2)
    return tri, pairs, len(tri), n_tf_only, motifs, n_pairs


def null_swap(G, rng, swap_factor):
    """Maslov-Sneppen double-edge swap with an exact SUCCESSFUL-swap budget.

    The loop runs until swap_factor x E_eff swaps have actually been applied
    (E_eff = number of non-self-loop edges), so rejected proposals do not eat
    into the budget. A proposal (a->b, c->d) -> (a->d, c->b) is rejected when
    it would create a self-loop or a parallel edge, or when the two edges share
    a node (then the swap is a no-op or creates a self-loop); every ACCEPTED
    swap therefore changes the graph. Self-loops (if any) stay untouched, so
    k_in and k_out are preserved exactly.

    Returns (R, info) with info = {swaps_target, swaps_done, swap_attempts}.
    Raises RuntimeError if the target cannot be reached within
    MAX_ATTEMPT_FACTOR x target attempts (graph too constrained to rewire).
    """
    R = G.copy()
    edges = [(u, v) for u, v in R.edges() if u != v]
    m = len(edges)
    target = swap_factor * m
    done = attempts = 0
    if m >= 2:
        max_attempts = MAX_ATTEMPT_FACTOR * target
        while done < target:
            if attempts >= max_attempts:
                raise RuntimeError(
                    f"null_swap: only {done}/{target} successful swaps after "
                    f"{attempts} attempts - graph too constrained")
            attempts += 1
            i, j = rng.randrange(m), rng.randrange(m)
            (a, b), (c, d) = edges[i], edges[j]
            if len({a, b, c, d}) < 4 or R.has_edge(a, d) or R.has_edge(c, b):
                continue
            R.remove_edge(a, b)
            R.remove_edge(c, d)
            R.add_edge(a, d)
            R.add_edge(c, b)
            edges[i], edges[j] = (a, d), (c, b)
            done += 1
    return R, {"swaps_target": target, "swaps_done": done,
               "swap_attempts": attempts}


def null_configuration(G, seed_val):
    """Stub-matching null, collapsed to a simple DiGraph (degree-lossy)."""
    nodes = list(G)
    in_deg = [G.in_degree(n) for n in nodes]
    out_deg = [G.out_degree(n) for n in nodes]
    R_multi = nx.directed_configuration_model(in_deg, out_deg, seed=seed_val)
    R = nx.DiGraph(R_multi)                       # drops parallel edges
    R.remove_edges_from(list(nx.selfloop_edges(R)))
    return nx.relabel_nodes(R, dict(enumerate(nodes)))


def make_null(G, rng, replicate):
    """Return (R, swap_info); swap_info is {} for the configuration model."""
    if NULL_MODEL == "configuration":
        return null_configuration(G, SEED + replicate), {}
    return null_swap(G, rng, SWAP_FACTOR)


def two_step_flags(d, hubcol):
    """Step 1 = statistical filter, Step 2 = Master Hub filter.

    Kept as three separate columns so the two stages can be counted
    independently: a motif that fails Step 1 is not evidence whatever hubs it
    contains, and a motif that passes Step 1 but holds no Master Hub is still
    a real enrichment, just not a hub-driven one.
    """
    # Z >= 2.0 is redundant in practice (q < 0.05 already implies Z of about 4.4+
    # with N = 1000) but is kept so the stated criterion is applied literally.
    z = pd.to_numeric(d["Z_score"], errors="coerce")
    d["passes_step1_stats"] = (d["q_value_BH"] < 0.05) & (z >= 2.0)
    d["contains_master_hub"] = d[hubcol] >= 1
    d["passes_significance"] = d["passes_step1_stats"] & d["contains_master_hub"]
    return d


def write_two_step(d, stem, keycols, family, summary):
    """Write the full flagged table plus the Step 1 and Step 2 survivor files."""
    d = d.sort_values(["condition"] + keycols)
    write(d, f"{stem}.csv", SUB)
    s1 = d[d["passes_step1_stats"]]
    s2 = d[d["passes_significance"]]
    write(s1, f"{stem}_step1_pass.csv", SUB)
    write(s2, f"{stem}_step2_pass.csv", SUB)
    for cond, g in d.groupby("condition"):
        summary.append({"motif_family": family, "condition": cond,
                        "n_tested": len(g),
                        "n_step1_pass": int(g["passes_step1_stats"].sum()),
                        "n_step2_pass_with_master_hub":
                            int(g["passes_significance"].sum()),
                        "n_step1_pass_no_master_hub":
                            int((g["passes_step1_stats"]
                                 & ~g["contains_master_hub"]).sum())})
    print(f"  {family}: {len(d)} tested, {len(s1)} pass Step 1, "
          f"{len(s2)} also contain a Master Hub (Step 2)")


def degree_error(G, R) -> dict:
    """Per-replicate fidelity of one null network against the real one.

    Both degree sequences are checked. The swap model preserves k_in and
    k_out exactly by construction; the configuration model can lose edges
    when it proposes a self-loop or a duplicate, and that loss shows up in
    BOTH directions - a dropped edge costs its source an out-edge and its
    target an in-edge - so reporting only k_out hides half the damage.
    """
    dout = [G.out_degree(n) - (R.out_degree(n) if n in R else 0) for n in G]
    din = [G.in_degree(n) - (R.in_degree(n) if n in R else 0) for n in G]
    return {"edges_real": G.number_of_edges(), "edges_null": R.number_of_edges(),
            "edges_lost": G.number_of_edges() - R.number_of_edges(),
            "n_nodes_kout_changed": int(sum(1 for d in dout if d != 0)),
            "n_nodes_kin_changed": int(sum(1 for d in din if d != 0)),
            "max_kout_loss": int(max(dout)) if dout else 0,
            "max_kin_loss": int(max(din)) if din else 0,
            "mean_abs_kout_error": round(float(np.mean(np.abs(dout))), 4),
            "mean_abs_kin_error": round(float(np.mean(np.abs(din))), 4)}


def summarise_quality(cond: str, records: list) -> dict:
    """Collapse the per-replicate fidelity records into one row per condition.

    Replicate 0 alone is enough for the swap model, where every replicate is
    degree-exact by construction, but not for the configuration model, where
    the number of edges lost varies from draw to draw. Reporting the mean and
    the worst case over all replicates makes that variation visible.
    """
    row = {"condition": cond, "null_model": NULL_MODEL,
           "n_replicates_checked": len(records),
           "edges_real": records[0]["edges_real"]}
    for k in ("edges_lost", "n_nodes_kout_changed", "n_nodes_kin_changed",
              "max_kout_loss", "max_kin_loss",
              "mean_abs_kout_error", "mean_abs_kin_error"):
        v = [r[k] for r in records]
        row[f"{k}_mean"] = round(float(np.mean(v)), 4)
        row[f"{k}_max"] = round(float(np.max(v)), 4)
    if records and "swaps_done" in records[0]:
        row["swaps_target_per_replicate"] = records[0]["swaps_target"]
        row["swaps_done_min"] = int(min(r["swaps_done"] for r in records))
        row["swap_attempts_mean"] = round(float(np.mean(
            [r["swap_attempts"] for r in records])), 1)
        row["swap_acceptance_rate"] = round(
            sum(r["swaps_done"] for r in records)
            / max(1, sum(r["swap_attempts"] for r in records)), 4)
        row["replicates_swap_target_met"] = int(sum(
            1 for r in records if r["swaps_done"] == r["swaps_target"]))
    row["replicates_degree_exact"] = int(sum(1 for r in records if r["edges_lost"] == 0
                                             and r["n_nodes_kout_changed"] == 0
                                             and r["n_nodes_kin_changed"] == 0))
    return row


def stats(obs, null) -> dict:
    null = np.asarray(null, dtype=float)
    mean = float(null.mean())
    sd = float(null.std(ddof=1)) if len(null) > 1 else 0.0
    return {"observed": obs, "null_mean": round(mean, 3), "null_sd": round(sd, 3),
            "null_min": float(null.min()), "null_max": float(null.max()),
            "fold_enrichment": round(obs / mean, 3) if mean else np.inf,
            "z_score": round((obs - mean) / sd, 2) if sd else np.nan,
            # 6 dp everywhere: rounding p before the FDR step changes the tie
            # ordering and makes the same correction return different q-values
            # in different tables
            "p_empirical_upper": round((int((null >= obs).sum()) + 1) / (len(null) + 1), 6)}


def bh_fdr(p) -> np.ndarray:
    """Benjamini-Hochberg q-values (step-up, monotonicity enforced).

    Note on resolution: with N randomisations the smallest attainable
    empirical p is 1/(N+1), so with N = 200 many motifs tie at p = 0.00498.
    BH handles ties correctly, but a q-value cannot be smaller than that
    floor - raise N_RANDOMISATIONS if you need finer discrimination among the
    strongest motifs.
    """
    p = np.asarray(p, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = np.minimum.accumulate((p[order] * n / (np.arange(n) + 1))[::-1])[::-1]
    q = np.empty(n)
    q[order] = np.clip(ranked, 0, 1)
    return q


def add_fdr(df: pd.DataFrame, pcol: str, by: str = "condition") -> pd.DataFrame:
    """Add BH q-values across the whole table and within each condition."""
    if df.empty or pcol not in df.columns:
        return df
    df = df.copy()
    df["q_value_BH_global"] = bh_fdr(df[pcol].values).round(6)
    if by in df.columns:
        df["q_value_BH_within_condition"] = (
            df.groupby(by)[pcol].transform(lambda s: bh_fdr(s.values)).round(6))
        df["significant_q05"] = df["q_value_BH_within_condition"] <= 0.05
    else:
        df["significant_q05"] = df["q_value_BH_global"] <= 0.05
    df["n_tests_in_group"] = (df.groupby(by)[pcol].transform("size")
                              if by in df.columns else len(df))
    return df


def bernoulli_z(n_hits, n_rep) -> float:
    """Z of an observed presence (=1) vs a 0/1 null indicator, sample SD (ddof=1).

    f = 1 -> 0.0 (exactly as common as the null); f = 0 -> inf.
    """
    if n_hits == n_rep:
        return 0.0
    freq = n_hits / n_rep
    var = freq * (1 - freq) * n_rep / (n_rep - 1) if n_rep > 1 else 0.0
    sd = sqrt(var)
    return float("inf") if sd == 0 else (1 - freq) / sd


def bernoulli_stats(n_hits, n_rep) -> dict:
    """For a motif that is simply present or absent in each replicate."""
    freq = n_hits / n_rep
    # freq == 1 (motif present in every replicate) is not undefined: the motif
    # is exactly as common as the null expects, so z = 0. Only freq == 0 has no
    # finite z, and that row is read on p and the hit count instead.
    if n_hits == n_rep:
        return {"n_null_replicates_containing": n_hits, "null_frequency": 1.0,
                "z_score": 0.0, "p_empirical": 1.0}
    return {"n_null_replicates_containing": n_hits,
            "null_frequency": round(freq, 4),
            "z_score": round(bernoulli_z(n_hits, n_rep), 2),
            # 6 dp, matching the workbook tables: rounding p to 4 dp before the
            # FDR step perturbs the tie ordering and makes the same correction
            # return slightly different q-values in the two files
            "p_empirical": round((n_hits + 1) / (n_rep + 1), 6)}


def main():
    graphs = load_all()
    print(f"  null model: {NULL_MODEL}, {N_RANDOMISATIONS} replicates")
    net_rows, hub_rows, raw, qual = [], [], [], []
    wb_ffl_rows, wb_bifan_rows = [], []      # workbook-schema tables
    tri_rows, pair_rows, bif_rows = [], [], []

    for cond, G in graphs.items():
        rng = random.Random(SEED)
        tfs = regulators(G)
        succ = adjacency(G)
        hubs = [h for h in hubs_for(cond) if h in G]

        obs_tri, obs_pairs, n_ffl, n_tf, n_bif, n_cop = summarise(succ, tfs)
        obs_tri_set = set(obs_tri)
        obs_seed = Counter((x, y) for x, y, _ in obs_tri)
        obs_seed_targets = defaultdict(set)
        for x, y, z in obs_tri:
            obs_seed_targets[(x, y)].add(z)
        obs_hub = Counter(h for t in obs_tri for h in set(t) & set(hubs))
        print(f"[step4] {cond}: {n_ffl} FFLs, {n_bif} bi-fans, "
              f"{len(obs_tri_set)} distinct triads")

        # --- workbook-format per-square bi-fan enumeration (observed)
        obs_shared = bifan_shared_sets(succ)
        squares = []                      # (X1, X2, Z1, Z2)
        for (x1, x2), zs in obs_shared.items():
            if len(zs) >= 2:
                for z1, z2 in combinations(sorted(zs), 2):
                    squares.append((x1, x2, z1, z2))
        sq_hits = Counter()
        pure_sq_hits = Counter()          # same square reappears AND is pure
        obs_pure_n = count_pure_bifans(succ)
        null_pure = []
        print(f"          {len(squares)} bi-fan squares to track per replicate")

        qual_reps = []
        null_ffl, null_tf, null_bif, null_cop = [], [], [], []
        null_hub = {h: [] for h in hubs}
        tri_hits = Counter()
        seed_counts = defaultdict(list)
        pair_shared = defaultdict(list)

        for i in range(N_RANDOMISATIONS):
            R, swap_info = make_null(G, rng, i)
            qual_reps.append({**degree_error(G, R), **swap_info})
            rsucc = adjacency(R)
            rtfs = {n for n in R if R.out_degree(n) > 0}
            rtri, rpairs, f, ftf, bf, bp = summarise(rsucc, rtfs)
            null_ffl.append(f); null_tf.append(ftf)
            null_bif.append(bf); null_cop.append(bp)
            n_pure_rep = count_pure_bifans(rsucc)
            null_pure.append(n_pure_rep)

            for t in obs_tri_set & set(rtri):
                tri_hits[t] += 1
            rseed = Counter((x, y) for x, y, _ in rtri)
            for key in obs_seed:
                seed_counts[key].append(rseed.get(key, 0))
            for key in obs_pairs:
                pair_shared[key].append(rpairs.get(key, 0))
            rshared = bifan_shared_sets(rsucc)
            for (x1, x2, z1, z2) in squares:
                zs = rshared.get((x1, x2))
                if zs is not None and z1 in zs and z2 in zs:
                    sq_hits[(x1, x2, z1, z2)] += 1
                    if is_pure_bifan_succ(rsucc, x1, x2, z1, z2):
                        pure_sq_hits[(x1, x2, z1, z2)] += 1
            rhub = Counter(h for t in rtri for h in set(t) & set(hubs))
            for h in hubs:
                null_hub[h].append(rhub.get(h, 0))
            raw.append({"condition": cond, "replicate": i + 1, "n_ffl": f,
                        "n_ffl_tf_only": ftf, "n_bifan_motifs": bf,
                        "n_pure_bifan_motifs": n_pure_rep,
                        "n_coregulating_pairs": bp})

        qual.append(summarise_quality(cond, qual_reps))

        for label, obs, null in [("FFL_total", n_ffl, null_ffl),
                                 ("FFL_TF_only", n_tf, null_tf),
                                 ("bifan_motifs", n_bif, null_bif),
                                 ("pure_bifan_motifs", obs_pure_n, null_pure),
                                 ("coregulating_pairs", n_cop, null_cop)]:
            net_rows.append({"condition": cond, "motif": label,
                             "null_model": NULL_MODEL,
                             "n_nodes": G.number_of_nodes(),
                             "n_edges": G.number_of_edges(),
                             "n_randomisations": N_RANDOMISATIONS,
                             **stats(obs, null)})

        for h in hubs:
            hub_rows.append({"condition": cond, "hub": h, "k_in": G.in_degree(h),
                             "k_out": G.out_degree(h), "motif": "FFL_involving_hub",
                             "null_model": NULL_MODEL,
                             "n_randomisations": N_RANDOMISATIONS,
                             **stats(obs_hub.get(h, 0), null_hub[h])})

        for (x1, x2, z1, z2) in squares:
            hits = sq_hits.get((x1, x2, z1, z2), 0)
            phits = pure_sq_hits.get((x1, x2, z1, z2), 0)
            wb_bifan_rows.append({
                "condition": cond, "X1": x1, "X2": x2, "Z1": z1, "Z2": z2,
                "X1_is_hub": x1 in MASTER_HUBS, "X2_is_hub": x2 in MASTER_HUBS,
                "Z1_is_hub": z1 in MASTER_HUBS, "Z2_is_hub": z2 in MASTER_HUBS,
                "n_hubs_in_bifan": sum(1 for n in (x1, x2, z1, z2)
                                       if n in MASTER_HUBS),
                "is_pure_ID204": is_pure_bifan(G, x1, x2, z1, z2),
                "rand_hits": hits,
                "Z_score": round(bernoulli_z(hits, N_RANDOMISATIONS), 4),
                "p_value": round((hits + 1) / (N_RANDOMISATIONS + 1), 6),
                "rand_hits_pure": phits,
                "Z_score_pure": round(bernoulli_z(phits, N_RANDOMISATIONS), 4),
                "p_value_pure": round((phits + 1) / (N_RANDOMISATIONS + 1), 6)})

        for (x, y, z) in sorted(obs_tri_set):
            tri_rows.append({"condition": cond, "X": x, "Y": y, "Z": z,
                             "triad": f"{x}->{y}->{z}",
                             "ffl_class": "TF_only" if z in tfs else "TF_target",
                             "involves_hub": bool({x, y, z} & set(hubs)),
                             "n_randomisations": N_RANDOMISATIONS,
                             **bernoulli_stats(tri_hits.get((x, y, z), 0),
                                               N_RANDOMISATIONS)})
            hits = tri_hits.get((x, y, z), 0)
            wb_ffl_rows.append({
                "condition": cond, "X": x, "Y": y, "Z": z,
                "X_is_hub": x in MASTER_HUBS, "Y_is_hub": y in MASTER_HUBS,
                "Z_is_hub": z in MASTER_HUBS,
                "n_hubs_in_triad": sum(1 for n in (x, y, z) if n in MASTER_HUBS),
                "rand_hits": hits,
                "Z_score": round(bernoulli_z(hits, N_RANDOMISATIONS), 4),
                "p_value": round((hits + 1) / (N_RANDOMISATIONS + 1), 6)})
        for (x, y), obs_n in obs_seed.items():
            if obs_n >= PER_MOTIF_MIN_OBS:
                st = stats(obs_n, seed_counts[(x, y)])
                st["observed_n_ffls"] = st.pop("observed")
                zs = sorted(obs_seed_targets[(x, y)])
                pair_rows.append({"condition": cond, "X": x, "Y": y,
                                  "involves_hub": bool({x, y} & set(hubs)),
                                  "motif": "FFLs_seeded_by_pair",
                                  "n_target_genes_Z": len(zs),
                                  "target_genes_Z": ";".join(zs), **st})
        for (a, b), obs_s in obs_pairs.items():
            if obs_s >= 2:
                st = stats(obs_s, pair_shared[(a, b)])
                st["observed_shared_targets"] = st.pop("observed")
                tg = sorted(obs_shared.get((a, b), ()))
                bif_rows.append({"condition": cond, "X1": a, "X2": b,
                                 "involves_hub": bool({a, b} & set(hubs)),
                                 "motif": "shared_targets_of_pair",
                                 "observed_bifan_motifs": comb(obs_s, 2),
                                 "shared_target_genes": ";".join(tg), **st})

    # Benjamini-Hochberg FDR on every table that carries a p-value. The
    # per-motif tables test hundreds of motifs per condition, so q - not p -
    # is the number to quote when naming an individual circuit.
    net_df = add_fdr(pd.DataFrame(net_rows), "p_empirical_upper")
    hub_df = add_fdr(pd.DataFrame(hub_rows), "p_empirical_upper")
    tri_df = add_fdr(pd.DataFrame(tri_rows), "p_empirical")
    pair_df = add_fdr(pd.DataFrame(pair_rows), "p_empirical_upper")
    bif_df = add_fdr(pd.DataFrame(bif_rows), "p_empirical_upper")

    write(net_df, "motif_enrichment_per_network.csv", SUB)
    write(hub_df, "hub_motif_enrichment.csv", SUB)
    write(pair_df.sort_values(["condition", "p_empirical_upper"]),
          "ffl_pair_significance.csv", SUB)
    for nm, d, pc in [("FFL triads", tri_df, "p_empirical"),
                      ("FFL seed pairs", pair_df, "p_empirical_upper"),
                      ("bi-fan pairs", bif_df, "p_empirical_upper")]:
        if len(d):
            print(f"  {nm}: {int((d[pc] <= 0.05).sum())} with p <= 0.05, "
                  f"{int(d['significant_q05'].sum())} surviving BH FDR "
                  f"(q <= 0.05) of {len(d)} tested")
    # ---- workbook-schema outputs (same null, same p; schema matches the
    #      All_FFLs / All_BiFans sheets so the two can be diffed row by row)
    two_step_summary = []
    for rows, stem, keycols, family in [
            (wb_ffl_rows, "ffl_triad_significance", ["X", "Y", "Z"], "FFL"),
            (wb_bifan_rows, "All_BiFans_Total", ["X1", "X2", "Z1", "Z2"],
             "BiFan_all_squares")]:
        d = pd.DataFrame(rows)
        if family != "FFL" and len(d):
            d = d.drop(columns=["Z_score_pure", "p_value_pure"])
        if len(d):
            d["q_value_BH"] = (d.groupby("condition")["p_value"]
                               .transform(lambda s: bh_fdr(s.values)).round(6)
                               if Q_SCOPE == "within_condition"
                               else bh_fdr(d["p_value"].values).round(6))
            d = two_step_flags(d, "n_hubs_in_triad" if "X" in d.columns
                               else "n_hubs_in_bifan")
            write_two_step(d, stem, keycols, family, two_step_summary)
    dbf = pd.DataFrame(wb_bifan_rows)
    if len(dbf):
        pure = dbf[dbf["is_pure_ID204"]].drop(columns=["is_pure_ID204"]).copy()
        # Pure family: a square counts as reproduced in a replicate only if it
        # reappears AND is pure there (induced ID 204), matching how the
        # observed motif was defined.
        pure["rand_hits"] = pure["rand_hits_pure"]
        pure["Z_score"] = pure["Z_score_pure"]
        pure["p_value"] = pure["p_value_pure"]
        pure = pure.drop(columns=["rand_hits_pure", "Z_score_pure",
                                  "p_value_pure"])
        if len(pure):
            pure["q_value_BH"] = (pure.groupby("condition")["p_value"]
                                  .transform(lambda s: bh_fdr(s.values)).round(6)
                                  if Q_SCOPE == "within_condition"
                                  else bh_fdr(pure["p_value"].values).round(6))
            pure = two_step_flags(pure, "n_hubs_in_bifan")
            write_two_step(pure, "bifan_pair_significance",
                           ["X1", "X2", "Z1", "Z2"], "BiFan_pure_ID204",
                           two_step_summary)
        print(f"  workbook tables: {len(wb_ffl_rows)} FFL triads, "
              f"{len(dbf)} bi-fan squares ({int(dbf['is_pure_ID204'].sum())} "
              f"pure ID204)")

    if two_step_summary:
        write(pd.DataFrame(two_step_summary), "step1_step2_summary.csv", SUB)
    write(pd.DataFrame(raw), "null_distributions.csv", SUB)
    write(pd.DataFrame(qual), "null_model_quality.csv", SUB)

    enr, hub, tri, bif = net_df, hub_df, tri_df, bif_df
    L = ["# Step 4 - motif enrichment against a degree-preserving null", "",
         f"Null model: {NULL_MODEL}, {N_RANDOMISATIONS} replicates per network.",
         "", "| condition | motif | observed | null mean +/- sd | fold | z | p | q |",
         "|---|---|---|---|---|---|---|---|"]
    for _, r in enr.iterrows():
        L.append(f"| {r['condition']} | {r['motif']} | {r['observed']} | "
                 f"{r['null_mean']} +/- {r['null_sd']} | {r['fold_enrichment']} | "
                 f"{r['z_score']} | {r['p_empirical_upper']} | "
                 f"{r['q_value_BH_within_condition']} |")
    L += ["", "## Hub-level FFL involvement", "",
          "| condition | hub | k_out | observed | null mean | fold | z | p |",
          "|---|---|---|---|---|---|---|---|"]
    for _, r in hub.iterrows():
        L.append(f"| {r['condition']} | {r['hub']} | {r['k_out']} | {r['observed']} | "
                 f"{r['null_mean']} | {r['fold_enrichment']} | {r['z_score']} | "
                 f"{r['p_empirical_upper']} |")
    if len(tri):
        sig = tri[tri["q_value_BH_within_condition"] <= 0.05]
        L += ["", "## Individually significant motifs", "",
              f"FFL triads passing BH FDR q <= 0.05: {len(sig)} of {len(tri)} "
              f"({100*len(sig)/len(tri):.1f}%).", ""]
        if len(bif):
            sb = bif[bif["q_value_BH_within_condition"] <= 0.05]
            L.append(f"Bi-fan pairs whose shared-target count passes q <= 0.05: "
                     f"{len(sb)} of {len(bif)}.")
    L += ["", "Reading: fold > 1 with small p means the motif is more frequent "
          "than the degree sequence alone produces. Per-motif rows name the "
          "specific circuits the null cannot reproduce. This rules out the "
          "degree-sequence explanation; it does not demonstrate selection.",
          "", "Null-model fidelity: see null_model_quality.csv - edges lost and "
          "both k_in and k_out error, measured on EVERY replicate and "
          "reported as _mean and _max. For NULL_MODEL='swap' these are zero "
          "by construction and replicates_degree_exact equals the replicate "
          "count.", ""]
    path = os.path.join(outdir(SUB), "step4_report.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(L))
    print(f"  wrote {path}")


if __name__ == "__main__":
    main()
