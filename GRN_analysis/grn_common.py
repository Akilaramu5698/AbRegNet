"""
grn_common.py -- shared loading / graph utilities for the condition-specific
Acinetobacter baumannii ATCC 19606 GRN analysis.

Edge files must have a header `source,target` (Regulator -> Target Gene).
Gene symbols are lower-cased and whitespace-stripped, so `FadR`, `fadr ` and
`fadr` collapse onto one node. Duplicate edges are removed.

Only the edge files that actually exist in DATA_DIR are loaded, so the same
code runs on 6 or 7 conditions without edits.
"""
from __future__ import annotations

import os
import pandas as pd
import networkx as nx

# condition label -> edge file name
GRN_FILES = {
    "amikacin":      "ami.grn.csv",
    "ciprofloxacin": "cip.grn.csv",
    "doripenem":     "dori.grn.csv",
    "meropenem":     "mero.grn.csv",
    "polymyxin":     "poly.grn.csv",
    "tigecycline":   "tige.grn.csv",
    "colistin_res":  "colis.resis.grn.csv",
}

# condition key -> short token used for loose (renamed-file) matching
GRN_TOKENS = {
    "amikacin": "ami", "ciprofloxacin": "cip", "doripenem": "dori",
    "meropenem": "mero", "polymyxin": "poly", "tigecycline": "tige",
    "colistin_res": "colis",
}

# Condition-specific hub panels: the regulators called as hubs in THAT
# network. Used by the hub-level tables in step1_dag_feedback.py and
# longest_paths.py. Lower-case, matching the parsed gene symbols.
HUBS_BY_CONDITION = {
    "amikacin":      ["beti", "rcsc", "fadr"],
    "ciprofloxacin": ["pgrr", "rcsc", "fadr"],
    "doripenem":     ["beti", "fadr", "dmlr"],
    "meropenem":     ["paax", "rcsc"],
    "polymyxin":     ["pgrr", "rcsc", "rpoh"],
    "tigecycline":   ["arsr2", "fadr"],
    "colistin_res":  ["lrp", "rcsc", "yybr", "rhas", "dmlr", "fadr"],
}


def hubs_for(condition: str) -> list:
    """Hub panel for one condition (empty list if the condition is unknown)."""
    return list(HUBS_BY_CONDITION.get(condition, []))


DATA_DIR = os.environ.get("GRN_DATA_DIR")        # None -> auto-search
OUT_DIR = os.environ.get("GRN_OUT_DIR", "results")

# Where to look for edge files when GRN_DATA_DIR is not set: the script's own
# folder, the working directory, and a data/ subfolder of either.
_HERE = os.path.dirname(os.path.abspath(__file__))
SEARCH_DIRS = [os.path.join(_HERE, "data"), _HERE, "data", "."]


def outdir(sub: str = "") -> str:
    """Create (if needed) and return an output directory."""
    p = os.path.join(OUT_DIR, sub) if sub else OUT_DIR
    os.makedirs(p, exist_ok=True)
    return p


def write(df: pd.DataFrame, filename: str, sub: str = "") -> str:
    """Save a table under the results directory; return its path."""
    path = os.path.join(outdir(sub), filename)
    df.to_csv(path, index=False)
    print(f"  wrote {path}  ({len(df)} rows)")
    return path


# ------------------------------------------------------------------ loading
def load_edges(path: str) -> pd.DataFrame:
    """Read one GRN csv -> tidy frame of unique directed edges."""
    df = pd.read_csv(path)
    cols = {str(c).lower().strip(): c for c in df.columns}
    if "source" in cols and "target" in cols:
        df = df[[cols["source"], cols["target"]]]
    else:                                   # fall back to first two columns
        df = df.iloc[:, :2]
    df.columns = ["source", "target"]
    for c in ("source", "target"):
        df[c] = df[c].astype(str).str.strip().str.lower()
    bad = {"", "nan", "none", "na"}
    df = df[~df["source"].isin(bad) & ~df["target"].isin(bad)]
    return df.drop_duplicates().reset_index(drop=True)


def build_graph(df: pd.DataFrame, name: str = "") -> nx.DiGraph:
    """Directed graph from an edge frame (self-loops kept, but counted)."""
    G = nx.DiGraph(name=name)
    G.add_edges_from(zip(df["source"], df["target"]))
    return G


def _norm(name: str) -> str:
    """Lower-case a filename and strip separators for loose matching."""
    return "".join(ch for ch in os.path.basename(name).lower() if ch.isalnum())


def resolve_files(data_dir: str, files: dict = None) -> dict:
    """{condition: path} for one directory.

    Exact filename matches win; a csv whose name merely starts with the
    condition token (e.g. `amikacin.csv`, `cip_grn.csv`) is accepted as a
    fallback, so renamed files are still picked up.
    """
    files = files or GRN_FILES
    if not os.path.isdir(data_dir):
        return {}
    present = [f for f in os.listdir(data_dir) if f.lower().endswith(".csv")]
    by_norm = {_norm(f): f for f in present}
    found, used = {}, set()
    for cond, fn in files.items():                       # pass 1: exact name
        if _norm(fn) in by_norm:
            f = by_norm[_norm(fn)]
            found[cond] = os.path.join(data_dir, f)
            used.add(f)
    for cond, fn in files.items():                       # pass 2: token prefix
        if cond in found:
            continue
        tok = GRN_TOKENS.get(cond, _norm(fn)[:4])
        hits = [f for f in present if f not in used and _norm(f).startswith(tok)]
        if len(hits) == 1:
            found[cond] = os.path.join(data_dir, hits[0])
            used.add(hits[0])
    return found


def find_data_dir(data_dir: str = None, files: dict = None) -> tuple:
    """Pick the directory holding the edge files; return (dir, {cond: path})."""
    cands = []
    for d in [data_dir, DATA_DIR] + SEARCH_DIRS:
        if d and os.path.abspath(d) not in [os.path.abspath(c) for c in cands]:
            cands.append(d)
    best, best_hits = None, {}
    for d in cands:
        hits = resolve_files(d, files)
        if len(hits) > len(best_hits):
            best, best_hits = d, hits
    if not best_hits:
        seen = []
        for d in cands:
            if os.path.isdir(d):
                seen += [f"{d}: " + (", ".join(f for f in os.listdir(d)
                                               if f.lower().endswith(".csv")) or "no csv files")]
        raise FileNotFoundError(
            "no GRN edge files found. Searched:\n  "
            + "\n  ".join(os.path.abspath(d) for d in cands)
            + ("\nCSV files seen:\n  " + "\n  ".join(seen) if seen else "")
            + "\nPut the edge files next to the scripts, or set GRN_DATA_DIR."
        )
    return best, best_hits


def load_all(data_dir: str = None, files: dict = None, verbose: bool = True) -> dict:
    """{condition: DiGraph} for every edge file found. Searches data_dir,
    $GRN_DATA_DIR, the script folder, the working directory and data/."""
    files = files or GRN_FILES
    chosen, paths = find_data_dir(data_dir, files)
    graphs = {cond: build_graph(load_edges(paths[cond]), name=cond)
              for cond in files if cond in paths}
    if verbose:
        print(f"  data directory: {os.path.abspath(chosen)}")
        for cond, G in graphs.items():
            print(f"  {cond:<14} {os.path.basename(paths[cond]):<22} "
                  f"{G.number_of_nodes():>5} nodes  {G.number_of_edges():>6} edges")
        missing = [c for c in files if c not in paths]
        if missing:
            print(f"  (no file matched, skipped: {', '.join(missing)})")
    return graphs


# ------------------------------------------------------------- node roles
def regulators(G: nx.DiGraph) -> set:
    """Active regulators: out-degree > 0 ignoring self-loops."""
    return {n for n in G if G.out_degree(n) - (1 if G.has_edge(n, n) else 0) > 0}


def terminals(G: nx.DiGraph) -> set:
    """Terminal target nodes: out-degree == 0."""
    return {n for n in G if G.out_degree(n) == 0}


def tf_tf_edges(G: nx.DiGraph) -> list:
    """Regulator -> regulator edges (both endpoints active, self-loops dropped)."""
    tfs = regulators(G)
    return [(u, v) for u, v in G.edges() if u in tfs and v in tfs and u != v]


def tf_subgraph(G: nx.DiGraph) -> nx.DiGraph:
    """Induced subgraph on active regulators = the TF->TF backbone."""
    H = G.subgraph(regulators(G)).copy()
    H.remove_edges_from(list(nx.selfloop_edges(H)))
    return H


def degree_frame(G: nx.DiGraph) -> pd.DataFrame:
    """Per-node k_in / k_out / k_tot, split by TF vs leaf out-edges."""
    tfs = regulators(G)
    rows = []
    for n in G.nodes():
        loop = 1 if G.has_edge(n, n) else 0
        kout_tf = sum(1 for v in G.successors(n) if v in tfs and v != n)
        rows.append({
            "node": n,
            "k_in": G.in_degree(n),
            "k_out": G.out_degree(n),
            "k_tot": G.in_degree(n) + G.out_degree(n),
            "k_out_to_tf": kout_tf,
            "k_out_to_leaf": G.out_degree(n) - kout_tf - loop,
            "is_regulator": n in tfs,
        })
    return (pd.DataFrame(rows)
              .sort_values(["k_out", "k_in"], ascending=False)
              .reset_index(drop=True))


def targets_of(G: nx.DiGraph, n: str) -> set:
    """Direct regulon of node n (excluding a self-loop)."""
    return set(G.successors(n)) - {n}
