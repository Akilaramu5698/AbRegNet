#!/usr/bin/env python3
"""
Tier architecture of the ENTIRE regulatory network  (all significant FFL triads + bifans, all antibiotics)
-----------------------------------------------------------------------------------------------------------

Inputs : ffl_triad_significance.csv, bifan_pair_significance.csv
Outputs: network_tiers.png / .pdf                  regulators only, 3 tiers
         network_tiers_with_targets.png / .pdf     + target genes in a band below Tier 3
         network_tiers_nodes.csv                   every gene: tier, class, in/out degree, #targets, antibiotics
         network_tiers_regulator_edges.csv         every regulator-regulator edge (FFL pair X->Y, bifan pair X1-X2)
         network_tiers_target_edges.csv            every regulator -> target edge
         network_tiers_summary.csv                 members of each tier

Tier definition (FFL arrows X -> Y of the significant FFL triads):
  k_in_hub = arrows received from OTHER HUBS,  k_out_tf = arrows sent to other regulators (hub or non-hub TF)
  Tier 1 = k_in_hub == 0                       (no other hub regulates it; upstream)
  Tier 2 = k_in_hub > 0 and k_out_tf > 0       (hub-controlled and regulates other regulators; intermediate)
  Tier 3 = k_in_hub > 0 and k_out_tf == 0      (hub-controlled, regulates no regulator; end-point)
  bifan-only regulators (no FFL arrow) therefore have k_in_hub = 0 -> Tier 1 (dashed outline)


usage:  python regulatory_network_tiers.py [strict|lenient|qonly] [--dir folder] [--min-abx 1] [--shared-min 5]
"""
import argparse, math, collections, re
from pathlib import Path
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

ap = argparse.ArgumentParser()
ap.add_argument('mode', nargs='?', default='strict', choices=['strict', 'lenient', 'qonly'])
ap.add_argument('--dir', default=str(Path(__file__).resolve().parent))
ap.add_argument('--min-abx', type=int, default=1, help='regulator pairs significant in >= this many antibiotics (1 = entire network)')
ap.add_argument('--shared-min', type=int, default=5, help='target is "shared" if linked to >= this many regulators')
ap.add_argument('--tier-labels', action='store_true', help='show the "Tier 1 upstream regulators ..." labels at the left of the bands (hidden by default)')
ap.add_argument('--dpi', type=int, default=200)
a = ap.parse_args()
MODE, Q_CUT, REQUIRE_MASTER_HUB = a.mode, 0.05, True

plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
                     "font.weight": "bold", "pdf.fonttype": 42, "ps.fonttype": 42})

# ----------------------------------------------------------------- load + filter (identical to the pipeline)
ffl = pd.read_csv(f"{a.dir}/ffl_triad_significance.csv")
bif = pd.read_csv(f"{a.dir}/bifan_pair_significance.csv")
ALL_ABX = set(ffl.condition) | set(bif.condition); N_ABX = len(ALL_ABX)      # total antibiotics in the data
hubs = set()
for g, h in [("X1", "X1_is_hub"), ("X2", "X2_is_hub"), ("Z1", "Z1_is_hub"), ("Z2", "Z2_is_hub")]:
    hubs |= set(bif.loc[bif[h], g])
hubs |= set(ffl.loc[ffl.X_is_hub, "X"]) | set(ffl.loc[ffl.Y_is_hub, "Y"]) | set(ffl.loc[ffl.Z_is_hub, "Z"])

def z_ok(z):
    fin = np.isfinite(z)
    if MODE == "strict":  return z >= 2
    if MODE == "lenient": return ~fin | (z > 2)
    return pd.Series(True, index=z.index)
def sig(d):
    ok = d.q_value_BH < Q_CUT
    return ok & d.contains_master_hub if REQUIRE_MASTER_HUB else ok
ffl = ffl[sig(ffl) & z_ok(ffl.Z_score)].copy()
bif = bif[sig(bif) & z_ok(bif.Z_score)].copy()
key = lambda x, y: "|".join(sorted([x, y]))
ffl["pair"] = [key(x, y) for x, y in zip(ffl.X, ffl.Y)]
bif["pair"] = [key(x, y) for x, y in zip(bif.X1, bif.X2)]

# ----------------------------------------------------------------- whole network (optionally >= min-abx antibiotics)
abx_of = collections.defaultdict(set)
for r in ffl.itertuples(): abx_of[r.pair].add(r.condition)
for r in bif.itertuples(): abx_of[r.pair].add(r.condition)
keep = {p for p, s in abx_of.items() if len(s) >= a.min_abx}
fsel, bsel = ffl[ffl.pair.isin(keep)], bif[bif.pair.isin(keep)]
print(f"[{MODE}] regulator pairs: {len(keep)}   FFL triads: {len(fsel)}   bifans: {len(bsel)}")

# regulator-regulator pair table
PR = {}                                                    # pair -> dict
for r in fsel.itertuples():
    d = PR.setdefault(r.pair, dict(src=None, dst=None, kinds=set(), n_ffl=0, n_bif=0, abx=set()))
    if d['src'] is None: d['src'], d['dst'] = r.X, r.Y
    d['kinds'].add('FFL'); d['n_ffl'] += 1; d['abx'].add(r.condition)
for r in bsel.itertuples():
    d = PR.setdefault(r.pair, dict(src=None, dst=None, kinds=set(), n_ffl=0, n_bif=0, abx=set()))
    d['kinds'].add('bifan'); d['n_bif'] += 1; d['abx'].add(r.condition)

regs = set()
for p in PR: regs |= set(p.split("|"))
tg_edges = collections.defaultdict(set)                    # (regulator, target) -> motif kinds
for r in fsel.itertuples():
    tg_edges[(r.X, r.Z)].add('FFL'); tg_edges[(r.Y, r.Z)].add('FFL')
for r in bsel.itertuples():
    for x in (r.X1, r.X2):
        for z in (r.Z1, r.Z2): tg_edges[(x, z)].add('bifan')
targets = sorted({z for (_, z) in tg_edges} - regs)
tlink = collections.defaultdict(set)                       # target -> regulators linked to it
for (x, z) in tg_edges:
    if z in targets: tlink[z].add(x)
shared = {g for g in targets if len(tlink[g]) >= a.shared_min}
ntg = collections.defaultdict(set)                         # all distinct target genes (Z) of each regulator (incl. targets that are regulators)
for (x, z) in tg_edges: ntg[x].add(z)
print(f"regulators {len(regs)} (hubs {len([g for g in regs if g in hubs])}), targets {len(targets)} (shared >= {a.shared_min}: {len(shared)})")

# ----------------------------------------------------------------- tiers
outn, inn = collections.defaultdict(set), collections.defaultdict(set)
for p, d in PR.items():
    if d['src'] is not None: outn[d['src']].add(d['dst']); inn[d['dst']].add(d['src'])
partners = collections.defaultdict(set)
for p in PR:
    x, y = p.split("|"); partners[x].add(y); partners[y].add(x)
# ---- tier rule (hub-based k_in, k_out_tf) -------------------------------------------------------------
#   k_in_hub = number of FFL arrows INTO the regulator that come from OTHER HUBS
#   k_out_tf = number of FFL arrows OUT of the regulator to other regulators (hub or non-hub TFs)
#   Tier 1 : k_in_hub == 0                      (not controlled by any other hub)
#   Tier 2 : k_in_hub  > 0  and  k_out_tf > 0   (controlled by a hub and regulates other regulators)
#   Tier 3 : k_in_hub  > 0  and  k_out_tf == 0  (controlled by a hub, regulates no regulator = end-point)
ARROWS = {(d['src'], d['dst']) for d in PR.values() if d['src'] is not None and d['src'] != d['dst']}
k_in_hub, k_out_tf = collections.Counter(), collections.Counter()
for s_, d_ in ARROWS:
    k_out_tf[s_] += 1
    if s_ in hubs: k_in_hub[d_] += 1
tier = {g: (1 if k_in_hub[g] == 0 else 2 if k_out_tf[g] > 0 else 3) for g in regs}
bifan_only = sorted(g for g in regs if g not in {x for e in ARROWS for x in e})   # regulators with no FFL arrow at all (bifan-only) -> dashed outline
print("tiers:", {k: len([g for g in regs if tier[g] == k]) for k in (1, 2, 3)}, "| bifan-only regulators (dashed):", len(bifan_only))

abx_reg = collections.defaultdict(set)                       # antibiotics in which each regulator occurs (in any significant pair)
for p, d in PR.items():
    for g in p.split("|"): abx_reg[g] |= d['abx']

# ----------------------------------------------------------------- layout
FW, W, ROWMAX, LANE = 32.0, 31.0, 20, 2.9
# publication-size text (figure is 32 in wide, so fonts are large)
FS_GENE_HUB, FS_GENE, FS_CAP, FS_LEG, FS_TITLE = 17, 13.5, 10, 20, 30
SZ = {'hub': 4800, 'reg': 3000, 'shared': 70, 'tgt': 45}
COL = {'hub': '#F5A623', 'reg': '#82E0AA', 'shared': '#8E44AD', 'tgt': '#D5D8DC'}
cls = lambda g: 'hub' if (g in regs and g in hubs) else 'reg' if g in regs else 'shared' if g in shared else 'tgt'
RED, BLK, BLU = '#C0392B', '#111111', '#2471A3'
NOMEN_MAP = {
    "beti": "betI", "dmlr": "dmlR", "arsr2": "arsR2",
    "fadr": "fadR", "lrp": "lrp", "pgrr": "pgrR", "rhas": "rhaS",
    "rcsc": "rcsC", "paax": "paaX", "yybr": "yybR", "rpoh": "rpoH",
    "mmgc": "mmgC", "cyda": "cydA", "cydb": "cydB", "fabg": "fabG",
    "fabz": "fabZ", "fadh": "fadH", "fadj": "fadJ", "catd": "catD",
    "btub": "btuB", "hmp": "hmp", "oprd": "oprD", "pilt": "pilT",
    "pilb": "pilB", "pile": "pilE", "pilv": "pilV", "pilo": "pilO",
    "epsf": "epsF", "epse": "epsE", "bla-oxa-256": "blaOXA-256",
    "lola": "lolA", "lolb": "lolB", "lold": "lolD", "tola": "tolA",
    "tolq": "tolQ", "tolc": "tolC", "mltb": "mltB", "slt": "slt",
    "loip": "loiP", "dsba": "dsbA", "olsb": "olsB", "mexb": "mexB",
    "emrb": "emrB", "emra": "emrA", "gstb": "gstB", "prse": "prsE",
    "mcpb": "mcpB", "tfpa": "tfpA", "cpda": "cpdA", "pgab": "pgaB",
    "grxd": "grxD", "pmpa": "pmpA", "uppp": "uppP", "deda": "dedA",
    "ftsa": "ftsA", "kdpb": "kdpB", "prop": "proP", "ahpc": "ahpC",
    "ahpf": "ahpF", "menf": "menF", "naip": "naiP", "iuta": "iutA",
    "rsua": "rsuA", "dpgd": "dpgD", "phac": "phaC", "recn": "recN",
    "sbcc": "sbcC", "xerd": "xerD", "para": "parA", "comc": "comC",
    "smf-1": "smf-1", "bepc": "bepC", "mdta": "mdtA", "czca": "czcA",
    "dhba": "dhbA", "irtb": "irtB", "frr": "frr", "iles": "ileS",
    "gata": "gatA", "gatb": "gatB", "quea": "queA", "secd": "secD",
    "secf": "secF", "lpxa": "lpxA", "plsc": "plsC", "kdsa": "kdsA",
    "acpp": "acpP", "mscl": "mscL", "mscs": "mscS", "mepa": "mepA",
    "baua": "bauA", "bauf": "bauF", "fhua": "fhuA", "bfd": "bfd",
    "rpla": "rplA", "nusb": "nusB", "rhlb": "rhlB", "rtcb": "rtcB",
    "cspv": "cspV", "sirb": "sirB", "pcaj": "pcaJ", "pyrd": "pyrD",
    "meth": "metH", "eno": "eno", "ppk": "ppk",
    # Partner & local TFs
    "adet": "adeT", "acrr": "acrR", "baer": "baeR", "tetc": "tetC",
    "phop": "phoP", "phor": "phoR", "pled": "pleD", "atoc": "atoC",
    "gltr": "gltR", "angr": "angR", "hcar": "hcaR", "kstr2": "kstR2",
    "fur": "fur", "cmpr": "cmpR", "gcva": "gcvA", "gntr": "gntR",
    "nimr": "nimR", "sutr": "sutR", "norg": "norG", "pcau": "pcaU",
    "pobr": "pobR", "decr": "decR", "acnr": "acnR", "hdfr": "hdfR",
    "yjie": "yjiE", "ompr": "ompR", "envz": "envZ", "bola": "bolA",
    "arac": "araC", "perr": "perR", "soxr": "soxR", "pmra": "pmrA",
    "daur": "dauR", "pcpr": "pcpR", "tas": "tas", "tusd": "tusD",
    "ydja": "ydjA", "ylbk": "ylbK", "tesa": "tesA", "mrcb": "mrcB",
}
def nomen(g):
    """gene label: your NOMEN_MAP first; otherwise bacterial nomenclature rule (3 lowercase letters + capital 4th letter, e.g. fadr -> fadR)"""
    if g in NOMEN_MAP: return NOMEN_MAP[g]
    m = re.fullmatch(r'([a-z]{3})([a-z])(\d*)', g) or re.fullmatch(r'([a-z]{4})([a-z])()', g)     # 4- or 5-letter symbols (e.g. pdtar -> pdtaR)
    return m.group(1) + m.group(2).upper() + m.group(3) if m else g
RED_STRONG, RED_LIGHT = '#C0392B', '#A6A6A6'     # hub-hub edges (strong red) / other regulator edges (grey)

def spread(nodes):
    nodes = list(nodes); n = len(nodes)
    return {g: (-W/2 + W*(i + .5)/n) for i, g in enumerate(nodes)}
tx = {}
for k in (1, 2, 3): tx.update(spread(sorted(g for g in regs if tier[g] == k)))
for _ in range(14):                                         # barycentre ordering (fewer crossings)
    for k in (2, 3):
        nodes = [g for g in regs if tier[g] == k]
        des = {g: (np.mean([tx[q] for q in inn[g]]) if inn[g] else np.mean([tx[q] for q in partners[g]]) + 0.3) for g in nodes}
        tx.update(spread(sorted(nodes, key=lambda g: des[g])))
    nodes = [g for g in regs if tier[g] == 1]
    des = {g: np.mean([tx[q] for q in outn[g]]) if outn[g] else np.mean([tx[q] for q in partners[g]]) for g in nodes}
    tx.update(spread(sorted(nodes, key=lambda g: des[g])))
if 'dmlr' in tx and 'fadr' in tx and tier['dmlr'] == tier['fadr'] and tx['dmlr'] > tx['fadr']: tx['dmlr'], tx['fadr'] = tx['fadr'], tx['dmlr']

rows_k = {k: max(1, math.ceil(sum(1 for g in regs if tier[g] == k) / ROWMAX)) for k in (1, 2, 3)}
bh = {k: 3.6 + (rows_k[k] - 1) * LANE for k in (1, 2, 3)}
TY, cur = {}, 0.0
for k in (1, 2, 3):
    TY[k] = cur - bh[k]/2; cur -= bh[k] + 2.0
BAND_Y, BAND_H = cur - 2.4 + 1.0, 4.4
tpos = {}
for k in (1, 2, 3):
    nodes = sorted((g for g in regs if tier[g] == k), key=lambda g: tx[g])
    for i, g in enumerate(nodes):
        r_ = i % rows_k[k]; tpos[g] = (tx[g], TY[k] + (r_ - (rows_k[k] - 1)/2) * LANE)

# target band
rng = np.random.default_rng(3)
X0 = np.array([np.mean([tx[o] for o in tlink[g]]) for g in targets]); Pt = np.c_[X0, BAND_Y + rng.uniform(-BAND_H/2.4, BAND_H/2.4, len(targets))]
dmin = 0.30
for it in range(900):
    Dd = Pt[:, None, :] - Pt[None, :, :]; dist = np.hypot(Dd[..., 0], Dd[..., 1]) + 1e-9; np.fill_diagonal(dist, 1e9)
    Pt += 0.5*((((dmin - dist)*(dist < dmin))/dist)[..., None]*Dd).sum(axis=1)*0.5
    Pt[:, 0] += 0.012*(X0 - Pt[:, 0]); Pt[:, 1] = np.clip(Pt[:, 1], BAND_Y - BAND_H/2, BAND_Y + BAND_H/2)
    Pt[:, 0] = np.clip(Pt[:, 0], -W/2 - 1.0, W/2 + 1.0)
tgpos = {g: tuple(p) for g, p in zip(targets, Pt)}

# ----------------------------------------------------------------- figures
def tier_figure(with_targets, fname):
    P = dict(tpos); P.update(tgpos if with_targets else {})
    ytop = 0.4; ybot = (BAND_Y - BAND_H/2 - 1.8) if with_targets else (TY[3] - bh[3]/2 - 1.2)
    fig, ax = plt.subplots(figsize=(FW, (ytop - ybot) * 0.92 + 3.5)); ax.axis('off')
    xl = (-W/2 - (5.2 if a.tier_labels else 1.8), W/2 + 2.2); ax.set_xlim(*xl); ax.set_ylim(ybot, ytop + 1.1)
    upi = (xl[1] - xl[0]) / FW
    rad_ = lambda g: math.sqrt(SZ[cls(g)]/math.pi)
    bands = [(1, "Tier 1\nupstream\nregulators", '#EEF3FA'), (2, "Tier 2\nintermediate\nregulators", '#F6F1E7'), (3, "Tier 3\nend-point\nregulators", '#EAF5EC')]
    for k, lab_, colb in bands:
        n_k = sum(1 for g in regs if tier[g] == k)
        ax.add_patch(Rectangle((-W/2 - 1.2, TY[k] - bh[k]/2), W + 2.4, bh[k], color=colb, zorder=0))
        if a.tier_labels: ax.text(-W/2 - 1.4, TY[k], f"{lab_}\n(n = {n_k})", ha='right', va='center', fontsize=17, fontweight='bold', color='#444444')
    if with_targets:
        ax.add_patch(Rectangle((-W/2 - 1.2, BAND_Y - BAND_H/2 - .4), W + 2.4, BAND_H + .8, color='#F4F4F4', zorder=0))
        ax.text(-W/2 - 1.4, BAND_Y, f"Targets\n(n = {len(targets)},\nunlabelled)", ha='right', va='center', fontsize=17, fontweight='bold', color='#444444')
    ax.set_title(f"Tier architecture of the entire regulatory network  ({len(regs)} regulators, {len(PR)} regulator pairs, "
                 f"{len(fsel)} FFL triads, {len(bsel)} bifans" + (f", {len(targets)} targets)" if with_targets else ")"),
                 fontsize=FS_TITLE, fontweight='bold', pad=22)
    # --- regulator target edges (thin, drawn first, only in the targets figure)
    if with_targets:
        segF, segB = [], []
        for (x, z), kinds in tg_edges.items():
            if z not in P or x not in P: continue
            (segF if 'FFL' in kinds else segB).append([P[x], P[z]])
        ax.add_collection(LineCollection(segB, colors=BLU, linewidths=0.35, alpha=0.10, zorder=1))
        ax.add_collection(LineCollection(segF, colors=BLK, linewidths=0.45, alpha=0.16, zorder=1))
    # --- regulator-regulator edges (FFL + bifan together, one edge per pair):
    #     hub-hub pair  -> strong red ;  any pair with a non-hub regulator -> light colour
    #     FFL pair -> arrow (X -> Y) ;  bifan-only pair -> plain line (no arrow).   All lines have the same width.
    HH = lambda p: all(g in hubs for g in p.split("|"))
    for p, d in sorted(PR.items(), key=lambda kv: HH(kv[0])):            # light edges first, strong hub-hub edges on top
        x, y = p.split("|"); isF = d['src'] is not None; hh = HH(p)
        u, v = (d['src'], d['dst']) if isF else (x, y)
        same = tier[u] == tier[v]; rr = (0.28 if tx[u] < tx[v] else -0.28) if same else 0.04
        col, al, z = (RED_STRONG, 0.95, 4) if hh else (RED_LIGHT, 0.50, 2)
        ax.add_patch(FancyArrowPatch(P[u], P[v], arrowstyle='-|>' if isF else '-', mutation_scale=20 if isF else 1, lw=1.9, color=col, alpha=al,
                     shrinkA=rad_(u) + 1.5, shrinkB=rad_(v) + 1.5, connectionstyle=f'arc3,rad={rr}', zorder=z))
    # --- nodes
    if with_targets:
        for g in targets:
            k_ = cls(g); ax.scatter(*P[g], s=SZ[k_], c=COL[k_], edgecolors='k', linewidths=0.4, zorder=3)
    for g in regs:
        k_ = cls(g)
        ax.scatter(*P[g], s=SZ[k_], c=COL[k_], edgecolors='k', linewidths=2.0, linestyle='--' if g in bifan_only else '-', zorder=5)
        ax.text(*P[g], nomen(g), ha='center', va='center', fontsize=FS_GENE_HUB if k_ == 'hub' else FS_GENE, fontstyle='italic', fontweight='bold', zorder=6)   # italic, gene nomenclature
        ax.text(P[g][0], P[g][1] - rad_(g)*upi/72 - 0.42, f"{len(ntg[g])} targets\n{len(abx_reg[g])}/{N_ABX} antibiotics", ha='center', va='center',
                fontsize=FS_CAP, fontweight='bold', color='#333333', zorder=6, linespacing=1.05,
                bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.85))
    n_hh_f = sum(1 for p, d in PR.items() if d['src'] is not None and HH(p)); n_hh_b = sum(1 for p, d in PR.items() if d['src'] is None and HH(p))
    n_o_f = sum(1 for p, d in PR.items() if d['src'] is not None and not HH(p)); n_o_b = sum(1 for p, d in PR.items() if d['src'] is None and not HH(p))
    leg = [Line2D([0], [0], color=RED_STRONG, lw=3, marker='>', markevery=[1], ms=11, label=f"hub–hub FFL edge X→Y (strong red arrow; {n_hh_f})"),
           Line2D([0], [0], color=RED_STRONG, lw=3, label=f"hub–hub bifan edge (strong red line, no arrow; {n_hh_b})"),
           Line2D([0], [0], color=RED_LIGHT, lw=3, marker='>', markevery=[1], ms=11, label=f"other FFL edge (grey arrow; {n_o_f})"),
           Line2D([0], [0], color=RED_LIGHT, lw=3, label=f"other bifan edge (grey line; {n_o_b})"),
           Line2D([0], [0], marker='o', ls='', mfc=COL['hub'], mec='k', ms=15, label=f"hub regulator ({len([g for g in regs if g in hubs])})"),
           Line2D([0], [0], marker='o', ls='', mfc=COL['reg'], mec='k', ms=12, label=f"other regulator ({len([g for g in regs if g not in hubs])})"),
           Line2D([0], [0], marker='o', ls='', mfc='white', mec='k', mew=2, ms=12, label=f"dashed outline = bifan-only regulator, no FFL arrow ({len(bifan_only)})")]
    if with_targets:
        leg += [Line2D([0], [0], color=BLK, lw=1.8, label='FFL triad edges X→Z, Y→Z (regulator → target)'), Line2D([0], [0], color=BLU, lw=1.8, alpha=.5, label='bifan edges to targets'),
                Line2D([0], [0], marker='o', ls='', mfc=COL['shared'], mec='k', ms=8, label=f'shared target (≥{a.shared_min} regulators; {len(shared)})'),
                Line2D([0], [0], marker='o', ls='', mfc=COL['tgt'], mec='k', ms=7, label=f'other target ({len(targets) - len(shared)})')]
    ax.legend(handles=leg, loc='upper center', bbox_to_anchor=(0.5, -0.005), ncol=3, frameon=True, prop={'weight': 'bold', 'size': FS_LEG}, markerscale=1.3, borderpad=0.9, labelspacing=0.7, handlelength=2.6)
    plt.tight_layout(); plt.savefig(fname + '.png', dpi=a.dpi, bbox_inches='tight'); plt.savefig(fname + '.pdf', bbox_inches='tight'); plt.close()

tier_figure(False, f"{a.dir}/network_tiers")
tier_figure(True,  f"{a.dir}/network_tiers_with_targets")

# ----------------------------------------------------------------- tables
out = f"{a.dir}/network_tiers"
nodes = []
for g in sorted(regs | set(targets)):
    isr = g in regs
    nodes.append(dict(node=g, tier=f"Tier {tier[g]}" if isr else 'target',
                      node_class={'hub': 'hub regulator', 'reg': 'other regulator', 'shared': 'shared target', 'tgt': 'other target'}[cls(g)],
                      bifan_only_regulator=(g in bifan_only) if isr else '', k_in_hub=k_in_hub[g] if isr else '', k_out_tf=k_out_tf[g] if isr else '',
                      n_out_FFL_arrows=len(outn[g]) if isr else '', n_in_FFL_arrows=len(inn[g]) if isr else '',
                      n_regulator_partners=len(partners[g]) if isr else '', n_targets=len(ntg[g]) if isr else '',
                      n_regulators_linked=len(tlink[g]) if not isr else '', antibiotics=';'.join(sorted(abx_reg[g])) if isr else ''))
pd.DataFrame(nodes).sort_values(['tier', 'node_class', 'node']).to_csv(out + '_nodes.csv', index=False)
pd.DataFrame([dict(source=(d['src'] or p.split('|')[0]), target=(d['dst'] or p.split('|')[1]), directed=d['src'] is not None, motif=';'.join(sorted(d['kinds'])),
                   n_ffl_triads=d['n_ffl'], n_bifans=d['n_bif'], n_antibiotics=len(d['abx']), antibiotics=';'.join(sorted(d['abx'])),
                   source_tier=tier[d['src'] or p.split('|')[0]], target_tier=tier[d['dst'] or p.split('|')[1]]) for p, d in PR.items()]
             ).sort_values(['n_antibiotics', 'n_ffl_triads', 'n_bifans'], ascending=False).to_csv(out + '_regulator_edges.csv', index=False)
pd.DataFrame([dict(regulator=x, target=z, motif=';'.join(sorted(k)), target_is_regulator=z in regs) for (x, z), k in tg_edges.items()]
             ).sort_values(['regulator', 'target']).to_csv(out + '_target_edges.csv', index=False)
pd.DataFrame([dict(tier=f"Tier {k}", n_regulators=len(m), n_hubs=len([g for g in m if g in hubs]), n_bifan_only=len([g for g in m if g in bifan_only]),
                   regulators=';'.join(sorted(m)))
              for k in (1, 2, 3) for m in [[g for g in regs if tier[g] == k]]]).to_csv(out + '_summary.csv', index=False)
print('saved', out + '.png/.pdf, _with_targets, _nodes.csv, _regulator_edges.csv, _target_edges.csv, _summary.csv')
