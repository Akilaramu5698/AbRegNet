#!/usr/bin/env python3
"""

Usage:  python build_fig2_curated.py [input_folder]
Outputs: results/fig2/{Cytoscape_Fig2_FullGRN_3Tier_CuratedClasses.xgmml,
         Fig2_Nodes_CuratedClasses.csv, Fig2_Edges_CuratedClasses.csv,
         Fig2_Class_Legend.csv, Fig2_render.png/.svg,
         Fig2_poster.png/.svg, Fig2_Outer_Labels.csv}

"""

import csv
import math
import sys
import xml.etree.ElementTree as ET
from collections import OrderedDict, defaultdict
from pathlib import Path

from xml.sax.saxutils import escape as xml_esc

# =========================================================================
# VENDORED verbatim from the user's fullgrn_data_and_xgmml.py
# (HUB_TIERS, SPECIAL_CASE, fmt_gene, write_xgmml) - one LOCAL ADDITION in
# write_xgmml, marked inline, to pass extra Cytoscape visual properties.
# =========================================================================
# Verbatim excerpts from fullgrn_data_and_xgmml.py (HUB_TIERS, SPECIAL_CASE, fmt_gene, write_xgmml)

HUB_TIERS = {
    "beti": ("Top_Initiator", 1, "betI"),
    "dmlr": ("Top_Initiator", 1, "dmlR"),
    "arsr2": ("Top_Initiator", 1, "arsR2"),
    "fadr": ("Mid_Integrator", 2, "fadR"),
    "lrp": ("Mid_Integrator", 2, "lrp"),
    "pgrr": ("Mid_Integrator", 2, "pgrR"),
    "rhas": ("Mid_Integrator", 2, "rhaS"),
    "rcsc": ("Bottom_Executor", 3, "rcsC"),
    "paax": ("Bottom_Executor", 3, "paaX"),
    "yybr": ("Bottom_Executor", 3, "yybR"),
    "rpoh": ("Bottom_Executor", 3, "rpoH"),
}

SPECIAL_CASE = {
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


def fmt_gene(g: str) -> str:
    gl = g.lower()
    if gl in SPECIAL_CASE:
        return SPECIAL_CASE[gl]
    if len(gl) == 4 and gl[3].isalpha():
        return gl[:3] + gl[3].upper()
    if len(gl) == 5 and gl[3].isalpha() and gl[4].isdigit():
        return gl[:3] + gl[3].upper() + gl[4]
    return gl


def write_xgmml(
    filepath: Path,
    title: str,
    nodes_spec: list,
    edges_spec: list,
    bg_color: str = "#F8FAFC",
    center_x: float = 0.0,
    center_y: float = 0.0,
    zoom: float = 0.65,
):
    """
    Write a standards-compliant Cytoscape 3.x XGMML 1.0 file with embedded (x, y) coordinates,
    node shapes/colors/transparencies, selective labels, and edge widths/transparencies/arrows.

    nodes_spec items: dict with keys:
      id (str), name (str), spotlight_label (str), x (float), y (float),
      w (float), h (float), shape (HEXAGON|OCTAGON|DIAMOND|ELLIPSE|ROUND_RECTANGLE),
      fill (hex), outline (hex), border_w (float), opacity (0..1),
      label_color (hex), font_size (int), attrs (dict of extra attributes)
    edges_spec items: dict with keys:
      source (str), target (str), interaction (str), stroke (hex),
      width (float), opacity (0..1), dashed (bool), target_arrow (bool), attrs (dict)
    """
    id_map = {}
    lines = [
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
        f'<graph id="{abs(hash(title)) % 900000 + 100000}" label="{xml_esc(title)}" '
        'directed="1" cy:documentVersion="3.0" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" '
        'xmlns:xlink="http://www.w3.org/1999/xlink" '
        'xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" '
        'xmlns:cy="http://www.cytoscape.org" '
        'xmlns="http://www.cs.rpi.edu/XGMML">',
        f'  <att name="shared name" value="{xml_esc(title)}" type="string"/>',
        f'  <att name="name" value="{xml_esc(title)}" type="string"/>',
        '  <graphics>',
        f'    <att name="NETWORK_BACKGROUND_PAINT" value="{bg_color}" type="string"/>',
        f'    <att name="NETWORK_CENTER_X_LOCATION" value="{center_x:.2f}" type="real"/>',
        f'    <att name="NETWORK_CENTER_Y_LOCATION" value="{center_y:.2f}" type="real"/>',
        f'    <att name="NETWORK_SCALE_FACTOR" value="{zoom:.4f}" type="real"/>',
        '  </graphics>',
    ]

    for idx, nd in enumerate(nodes_spec, start=1):
        nid = str(idx)
        id_map[nd["id"]] = nid
        slbl = nd.get("spotlight_label", "")
        alpha_255 = max(5, min(255, int(round(nd.get("opacity", 1.0) * 255))))
        fsz = nd.get("font_size", 12)
        lcol = nd.get("label_color", "#0F172A")
        # Notice: setting label to spotlight_label keeps ghosted background nodes completely label-free on default import,
        # while gene_symbol & shared_name preserve full searchability of every gene in Cytoscape!
        node_xml_label = slbl if slbl else ""
        lines.append(f'  <node id="{nid}" label="{xml_esc(node_xml_label)}">')
        lines.append(f'    <att name="gene_id" value="{xml_esc(nd["id"])}" type="string"/>')
        lines.append(f'    <att name="gene_symbol" value="{xml_esc(nd["name"])}" type="string"/>')
        lines.append(f'    <att name="spotlight_label" value="{xml_esc(slbl)}" type="string"/>')
        for k, v in nd.get("attrs", {}).items():
            if isinstance(v, bool):
                lines.append(f'    <att name="{xml_esc(k)}" value="{str(v).lower()}" type="boolean"/>')
            elif isinstance(v, int):
                lines.append(f'    <att name="{xml_esc(k)}" value="{v}" type="integer"/>')
            elif isinstance(v, float):
                lines.append(f'    <att name="{xml_esc(k)}" value="{v:.4f}" type="real"/>')
            else:
                lines.append(f'    <att name="{xml_esc(k)}" value="{xml_esc(str(v))}" type="string"/>')
        lines.append(
            f'    <graphics type="{nd.get("shape", "ELLIPSE")}" '
            f'x="{nd["x"]:.2f}" y="{nd["y"]:.2f}" '
            f'w="{nd.get("w", 12.0):.2f}" h="{nd.get("h", 12.0):.2f}" '
            f'fill="{nd.get("fill", "#CBD5E1")}" '
            f'width="{nd.get("border_w", 1.0):.2f}" '
            f'outline="{nd.get("outline", "#64748B")}">'
        )
        lines.append(f'      <att name="NODE_TRANSPARENCY" value="{alpha_255}" type="integer"/>')
        lines.append(f'      <att name="NODE_BORDER_TRANSPARENCY" value="{alpha_255}" type="integer"/>')
        lines.append(f'      <att name="NODE_LABEL" value="{xml_esc(slbl)}" type="string"/>')
        lines.append(f'      <att name="NODE_LABEL_COLOR" value="{lcol}" type="string"/>')
        lines.append(f'      <att name="NODE_LABEL_FONT_SIZE" value="{fsz}" type="integer"/>')
        # LOCAL ADDITION (not present in fullgrn_data_and_xgmml.py): pass through
        # extra Cytoscape visual properties, e.g. NODE_LABEL_POSITION, so a label
        # can be placed outside its node shape.
        for _k, _v in nd.get('graphics_atts', {}).items():
            lines.append(f'      <att name="{xml_esc(_k)}" value="{xml_esc(str(_v))}" type="string"/>')
        lines.append('    </graphics>')
        lines.append('  </node>')

    for eidx, ed in enumerate(edges_spec, start=1):
        su = id_map.get(ed["source"])
        sv = id_map.get(ed["target"])
        if not su or not sv:
            continue
        inter = ed.get("interaction", "regulates")
        elbl = f'{ed["source"]} ({inter}) {ed["target"]}'
        alpha_255 = max(5, min(255, int(round(ed.get("opacity", 0.5) * 255))))
        ltype = "EQUAL_DASH" if ed.get("dashed", False) else "SOLID"
        arrow = "DELTA" if ed.get("target_arrow", True) else "NONE"
        stroke = ed.get("stroke", "#94A3B8")
        ew = ed.get("width", 1.0)

        lines.append(f'  <edge id="e{eidx}" label="{xml_esc(elbl)}" source="{su}" target="{sv}" cy:directed="1">')
        lines.append(f'    <att name="interaction" value="{xml_esc(inter)}" type="string"/>')
        for k, v in ed.get("attrs", {}).items():
            if isinstance(v, int):
                lines.append(f'    <att name="{xml_esc(k)}" value="{v}" type="integer"/>')
            elif isinstance(v, float):
                lines.append(f'    <att name="{xml_esc(k)}" value="{v:.4f}" type="real"/>')
            else:
                lines.append(f'    <att name="{xml_esc(k)}" value="{xml_esc(str(v))}" type="string"/>')
        lines.append(f'    <graphics width="{ew:.2f}" fill="{stroke}">')
        lines.append(f'      <att name="EDGE_TRANSPARENCY" value="{alpha_255}" type="integer"/>')
        lines.append(f'      <att name="EDGE_LINE_TYPE" value="{ltype}" type="string"/>')
        lines.append(f'      <att name="EDGE_TARGET_ARROW_SHAPE" value="{arrow}" type="string"/>')
        lines.append(f'      <att name="EDGE_STROKE_UNSELECTED_PAINT" value="{stroke}" type="string"/>')
        lines.append('    </graphics>')
        lines.append('  </edge>')

    lines.append('</graph>')
    filepath.parent.mkdir(parents=True, exist_ok=True)
    filepath.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return filepath

# ----------------------------------------------------------------- inputs ---
IN = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(".")
OUT = Path("results/fig2")
FFL_CSV = IN / "ffl_triad_significance.csv"
BIF_CSV = IN / "bifan_pair_significance.csv"
FUN_CSV = IN / "function.with.classification.csv"
XGMML_IN = IN / "Cytoscape_Fig2_FullGRN_3Tier_Architecture_Spotlight.xgmml"

# condition key (as spelled in the Step 4 tables) -> accepted edge-file names
GRN_FILES = {
    "amikacin":      ["ami.grn.csv", "amikacin.grn.csv"],
    "ciprofloxacin": ["cip.grn.csv", "ciprofloxacin.grn.csv"],
    "doripenem":     ["dori.grn.csv", "doripenem.grn.csv"],
    "meropenem":     ["mero.grn.csv", "meropenem.grn.csv"],
    "polymyxin":     ["poly.grn.csv", "polymyxin.grn.csv"],
    "tigecycline":   ["tige.grn.csv", "tigecycline.grn.csv"],
    "colistin_res":  ["colis.resis.grn.csv", "colistin_res.grn.csv", "colis.grn.csv"],
}

COND_DISPLAY = {
    "amikacin": "Amikacin", "ciprofloxacin": "Ciprofloxacin", "doripenem": "Doripenem",
    "meropenem": "Meropenem", "polymyxin": "Polymyxin B", "tigecycline": "Tigecycline",
    "colistin_res": "Colistin Resistance",
}

# Used ONLY in XGMML fallback mode, where per-condition membership is unknown.
HUB_ACTIVE_BADGES_FALLBACK = {
    "beti": "7/7", "dmlr": "7/7", "arsr2": "5/7",
    "fadr": "7/7", "lrp": "7/7", "pgrr": "6/7", "rhas": "6/7",
    "rcsc": "7/7", "paax": "5/7", "yybr": "3/7", "rpoh": "2/7",
}


# =========================================================================
# STEP 1 - TOPOLOGY
# =========================================================================

def find_grn_files(folder):
    """Return {condition_key: Path} for every per-condition edge file present."""
    found = {}
    for cond, names in GRN_FILES.items():
        for nm in names:
            p = folder / nm
            if p.exists():
                found[cond] = p
                break
    return found


def load_pan_grn_from_csvs(files):
    """Seven (or fewer) `source,target` edge files -> pan-union topology.

    Returns
      all_nodes    set of genes
      pan_conds    {(u, v): set of condition keys the edge occurs in}
      cond_sources {condition: set of genes acting as a regulator there}
    """
    pan_conds = defaultdict(set)
    cond_sources = defaultdict(set)
    all_nodes = set()
    for cond, path in sorted(files.items()):
        n = 0
        with open(path, newline="", encoding="utf-8-sig") as fh:
            rdr = csv.DictReader(fh)
            cols = {c.strip().lower(): c for c in (rdr.fieldnames or [])}
            if "source" not in cols or "target" not in cols:
                raise SystemExit(f"{path.name}: expected a 'source,target' header, "
                                 f"found {rdr.fieldnames}")
            for r in rdr:
                u = (r[cols["source"]] or "").strip().lower()
                v = (r[cols["target"]] or "").strip().lower()
                if not u or not v:
                    continue
                pan_conds[(u, v)].add(cond)
                cond_sources[cond].add(u)
                all_nodes.update((u, v))
                n += 1
        print(f"  {path.name:<24s} {cond:<14s} {n:>5d} edges")
    return all_nodes, dict(pan_conds), dict(cond_sources)


def load_pan_grn_from_xgmml(path):
    """Fallback: topology from the published figure. Per-edge condition COUNT
    only, so pan_conds values are None and cond_sources is empty."""
    root = ET.parse(path).getroot()
    nid2gene, counts = {}, {}
    for el in root.iter():
        if el.tag.endswith("node"):
            atts = {c.attrib["name"]: c.attrib.get("value", "")
                    for c in el if c.tag.endswith("att")}
            nid2gene[el.attrib["id"]] = atts["gene_id"]
    for el in root.iter():
        if el.tag.endswith("edge"):
            atts = {c.attrib["name"]: c.attrib.get("value", "")
                    for c in el if c.tag.endswith("att")}
            counts[(nid2gene[el.attrib["source"]], nid2gene[el.attrib["target"]])] = \
                int(atts.get("n_conditions", 1))
    return set(nid2gene.values()), counts, {}


# =========================================================================
# STEP 2-3 - MOTIFS AND CURATED CLASSES
# =========================================================================

UNKNOWN = "Unknown / uncharacterized"

# --- labelling -------------------------------------------------------------
# Which validated effectors on the outer rim carry a visible label.
#   "mmgc_only"       only mmgC
#   "recurrence3plus" every effector found in >= 3 conditions
#   "sentinels"       recurrence >= 3 plus the curated FIG2_SENTINEL_GENES
#   "none"            no rim labels at all
# --- hub tier classification ----------------------------------------------
# "computed": derive each Master Hub's tier from the arrows of the
#             passes_significance FFL triads (see compute_hub_tiers).
# "fixed":    use the hardcoded HUB_TIERS lookup from the original code.
TIER_SOURCE = "computed"

# Which arrows of a triad X->Y->Z count. "all" uses X->Y, Y->Z and X->Z;
# "xy_only" uses just X->Y. Both give the same 11 tiers on this dataset.
FFL_ARROWS = "all"

# Place each hub on the ring of its COMPUTED tier (hubs within a tier are
# spread evenly, preserving their published angular order). False keeps the
# published radii, which can leave a hub's colour disagreeing with its ring.
TIER_DRIVES_RADIUS = True

RIM_LABELS = "mmgc_only"

# Write Fig2_render.png / .svg as well as the XGMML and CSVs.
# Needs matplotlib; set False to skip (nothing else requires it).
RENDER_IMAGE = True
RIM_LABEL_GENES = {"mmgc"}

# Cytoscape label placement: "<node anchor>,<label anchor>,<justify>,<dx>,<dy>".
# "S,N,c,0.00,6.00" hangs the text below the shape instead of printing over it.
LABEL_OUTSIDE_HUB = "S,N,c,0.00,8.00"
LABEL_OUTSIDE_RIM = "S,N,c,0.00,5.00"

CLASS_ORDER = [
    "Cell envelope, membrane & cell-wall maintenance",
    "Metabolism & energy production",
    "Transport, efflux",
    "Antibiotic/metal/toxin resistance & protection",
    "Stress response & detoxification",
    "Information processing: DNA/RNA/protein synthesis & repair",
    "Nutrient acquisition & iron/heme metabolism",
    "Virulence, adhesion, motility & secretion",
    "Regulatory & signaling functions",
    UNKNOWN,
]

CLASS_STYLE = {
    "Cell envelope, membrane & cell-wall maintenance":
        {"short": "Cell Envelope & Cell Wall", "hex": "#0D9488", "light": "#CCFBF1"},
    "Metabolism & energy production":
        {"short": "Metabolism & Energy", "hex": "#D97706", "light": "#FEF3C7"},
    "Transport, efflux":
        {"short": "Transport & Efflux", "hex": "#DC2626", "light": "#FEE2E2"},
    "Antibiotic/metal/toxin resistance & protection":
        {"short": "Antibiotic/Metal Resistance", "hex": "#BE185D", "light": "#FCE7F3"},
    "Stress response & detoxification":
        {"short": "Stress & Detoxification", "hex": "#EA580C", "light": "#FFEDD5"},
    "Information processing: DNA/RNA/protein synthesis & repair":
        {"short": "DNA/RNA/Protein Processing", "hex": "#7C3AED", "light": "#EDE9FE"},
    "Nutrient acquisition & iron/heme metabolism":
        {"short": "Nutrient & Iron Acquisition", "hex": "#9333EA", "light": "#F3E8FF"},
    "Virulence, adhesion, motility & secretion":
        {"short": "Virulence, Motility & Secretion", "hex": "#059669", "light": "#D1FAE5"},
    "Regulatory & signaling functions":
        {"short": "Regulatory & Signaling", "hex": "#2563EB", "light": "#DBEAFE"},
    UNKNOWN:
        {"short": "Unknown / Uncharacterized", "hex": "#475569", "light": "#E2E8F0"},
}

FIG2_SENTINEL_GENES = {
    "lpxa", "kdsa", "olsb", "uppp", "deda", "ftsa", "lola", "lolb", "mltb", "fabg", "acpp",
    "mmgc", "fadh", "catd", "pcaj", "paaf",
    "mexb", "bepc", "mdta", "mepa", "tolc", "oprd", "bla-oxa-256",
    "recn", "sbcc", "comc", "smf-1", "xerd",
    "frr", "iles", "secd", "secf", "rtcb", "rpla",
    "hmp", "ahpc", "grxd", "mscl", "kdpb", "dsba",
    "cyda", "cydb", "menf",
    "baua", "btub", "dhba", "fhua",
    "cpda", "mcpb", "tfpa", "pilt", "pgab",
    "naip", "meth", "pyrd",
}


def load_motifs(all_sources):
    """Significant Master-Hub FFLs and pure bi-fans -> motif edges, partner TFs,
    and the per-condition recurrence of every validated B_Combined effector."""
    hub_motif_edges, partner_tfs = set(), set()
    b_cond = defaultdict(set)
    n_ffl = n_bif = 0

    with open(FFL_CSV, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r["passes_significance"] != "True":
                continue
            n_ffl += 1
            x, y, z = (r["X"].strip().lower(), r["Y"].strip().lower(), r["Z"].strip().lower())
            hub_motif_edges |= {(x, y), (y, z), (x, z)}
            for g in (x, y):
                if g not in HUB_TIERS and g in all_sources:
                    partner_tfs.add(g)
            if z not in all_sources:
                b_cond[z].add(r["condition"].strip().lower())

    with open(BIF_CSV, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r["passes_significance"] != "True":
                continue
            n_bif += 1
            c = r["condition"].strip().lower()
            xs = [r["X1"].strip().lower(), r["X2"].strip().lower()]
            zs = [r["Z1"].strip().lower(), r["Z2"].strip().lower()]
            for u in xs:
                if u not in HUB_TIERS and u in all_sources:
                    partner_tfs.add(u)
                for v in zs:
                    hub_motif_edges.add((u, v))
            for v in zs:
                if v not in all_sources:
                    b_cond[v].add(c)

    return hub_motif_edges, partner_tfs, b_cond, n_ffl, n_bif


def load_classes():
    cls = {}
    with open(FUN_CSV, newline="", encoding="utf-8-sig") as fh:
        rdr = csv.DictReader(fh)
        rdr.fieldnames = [f.strip() for f in rdr.fieldnames]
        for r in rdr:
            g = (r.get("Gene") or "").strip().lower()
            k = (r.get("Your class name") or "").strip()
            if g:
                cls[g] = k if k in CLASS_STYLE else UNKNOWN
    return cls



def compute_hub_tiers(all_sources):
    """Hub tier from the arrows of the significant FFL triads.

      k_in_hub = arrows received from OTHER HUBS
      k_out_tf = arrows sent to other regulators (hub or non-hub TF)

      Tier 1  k_in_hub == 0                  no hub regulates it; upstream
      Tier 2  k_in_hub > 0 and k_out_tf > 0  hub-controlled, regulates regulators
      Tier 3  k_in_hub > 0 and k_out_tf == 0 hub-controlled, end-point

    Returns (tier_map, diagnostics) where tier_map mirrors HUB_TIERS:
    {hub: (tier_name, tier_number, display_label)}.
    """
    arrows = set()
    with open(FFL_CSV, newline="", encoding="utf-8-sig") as fh:
        for r in csv.DictReader(fh):
            if r["passes_significance"] != "True":
                continue
            x, y, z = (r["X"].strip().lower(), r["Y"].strip().lower(), r["Z"].strip().lower())
            arrows |= {(x, y)} if FFL_ARROWS == "xy_only" else {(x, y), (y, z), (x, z)}
    arrows = {(u, v) for u, v in arrows if u != v}

    names = {1: "Top_Initiator", 2: "Mid_Integrator", 3: "Bottom_Executor"}
    tier_map, diag = {}, []
    for h in HUB_TIERS:
        k_in_hub = sum(1 for u, v in arrows if v == h and u in HUB_TIERS)
        k_out_tf = sum(1 for u, v in arrows if u == h and v in all_sources)
        tier = 1 if k_in_hub == 0 else (2 if k_out_tf > 0 else 3)
        tier_map[h] = (names[tier], tier, HUB_TIERS[h][2])
        diag.append((HUB_TIERS[h][2], k_in_hub, k_out_tf, tier, HUB_TIERS[h][1]))
    return tier_map, diag


# =========================================================================
# STEP 4 - GEOMETRY
# =========================================================================

CX, CY = 1130.0, 1140.0
R_T1, R_T2, R_T3 = 152.0, 298.0, 442.0
R_SPOT = 632.0
R_BG_0, R_BG_2 = 678.0, 736.0

# Hub angles and radii as published (build_fig2_3tier_fullgrn).
HUB_ANGLES_DEG = {
    "dmlr": -90.0, "beti": -205.0, "arsr2": 25.0,
    "fadr": -128.0, "lrp": -38.0, "pgrr": 145.0, "rhas": 52.0,
    "rcsc": -90.0, "paax": -8.0, "yybr": 95.0, "rpoh": -172.0,
}
HUB_RADII = {
    "dmlr": R_T1, "beti": R_T1, "arsr2": R_T1,
    "fadr": R_T2, "lrp": R_T2, "pgrr": R_T2, "rhas": R_T2,
    "rcsc": R_T3, "paax": R_T3, "yybr": R_T3, "rpoh": R_T3,
}


def place_ring_tfs(coords, tf_list, radius, hubs_on_ring, min_hub_sep_deg):
    """Spread local TFs around a ring, keeping an angular buffer around the
    Master Hubs sitting on that same ring (as published)."""
    if not tf_list:
        return
    avail = []
    step = 360.0 / max(180, len(tf_list) * 3)
    a = -180.0
    while a < 180.0:
        if all(abs((a - HUB_ANGLES_DEG[h] + 180.0) % 360.0 - 180.0) >= min_hub_sep_deg
               for h in hubs_on_ring):
            avail.append(a)
        a += step
    for idx, g in enumerate(tf_list):
        pick = avail[int(idx * len(avail) / len(tf_list)) % len(avail)]
        rad = math.radians(pick)
        coords[g] = (CX + radius * math.cos(rad), CY + radius * math.sin(rad))


def build_layout(all_nodes, all_sources, pan_edges, partner_tfs, b_cond, gene_class,
                 hub_tier_map):
    """Full concentric layout: hubs at fixed angles, local TFs on the three
    inner rings, validated effectors in their 10 class sectors, unvalidated
    leaves as an unsectored halo."""
    coords = {}
    tier_radius = {1: R_T1, 2: R_T2, 3: R_T3}
    if TIER_DRIVES_RADIUS:
        # Spread each tier's hubs evenly on its ring, keeping published order.
        for tier, rr in tier_radius.items():
            members = sorted([h for h in HUB_ANGLES_DEG if h in all_nodes
                              and hub_tier_map[h][1] == tier],
                             key=lambda g: HUB_ANGLES_DEG[g])
            for i, h in enumerate(members):
                ang = -90.0 + 360.0 * i / max(1, len(members))
                coords[h] = (CX + rr * math.cos(math.radians(ang)),
                             CY + rr * math.sin(math.radians(ang)))
    else:
        for h, ang in HUB_ANGLES_DEG.items():
            if h in all_nodes:
                rr = HUB_RADII[h]
                coords[h] = (CX + rr * math.cos(math.radians(ang)),
                             CY + rr * math.sin(math.radians(ang)))

    # Non-hub TF tier from TF->TF in/out degree over the pan-union (as published)
    indeg, outdeg = defaultdict(int), defaultdict(int)
    for u, v in pan_edges:
        if u in all_sources and v in all_sources and u != v:
            outdeg[u] += 1
            indeg[v] += 1
    t1, t2, t3 = [], [], []
    for g in sorted(g for g in all_sources if g not in HUB_TIERS):
        if indeg[g] == 0:
            t1.append(g)
        elif outdeg[g] >= 2 or g in partner_tfs:
            t2.append(g)
        else:
            t3.append(g)
    place_ring_tfs(coords, t1, R_T1, ["dmlr", "beti", "arsr2"], 18.0)
    place_ring_tfs(coords, t2, R_T2, ["fadr", "lrp", "pgrr", "rhas"], 15.0)
    place_ring_tfs(coords, t3, R_T3, ["rcsc", "paax", "yybr", "rpoh"], 13.0)

    # 10-class horizon, sector width weighted by validated-effector count
    by_class = defaultdict(list)
    for g in b_cond:
        by_class[gene_class[g]].append(g)
    gap_deg = 2.4
    usable = 360.0 - gap_deg * len(CLASS_ORDER)
    w = {k: (95.0 if k == UNKNOWN else min(130.0, len(by_class[k]) * 1.35 + 24.0))
         for k in CLASS_ORDER}
    sum_w = sum(w.values())

    cur = -112.0
    sector_spans = {}
    for k in CLASS_ORDER:
        span = usable * w[k] / sum_w
        a0, a1 = cur, cur + span
        sector_spans[k] = (a0, a1, len(by_class[k]))
        for i, g in enumerate(sorted(by_class[k], key=lambda x: (-len(b_cond[x]), x))):
            frac = (i + 0.5) / max(1, len(by_class[k]))
            ang = math.radians(a0 + 1.0 + frac * max(0.5, span - 2.0))
            rec = len(b_cond[g])
            rr = R_SPOT - (12.0 if rec >= 3 else (5.0 if rec == 2 else 0.0))
            coords[g] = (CX + rr * math.cos(ang), CY + rr * math.sin(ang))
        cur = a1 + gap_deg

    # Unvalidated background leaves: uniform halo over 8 staggered sub-tracks
    ghosts = sorted(all_nodes - all_sources - set(b_cond))
    n_sub = 8
    for i, g in enumerate(ghosts):
        ang = math.radians(-180.0 + 360.0 * (i + 0.5) / len(ghosts))
        rr = R_BG_0 + (i % n_sub) * ((R_BG_2 - R_BG_0) / (n_sub - 1))
        coords[g] = (CX + rr * math.cos(ang), CY + rr * math.sin(ang))

    return coords, sector_spans, by_class, ghosts


# =========================================================================
# STEP 5-7 - STYLE, WRITE
# =========================================================================

def main():
    # ---- Step 1: topology -------------------------------------------------
    grn_files = find_grn_files(IN)
    if grn_files:
        mode = "grn"
        print(f"Topology from {len(grn_files)} per-condition GRN edge file(s):")
        all_nodes, pan_conds, cond_sources = load_pan_grn_from_csvs(grn_files)
        pan_n = {e: len(c) for e, c in pan_conds.items()}
        n_conditions = len(grn_files)
    elif XGMML_IN.exists():
        mode = "xgmml"
        print(f"No GRN edge files found in {IN.resolve()} - falling back to "
              f"{XGMML_IN.name} (per-edge condition COUNT only).")
        all_nodes, pan_n, cond_sources = load_pan_grn_from_xgmml(XGMML_IN)
        pan_conds = {}
        n_conditions = 7
    else:
        raise SystemExit(
            f"No topology in {IN.resolve()}.  Provide either the per-condition "
            f"edge files ({', '.join(v[0] for v in GRN_FILES.values())}) or "
            f"{XGMML_IN.name}.")

    all_sources = {u for u, _ in pan_n}

    # ---- Steps 2-3: motifs and classes ------------------------------------
    hub_motif_edges, partner_tfs, b_cond, n_ffl, n_bif = load_motifs(all_sources)
    curated = load_classes()
    gene_class = {g: curated.get(g, UNKNOWN) for g in b_cond}

    # Hub activity badges: computed when per-condition membership is known
    if mode == "grn":
        hub_badges = {h: f"{sum(1 for c in cond_sources if h in cond_sources[c])}"
                         f"/{n_conditions}" for h in HUB_TIERS}
    else:
        hub_badges = dict(HUB_ACTIVE_BADGES_FALLBACK)

    # ---- Hub tiers ---------------------------------------------------------
    if TIER_SOURCE == "computed":
        hub_tier_map, tier_diag = compute_hub_tiers(all_sources)
    else:
        hub_tier_map, tier_diag = dict(HUB_TIERS), []

    # ---- Step 4: layout ----------------------------------------------------
    coords, sector_spans, by_class, ghosts = build_layout(
        all_nodes, all_sources, pan_n.keys(), partner_tfs, b_cond, gene_class,
        hub_tier_map)

    val_hub_leaf = {(u, v) for (u, v) in hub_motif_edges
                    if u in HUB_TIERS and v in b_cond}

    # ---- Step 5: node styling ---------------------------------------------
    xg_nodes, node_rows = [], []
    for g in sorted(all_nodes):
        x, y = coords[g]
        lbl = fmt_gene(g)
        rec = len(b_cond[g]) if g in b_cond else 0
        if g in HUB_TIERS:
            tname, tnum, hlbl = hub_tier_map[g]
            fill, outl, shp, sz = (("#0F172A", "#38BDF8", "HEXAGON", 52.0) if tnum == 1 else
                                   ("#1D4ED8", "#93C5FD", "HEXAGON", 48.0) if tnum == 2 else
                                   ("#0F766E", "#5EEAD4", "OCTAGON", 48.0))
            slbl = f"{hlbl} [{hub_badges[g]}]"
            ncls = f"Master_Hub_Tier{tnum}_{tname}"
            klass, sector = "Master_Hub_Command_Core", f"00_Master_Hub_Tier{tnum}"
            opacity, bw, fsz, lblcol = 1.0, 3.2, 16, "#0F172A"
            gatts = {"NODE_LABEL_POSITION": LABEL_OUTSIDE_HUB}
        elif g in all_sources:
            is_ptf = g in partner_tfs
            fill = "#6366F1" if is_ptf else "#94A3B8"
            outl = "#4F46E5" if is_ptf else "#64748B"
            shp = "DIAMOND" if is_ptf else "ELLIPSE"
            sz, slbl = (9.0 if is_ptf else 6.0), ""
            ncls = "Partner_TF_Motif_Relay" if is_ptf else "Background_Local_TF"
            klass, sector = "Local_Regulator", "01_Transcriptional_Regulator"
            opacity, bw, fsz, lblcol = (0.42 if is_ptf else 0.22), 0.8, 9, "#475569"
            gatts = {}
        elif g in b_cond:
            klass = gene_class[g]
            fill = CLASS_STYLE[klass]["hex"]
            outl = "#FFFFFF" if rec == 1 else "#0F172A"
            shp = "DIAMOND" if rec >= 3 else "ELLIPSE"
            sz = 24.0 if rec >= 5 else (17.0 if rec >= 3 else (13.0 if rec == 2 else 9.5))
            if RIM_LABELS == "mmgc_only":
                show = g in RIM_LABEL_GENES
            elif RIM_LABELS == "recurrence3plus":
                show = rec >= 3
            elif RIM_LABELS == "sentinels":
                show = rec >= 3 or g in FIG2_SENTINEL_GENES
            else:
                show = False
            slbl = f"{lbl} ({rec}/{n_conditions})" if show else ""
            gatts = {"NODE_LABEL_POSITION": LABEL_OUTSIDE_RIM} if show else {}
            ncls = ("Spotlight_MultiDrug_Effector_3plus" if rec >= 3 else
                    "Spotlight_TwoDrug_Effector" if rec == 2 else "Spotlight_SingleDrug_Effector")
            sector = f"{CLASS_ORDER.index(klass) + 2:02d}_{CLASS_STYLE[klass]['short']}"
            opacity, bw = 0.96, (2.0 if rec >= 3 else 1.0)
            fsz, lblcol = (12 if rec >= 3 else 10), "#0F172A"
        else:
            fill, outl, shp, sz, slbl = "#CBD5E1", "#94A3B8", "ELLIPSE", 4.5, ""
            ncls = "Ghosted_Unvalidated_GRN_Leaf"
            klass, sector = "Unclassified_background", "99_Unclassified_Background"
            opacity, bw, fsz, lblcol = 0.20, 0.3, 8, "#94A3B8"
            gatts = {}

        xg_nodes.append({
            "id": g, "name": lbl, "spotlight_label": slbl, "x": x, "y": y,
            "w": sz, "h": sz, "shape": shp, "fill": fill, "outline": outl,
            "border_w": bw, "opacity": opacity, "label_color": lblcol, "font_size": fsz,
            "graphics_atts": gatts,
            "attrs": {"role": ncls, "functional_class": klass, "recurrence": rec},
        })
        node_rows.append([g, lbl, slbl, ncls, klass, sector, rec, fill, outl, sz,
                          int(round(opacity * 255)), round(x, 2), round(y, 2)])

    # ---- Step 5: edge styling ---------------------------------------------
    xg_edges, edge_rows = [], []
    for (u, v), n_c in sorted(pan_n.items()):
        cond_str = "; ".join(COND_DISPLAY.get(c, c) for c in sorted(pan_conds.get((u, v), [])))
        if u in HUB_TIERS and v in HUB_TIERS:
            if (u, v) == ("dmlr", "fadr"):
                ecls, col, ew, op = "Constitutive_3Tier_Backbone", "#DC2626", 6.2, 0.98
            elif n_c >= 2:
                ecls, col, ew, op = "Recurrent_3Tier_Hub_Link", "#D97706", 3.4, 0.92
            else:
                ecls, col, ew, op = "SingleCond_3Tier_Hub_Link", "#0284C7", 1.9, 0.78
            inter, arrow = "3Tier_Hub_Backbone", True
        elif (u, v) in val_hub_leaf:
            rec_v = len(b_cond[v])
            col = CLASS_STYLE[gene_class[v]]["hex"]
            ecls = ("Validated_Hub_to_MultiDrug_Effector" if rec_v >= 3
                    else "Validated_Hub_to_Class_Effector")
            ew, op = (2.1, 0.75) if rec_v >= 3 else (1.05, 0.34)
            inter, arrow = "Validated_Hub_Effector", False
        else:
            ecls, col, ew, op = "Ghosted_FullGRN_Background_Edge", "#94A3B8", 0.45, 0.08
            inter, arrow = "Ghosted_GRN_Edge", False

        xg_edges.append({"source": u, "target": v, "interaction": inter, "stroke": col,
                         "width": ew, "opacity": op, "target_arrow": arrow,
                         "attrs": {"n_conditions": n_c, "conditions": cond_str,
                                   "tier_class": ecls}})
        edge_rows.append([u, v, inter, n_c, cond_str, ecls, col, ew, int(round(op * 255))])

    # ---- Steps 6-7: write --------------------------------------------------
    OUT.mkdir(parents=True, exist_ok=True)
    write_xgmml(
        OUT / "Cytoscape_Fig2_FullGRN_3Tier_CuratedClasses.xgmml",
        f"Fig2_Full_PanAntibiotic_GRN_3Tier_CuratedClasses (N={len(all_nodes)}, E={len(pan_n)})",
        xg_nodes, xg_edges, bg_color="#F8FAFC", center_x=CX, center_y=CY, zoom=0.55)

    with open(OUT / "Fig2_Nodes_CuratedClasses.csv", "w", newline="", encoding="utf-8") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["id", "gene_symbol", "fig2_spotlight_label", "fig2_node_class",
                     "functional_class", "horizon_sector", "effector_recurrence",
                     "fill_hex", "outline_hex", "node_size_px", "node_opacity_255",
                     "preset_x", "preset_y"])
        w_.writerows(node_rows)

    with open(OUT / "Fig2_Edges_CuratedClasses.csv", "w", newline="", encoding="utf-8") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["source", "target", "interaction", "n_conditions", "conditions_list",
                     "fig2_edge_class", "edge_color_hex", "edge_width_pt", "edge_opacity_255"])
        w_.writerows(edge_rows)

    with open(OUT / "Fig2_Class_Legend.csv", "w", newline="", encoding="utf-8") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["sector_index", "functional_class", "short_label", "hex", "light_hex",
                     "n_validated_effectors", "sector_start_deg", "sector_end_deg"])
        for i, k in enumerate(CLASS_ORDER, start=2):
            a0, a1, n = sector_spans[k]
            w_.writerow([f"{i:02d}", k, CLASS_STYLE[k]["short"], CLASS_STYLE[k]["hex"],
                         CLASS_STYLE[k]["light"], n, round(a0, 2), round(a1, 2)])

    # ---- report ------------------------------------------------------------
    print(f"\ntopology mode        : {mode}  ({n_conditions} condition(s))")
    print(f"significant motifs   : {n_ffl} FFL triads | {n_bif} pure bi-fans")
    print(f"nodes {len(xg_nodes)} | edges {len(xg_edges)} | effectors {len(b_cond)} | "
          f"partner TFs {len(partner_tfs)} | ghosted leaves {len(ghosts)}")
    print(f"validated hub->effector edges present in the GRN: "
          f"{len(val_hub_leaf & set(pan_n))} of {len(val_hub_leaf)}")
    if tier_diag:
        print(f"hub tiers (computed from {len(set(h for h in HUB_TIERS))} hubs, "
              f"FFL_ARROWS={FFL_ARROWS!r}):")
        print(f"  {'hub':<8s}{'k_in_hub':>9s}{'k_out_tf':>9s}{'tier':>6s}   hardcoded")
        for lbl, kin, kout, tier, old in sorted(tier_diag, key=lambda r: (r[3], r[0])):
            flag = "" if tier == old else f"  <- was T{old}"
            print(f"  {lbl:<8s}{kin:>9d}{kout:>9d}{'T'+str(tier):>6s}{flag}")
    print(f"hub activity badges  : "
          + ", ".join(f"{HUB_TIERS[h][2]} {hub_badges[h]}" for h in HUB_TIERS))
    for k in CLASS_ORDER:
        print(f"  {CLASS_STYLE[k]['short']:<32s} {len(by_class[k]):>4d}")

    if RENDER_IMAGE:
        try:
            render(OUT / "Cytoscape_Fig2_FullGRN_3Tier_CuratedClasses.xgmml")
        except ImportError:
            print("\nmatplotlib not installed - XGMML and CSVs written, image skipped.")

    if RENDER_POSTER:
        try:
            render_poster(
                xg_nodes=xg_nodes, xg_edges=xg_edges, sector_spans=sector_spans,
                by_class=by_class, b_cond=b_cond, gene_class=gene_class,
                hub_tier_map=hub_tier_map, hub_badges=hub_badges,
                partner_tfs=partner_tfs, all_sources=all_sources, all_nodes=all_nodes,
                pan_n=pan_n, ghosts=ghosts, n_conditions=n_conditions, out_dir=OUT)
        except ImportError:
            print("\nmatplotlib not installed - poster image skipped.")


# =========================================================================
# STEP 8 - RENDER (merged from render_fig2_png.py)
# =========================================================================

DPI = 400
FIGSIZE = (13.0, 13.0)
BG = "#F8FAFC"

# Marker per node shape written by the builder
SHAPE_MARKER = {"HEXAGON": "H", "OCTAGON": "8", "DIAMOND": "D",
                "ROUND_RECTANGLE": "s", "ELLIPSE": "o"}

# Draw order: background first, spotlight last
ROLE_Z = {
    # Nodes occupy a z-band strictly ABOVE every edge, so no line is ever
    # painted across a node shape.
    "Ghosted_Unvalidated_GRN_Leaf": 10,
    "Background_Local_TF": 11,
    "Partner_TF_Motif_Relay": 12,
    "Spotlight_SingleDrug_Effector": 13,
    "Spotlight_TwoDrug_Effector": 14,
    "Spotlight_MultiDrug_Effector_3plus": 15,
    "Master_Hub_Tier3_Bottom_Executor": 16,
    "Master_Hub_Tier2_Mid_Integrator": 16,
    "Master_Hub_Tier1_Top_Initiator": 16,
}
EDGE_Z = {
    # Edges occupy z 1-6; see ROLE_Z for the node band above them.
    "Ghosted_FullGRN_Background_Edge": 1,
    "Validated_Hub_to_Class_Effector": 2,
    "Validated_Hub_to_MultiDrug_Effector": 3,
    "SingleCond_3Tier_Hub_Link": 4,
    "Recurrent_3Tier_Hub_Link": 5,
    "Constitutive_3Tier_Backbone": 6,
}


def parse_xgmml_for_render(path):
    root = ET.parse(path).getroot()
    nid, nodes = {}, []
    for el in root.iter():
        if not el.tag.endswith("node"):
            continue
        att = {c.attrib["name"]: c.attrib.get("value", "")
               for c in el if c.tag.endswith("att")}
        gr = next(c for c in el if c.tag.endswith("graphics"))
        gatt = {c.attrib["name"]: c.attrib.get("value", "")
                for c in gr if c.tag.endswith("att")}
        nid[el.attrib["id"]] = att["gene_id"]
        nodes.append({
            "id": att["gene_id"], "label": att.get("spotlight_label", ""),
            "role": att.get("role", ""), "klass": att.get("functional_class", ""),
            "x": float(gr.attrib["x"]), "y": float(gr.attrib["y"]),
            "w": float(gr.attrib["w"]), "shape": gr.attrib.get("type", "ELLIPSE"),
            "fill": gr.attrib["fill"], "outline": gr.attrib.get("outline", "#64748B"),
            "border_w": float(gr.attrib.get("width", 1.0)),
            "alpha": int(gatt.get("NODE_TRANSPARENCY", 255)) / 255.0,
            "font": int(gatt.get("NODE_LABEL_FONT_SIZE", 10)),
            "label_pos": gatt.get("NODE_LABEL_POSITION", ""),
        })
    edges = []
    for el in root.iter():
        if not el.tag.endswith("edge"):
            continue
        att = {c.attrib["name"]: c.attrib.get("value", "")
               for c in el if c.tag.endswith("att")}
        gr = next(c for c in el if c.tag.endswith("graphics"))
        gatt = {c.attrib["name"]: c.attrib.get("value", "")
                for c in gr if c.tag.endswith("att")}
        edges.append({
            "u": nid[el.attrib["source"]], "v": nid[el.attrib["target"]],
            "tier": att.get("tier_class", ""), "w": float(gr.attrib.get("width", 1.0)),
            "fill": gr.attrib["fill"],
            "alpha": int(gatt.get("EDGE_TRANSPARENCY", 128)) / 255.0,
        })
    return nodes, edges, root.attrib.get("label", "Figure 2")


def render(xgmml_path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.lines import Line2D

    nodes, edges, title = parse_xgmml_for_render(xgmml_path)
    pos = {n["id"]: (n["x"], -n["y"]) for n in nodes}   # XGMML y grows downward

    fig, ax = plt.subplots(figsize=FIGSIZE)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)

    # --- edges, grouped by tier so each group gets one collection ----------
    by_tier = OrderedDict()
    for e in edges:
        by_tier.setdefault(e["tier"], []).append(e)
    for tier, grp in sorted(by_tier.items(), key=lambda kv: EDGE_Z.get(kv[0], 5)):
        segs = [(pos[e["u"]], pos[e["v"]]) for e in grp]
        ax.add_collection(LineCollection(
            segs, colors=[e["fill"] for e in grp],
            linewidths=[e["w"] * 0.55 for e in grp],
            alpha=max(0.05, sum(e["alpha"] for e in grp) / len(grp)),
            zorder=EDGE_Z.get(tier, 3)))

    # --- nodes, grouped by (role, shape) ------------------------------------
    groups = OrderedDict()
    for n in nodes:
        groups.setdefault((n["role"], n["shape"]), []).append(n)
    for (role, shape), grp in sorted(groups.items(),
                                     key=lambda kv: ROLE_Z.get(kv[0][0], 16)):
        ax.scatter([pos[n["id"]][0] for n in grp], [pos[n["id"]][1] for n in grp],
                   s=[(n["w"] ** 1.45) / 4.0 for n in grp],
                   c=[n["fill"] for n in grp],
                   edgecolors=[n["outline"] for n in grp],
                   linewidths=[n["border_w"] * 0.45 for n in grp],
                   marker=SHAPE_MARKER.get(shape, "o"),
                   alpha=max(0.05, sum(n["alpha"] for n in grp) / len(grp)),
                   zorder=ROLE_Z.get(role, 16))

    # --- labels, hung outside the shape when the XGMML says so --------------
    pill = dict(boxstyle="round,pad=0.18", fc=BG, ec="none", alpha=0.85)
    for n in nodes:
        if not n["label"]:
            continue
        dy = -17 if n["label_pos"].startswith("S,N") else 0
        ax.annotate(n["label"], pos[n["id"]], textcoords="offset points",
                    xytext=(0, dy), ha="center",
                    va="top" if dy else "center",
                    fontsize=max(7, n["font"] * 0.55), style="italic",
                    fontweight="bold", color="#0F172A", zorder=20, bbox=pill)

    # --- legend built from the spotlight nodes themselves -------------------
    seen = {}
    for n in nodes:
        if n["role"].startswith("Spotlight") and n["klass"]:
            seen.setdefault(n["klass"], {"hex": n["fill"], "n": 0})["n"] += 1
    handles = [Line2D([], [], marker="o", ls="", ms=5, color=v["hex"],
                      label=f"{k}  ({v['n']})")
               for k, v in sorted(seen.items(), key=lambda kv: -kv[1]["n"])]
    leg = ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(-0.02, 1.0),
                    frameon=False, fontsize=6.5, title="Curated functional class",
                    title_fontsize=7.5)
    leg.get_title().set_fontweight("bold")
    for txt in leg.get_texts():
        txt.set_fontweight("bold")

    ax.set_aspect("equal")
    ax.axis("off")
    ax.margins(0.05)
    ax.set_title(title, fontsize=10, loc="left")

    out_png = xgmml_path.with_name("Fig2_render.png")
    out_svg = xgmml_path.with_name("Fig2_render.svg")
    fig.savefig(out_png, dpi=DPI, bbox_inches="tight", facecolor=BG)
    fig.savefig(out_svg, bbox_inches="tight", facecolor=BG)
    print(f"{len(nodes)} nodes, {len(edges)} edges")
    print(f"wrote {out_png}")
    print(f"wrote {out_svg}")


# =========================================================================
# STEP 9 - POSTER FIGURE
#   header band + 10 functional-class arcs (with the class name written on
#   each arc) + outside gene labels with leader lines + right-hand A/B/C panel.
#   Every number, name and label is computed from the data above.
# =========================================================================

RENDER_POSTER = True            # write Fig2_poster.png / .svg
POSTER_DPI = 200
POSTER_FIGSIZE = (23.0, 15.85)
POSTER_BG = "#FFFFFF"

TIER_WORD = "Layer"             # wording used for the hub hierarchy ("Tier" / "Layer")
SHOW_SUMMARY_BOX = False        # True brings back the yellow "No universal effector set" box
SHOW_CLASS_NAMES = True         # write each functional-class name on its outer arc

LABEL_MIN_REC = 2               # effector gets an outside label if found in >= this many conditions
R_ARC_IN, R_ARC_OUT = 742.0, 798.0   # radial band of the class arcs
R_LABEL = 1025.0                # x-distance of the label columns from the centre
LABEL_GAP = 34.0                # minimum vertical spacing between stacked labels

# Ghosted full-GRN background: faint, but visible (poster only)
GHOST_EDGE_ALPHA = 0.16
GHOST_EDGE_COLOR = "#94A3B8"
GHOST_EDGE_LW = 0.40
GHOST_NODE_ALPHA = 0.60

TIER_TITLE = {1: f"{TIER_WORD} 1 - Top Initiators (Shell 1)",
              2: f"{TIER_WORD} 2 - Mid Integrators (Shell 2)",
              3: f"{TIER_WORD} 3 - Bottom Executors (Shell 3)"}
NUM_WORD = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six",
            7: "seven", 8: "eight", 9: "nine", 10: "ten"}


def _spread_labels(targets, gap):
    """Positions as close as possible to `targets` (listed top -> bottom) while
    keeping every neighbouring pair at least `gap` apart and the order unchanged
    (isotonic regression / pool-adjacent-violators)."""
    z = [t + i * gap for i, t in enumerate(targets)]
    blocks = []
    for v in z:
        blocks.append([v, 1])
        while len(blocks) > 1 and blocks[-2][0] / blocks[-2][1] < blocks[-1][0] / blocks[-1][1]:
            s, c = blocks.pop()
            blocks[-1][0] += s
            blocks[-1][1] += c
    q = []
    for s, c in blocks:
        q += [s / c] * c
    return [q[i] - i * gap for i in range(len(targets))]


def _fit_text(fig, t, right_px, min_size=7.0):
    """Shrink a Text object until its right edge is inside `right_px`."""
    rend = fig.canvas.get_renderer()
    while t.get_window_extent(rend).x1 > right_px and t.get_fontsize() > min_size:
        t.set_fontsize(t.get_fontsize() - 0.4)
    return t


def _arc_text(fig, ax, text, cx, cy, r, phi_deg, arc_len, color,
              max_size=11.0, min_size=6.5):
    """Write `text` along a circle of radius r (data coords, y up), centred on the
    visual angle phi_deg, one character at a time so the text follows the arc.
    Font size is reduced so the string fits inside `arc_len` (data units)."""
    from matplotlib import font_manager as fm
    font = fm.get_font(fm.findfont(fm.FontProperties(weight="bold")))
    ax.figure.canvas.draw()                         # make transData final
    px_per_unit = abs(ax.transData.transform((1, 0))[0] - ax.transData.transform((0, 0))[0])
    unit_per_pt = (fig.dpi / 72.0) / px_per_unit

    def advances(sz):
        font.set_size(sz, 72)
        return [font.load_char(ord(c)).linearHoriAdvance / 65536.0 * unit_per_pt for c in text]

    adv = advances(10.0)
    size = max(min_size, min(max_size, 10.0 * 0.90 * arc_len / max(sum(adv), 1e-6)))
    adv = advances(size)
    total = sum(adv)
    top = math.sin(math.radians(phi_deg)) >= 0      # top half reads clockwise
    pos_c = 0.0
    for ch, w in zip(text, adv):
        mid = pos_c + w / 2.0
        pos_c += w
        if ch == " ":
            continue
        d = (total / 2.0 - mid) / r                 # radians from the centre line
        ang = math.radians(phi_deg) + (d if top else -d)
        rot = math.degrees(ang) - 90.0 if top else math.degrees(ang) + 90.0
        ax.text(cx + r * math.cos(ang), cy + r * math.sin(ang), ch, fontsize=size,
                fontweight="bold", color=color, rotation=rot, rotation_mode="anchor",
                ha="center", va="center", zorder=27)


def render_poster(*, xg_nodes, xg_edges, sector_spans, by_class, b_cond, gene_class,
                  hub_tier_map, hub_badges, partner_tfs, all_sources, all_nodes,
                  pan_n, ghosts, n_conditions, out_dir):
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection
    from matplotlib.colors import to_rgba
    from matplotlib.patches import Wedge, Rectangle, FancyBboxPatch

    # ------------------------------------------------------------ numbers --
    n_nodes, n_edges = len(all_nodes), len(pan_n)
    hubs = [h for h in HUB_TIERS if h in all_nodes]
    n_regs = len(all_sources)
    n_terminal = n_nodes - n_regs
    n_eff = len(b_cond)
    n_partner = len(partner_tfs)
    n_bg_tf = len([g for g in all_sources if g not in HUB_TIERS and g not in partner_tfs])

    def tier_counts(cls):
        return [e["attrs"]["n_conditions"] for e in xg_edges if e["attrs"]["tier_class"] == cls]

    spine_bb = tier_counts("Constitutive_3Tier_Backbone")
    spine_rec = tier_counts("Recurrent_3Tier_Hub_Link")
    spine_one = tier_counts("SingleCond_3Tier_Hub_Link")
    n_spine = len(spine_bb) + len(spine_rec) + len(spine_one)
    n_ghost_edges = len(tier_counts("Ghosted_FullGRN_Background_Edge"))

    def rng(v):
        if not v:
            return "n/a"
        lo, hi = min(v), max(v)
        return f"{lo}/{n_conditions}" if lo == hi else f"{lo}-{hi}/{n_conditions}"

    rec_of = {g: len(c) for g, c in b_cond.items()}
    k_all = sum(1 for r in rec_of.values() if r == n_conditions)
    k_one = sum(1 for r in rec_of.values() if r == 1)
    partial = {g: r for g, r in rec_of.items() if 2 <= r < n_conditions}
    pct = lambda k: (100.0 * k / n_eff) if n_eff else 0.0
    lead = sorted((g for g, r in partial.items() if r >= 3), key=lambda g: (-partial[g], g))[:8]

    pos = {n["id"]: (n["x"], -n["y"]) for n in xg_nodes}      # XGMML y grows downward
    cx, cy = CX, -CY
    n_layers = len(TIER_TITLE)                      # 3 layers by design

    fig = plt.figure(figsize=POSTER_FIGSIZE, facecolor=POSTER_BG)

    # ------------------------------------------------------------- header --
    axh = fig.add_axes([0.0, 0.0, 1.0, 1.0], facecolor="none")
    axh.set_xlim(0, 1); axh.set_ylim(0, 1); axh.axis("off")
    axh.add_patch(Rectangle((0.012, 0.930), 0.976, 0.062, fc="#111827", ec="none"))
    right_hdr = axh.transData.transform((0.980, 0))[0]
    t = axh.text(0.024, 0.968, f"Figure 2 | Three-{TIER_WORD} Master-Hub Command Hierarchy & Validated "
                 f"Functional Classes Inside the Full Pan-Antibiotic GRN",
                 color="white", fontsize=25, fontweight="bold", va="center", ha="left")
    _fit_text(fig, t, right_hdr)
    t = axh.text(0.024, 0.944,
                 f"Complete consensus GRN manifold (N = {n_nodes:,} genes; E = {n_edges:,} unique "
                 f"directed edges) - ghosting unvalidated background edges to spotlight the "
                 f"{len(hubs)} Master Hubs, {n_spine} inter-hub spine links and {n_eff} "
                 f"motif-validated effectors",
                 color="#CBD5E1", fontsize=13.5, va="center", ha="left")
    _fit_text(fig, t, right_hdr)

    # --------------------------------------------------------- main circle --
    ax = fig.add_axes([0.0, 0.0, 0.705, 0.915], facecolor=POSTER_BG)

    # class tints (behind everything)
    for k in CLASS_ORDER:
        a0, a1, _n = sector_spans[k]
        ax.add_patch(Wedge((cx, cy), R_ARC_IN, -(a1 - 0.3), -(a0 + 0.3),
                           fc=CLASS_STYLE[k]["light"], ec="none", alpha=0.55, zorder=0))

    # edges, one collection per edge class, per-edge RGBA
    by_tier = OrderedDict()
    for e in xg_edges:
        by_tier.setdefault(e["attrs"]["tier_class"], []).append(e)
    for tier, grp in sorted(by_tier.items(), key=lambda kv: EDGE_Z.get(kv[0], 3)):
        if tier == "Ghosted_FullGRN_Background_Edge":      # faint but visible
            cols = [to_rgba(GHOST_EDGE_COLOR, GHOST_EDGE_ALPHA)] * len(grp)
            lws = [GHOST_EDGE_LW] * len(grp)
        else:
            cols = [to_rgba(e["stroke"], e["opacity"]) for e in grp]
            lws = [e["width"] * 0.6 for e in grp]
        ax.add_collection(LineCollection(
            [(pos[e["source"]], pos[e["target"]]) for e in grp],
            colors=cols, linewidths=lws, zorder=EDGE_Z.get(tier, 3)))

    # nodes, one scatter per (role, shape), per-node RGBA
    groups = OrderedDict()
    for n in xg_nodes:
        groups.setdefault((n["attrs"]["role"], n["shape"]), []).append(n)
    for (role, shape), grp in sorted(groups.items(), key=lambda kv: ROLE_Z.get(kv[0][0], 16)):
        ghost = role == "Ghosted_Unvalidated_GRN_Leaf"
        ax.scatter([pos[n["id"]][0] for n in grp], [pos[n["id"]][1] for n in grp],
                   s=[(n["w"] ** 1.45) / 4.0 * (1.5 if ghost else 1.15) for n in grp],
                   c=[to_rgba(n["fill"], GHOST_NODE_ALPHA if ghost else n["opacity"]) for n in grp],
                   edgecolors=[to_rgba(n["outline"], GHOST_NODE_ALPHA if ghost else n["opacity"])
                               for n in grp],
                   linewidths=[n["border_w"] * 0.5 for n in grp],
                   marker=SHAPE_MARKER.get(shape, "o"), zorder=ROLE_Z.get(role, 16))

    # hub labels:  name [active in k/n conditions]
    pill = dict(boxstyle="round,pad=0.18", fc="white", ec="none", alpha=0.88)
    for n in xg_nodes:
        if n["attrs"]["role"].startswith("Master_Hub") and n["spotlight_label"]:
            ax.annotate(n["spotlight_label"], pos[n["id"]], textcoords="offset points",
                        xytext=(0, -13), ha="center", va="top", fontsize=12.5,
                        fontweight="bold", style="italic", color="#0F172A", zorder=30, bbox=pill)

    # ---- outside labels: every effector found in >= LABEL_MIN_REC conditions
    side_groups = {"L": [], "R": []}
    for k in CLASS_ORDER:
        genes = [g for g in by_class[k] if g in pos and len(b_cond[g]) >= LABEL_MIN_REC]
        if not genes:
            continue
        a0, a1, _n = sector_spans[k]
        side = "R" if math.cos(math.radians((a0 + a1) / 2.0)) >= 0 else "L"
        genes.sort(key=lambda g: -pos[g][1])                       # top -> bottom
        side_groups[side].append((sum(pos[g][1] for g in genes) / len(genes), k, genes))

    label_rows, y_extent = [], 880.0
    for side, groups_ in side_groups.items():
        groups_.sort(key=lambda t: -t[0])
        items = [(k, g) for _, k, genes in groups_ for g in genes]
        if not items:
            continue
        ys = _spread_labels([pos[g][1] for _, g in items], LABEL_GAP)
        y_extent = max(y_extent, max(abs(y - cy) for y in ys) + 60.0)
        xl = cx + (R_LABEL if side == "R" else -R_LABEL)
        for (k, g), y in zip(items, ys):
            hexc = CLASS_STYLE[k]["hex"]
            gx, gy = pos[g]
            ax.plot([gx, xl - (12 if side == "R" else -12)], [gy, y], color=hexc,
                    lw=0.9, alpha=0.65, zorder=8, solid_capstyle="round")
            txt = f"{fmt_gene(g)} ({len(b_cond[g])}/{n_conditions})"
            ax.text(xl, y, txt, color=hexc, fontsize=11, style="italic",
                    ha="left" if side == "R" else "right", va="center", zorder=25)
            label_rows.append([g, fmt_gene(g), k, len(b_cond[g]),
                               "; ".join(COND_DISPLAY.get(c, c) for c in sorted(b_cond[g])),
                               "right" if side == "R" else "left"])

    ax.set_xlim(cx - 1230, cx + 1230)
    ax.set_ylim(cy - y_extent, cy + y_extent)
    ax.set_aspect("equal")
    ax.axis("off")

    # class arcs (drawn above the leader lines) + class name written on each arc
    for k in CLASS_ORDER:
        a0, a1, _n = sector_spans[k]
        ax.add_patch(Wedge((cx, cy), R_ARC_OUT, -(a1 - 0.3), -(a0 + 0.3),
                           width=R_ARC_OUT - R_ARC_IN, fc=CLASS_STYLE[k]["hex"], ec="none", zorder=9))
    if SHOW_CLASS_NAMES:
        r_mid = (R_ARC_IN + R_ARC_OUT) / 2.0
        for k in CLASS_ORDER:
            a0, a1, _n = sector_spans[k]
            span = math.radians(a1 - a0 - 0.6)
            _arc_text(fig, ax, CLASS_STYLE[k]["short"], cx, cy, r_mid,
                      -(a0 + a1) / 2.0, r_mid * span, "white")

    # ------------------------------------------------------- right panel ---
    axp = fig.add_axes([0.708, 0.012, 0.284, 0.906], facecolor="none")
    axp.set_xlim(0, 100); axp.set_ylim(0, 100); axp.axis("off")
    axp.add_patch(FancyBboxPatch((0.4, 0.4), 99.2, 99.2,
                  boxstyle="round,pad=0,rounding_size=2.2", mutation_aspect=0.455,
                  fc="white", ec="#E2E8F0", lw=1.6))
    right_px = axp.transData.transform((96.0, 0))[0]

    def put(x, y, s, size, **kw):
        t = axp.text(x, y, s, fontsize=size, va="center", **kw)
        return _fit_text(fig, t, right_px)

    def head(y, text):
        """Section caption: large, bold, near-black, with a dark rule under it."""
        put(3.9, y, text, 18.5, fontweight="bold", color="#020617")
        axp.plot([3.9, 96.1], [y - 2.1, y - 2.1], color="#334155", lw=2.2, solid_capstyle="butt")

    K = 1.0 if SHOW_SUMMARY_BOX else 1.12      # a little more air when the box is gone

    # A. layers
    head(97.0, f"A. {n_layers}-{TIER_WORD} Master-Hub Command Core ({len(hubs)}/{n_regs} regulators)")
    tier_style = {1: ("H", "#0F172A", "#38BDF8"), 2: ("H", "#1D4ED8", "#93C5FD"),
                  3: ("8", "#0F766E", "#5EEAD4")}
    y = 93.0
    for t_ in (1, 2, 3):
        members = sorted((h for h in hubs if hub_tier_map[h][1] == t_),
                         key=lambda h: HUB_TIERS[h][2].lower())
        mk, fc, ec = tier_style[t_]
        axp.scatter([6.3], [y], marker=mk, s=230, c=fc, edgecolors=ec, linewidths=2.2, zorder=5)
        put(11.5, y, TIER_TITLE[t_], 13, fontweight="bold", color="#0F172A")
        sub = textwrap.wrap(", ".join(f"{HUB_TIERS[h][2]} [{hub_badges[h]}]" for h in members)
                            or "(none)", 62)
        for j, ln in enumerate(sub):
            put(11.5, y - 2.4 - 1.9 * j, ln, 10.3, color="#475569")
        y -= (5.1 + 1.9 * (len(sub) - 1)) * K
    axp.scatter([6.3], [y], marker="D", s=200, c="#6366F1", edgecolors="#4F46E5", linewidths=1.2, zorder=5)
    put(11.5, y, f"Combinatorial Partner TFs ({n_partner})", 13, fontweight="bold", color="#0F172A")
    put(11.5, y - 2.4, "non-hub regulators inside a significant motif", 10.3, color="#475569")
    y -= 5.1 * K
    axp.scatter([6.3], [y], marker="o", s=70, c="#CBD5E1", edgecolors="#94A3B8", linewidths=1, zorder=5)
    put(11.5, y, "Ghosted Full-GRN Background", 13, fontweight="bold", color="#0F172A")
    put(11.5, y - 2.4, f"{len(ghosts):,} unvalidated leaves, {n_bg_tf:,} background TFs, "
        f"{n_ghost_edges:,} non-motif edges", 10.3, color="#475569")

    # B. spine
    yb = y - 6.2 * K
    head(yb, f"B. Inter-Hub {n_layers}-{TIER_WORD} Directed Spine ({n_spine} links)")
    rows = [("#DC2626", 4.2, "dmlR -> fadR conserved primary artery", spine_bb),
            ("#D97706", 2.6, "recurrent multi-drug hub -> hub spines", spine_rec),
            ("#0284C7", 1.4, "condition-specific hub -> hub relays", spine_one)]
    yy = yb - 4.8 * K
    for col, lw, text, vals in rows:
        axp.plot([4.2, 9.8], [yy, yy], color=col, lw=lw, solid_capstyle="butt")
        put(11.5, yy, f"{text}  ({len(vals)}; {rng(vals)} conditions)", 11.3, color="#1E293B")
        yy -= 3.7 * K

    # C. classes
    yc = yy - 1.6 * K
    head(yc, "C. Validated Effectors by Functional Class")
    put(3.9, yc - 3.7, f"({n_eff} of {n_terminal:,} terminal genes)", 11, color="#475569")
    maxc = max([len(by_class[k]) for k in CLASS_ORDER] + [1])
    yy = yc - 7.6
    lower, max_gap = (19.6, 2.95) if SHOW_SUMMARY_BOX else (6.5, 4.4)
    row_gap = min(max_gap, max(1.9, (yy - lower) / max(1, len(CLASS_ORDER) - 1)))
    for k in CLASS_ORDER:
        n = len(by_class[k])
        axp.scatter([6.3], [yy], marker="o", s=80, c=CLASS_STYLE[k]["hex"], zorder=5)
        put(10.4, yy, CLASS_STYLE[k]["short"], 11.8, color="#1E293B")
        axp.add_patch(Rectangle((62, yy - 0.8), 26.0, 1.6, fc="#E2E8F0", ec="none"))
        axp.add_patch(Rectangle((62, yy - 0.8), 26.0 * n / maxc, 1.6, fc=CLASS_STYLE[k]["hex"], ec="none"))
        axp.text(90.8, yy, f"{n}", fontsize=11.5, color="#1E293B", va="center")
        yy -= row_gap

    # optional recurrence summary box
    if SHOW_SUMMARY_BOX:
        axp.add_patch(FancyBboxPatch((3.0, 1.4), 94.0, 15.6,
                      boxstyle="round,pad=0,rounding_size=1.6", mutation_aspect=0.455,
                      fc="#FFFBEB", ec="#FCD34D", lw=1.5))
        nw = NUM_WORD.get(n_conditions, str(n_conditions))
        title = (f"No universal effector set across the {nw} antibiotics" if k_all == 0
                 else f"{k_all} effectors shared by all {nw} antibiotics")
        put(5.2, 14.6, title, 12.2, fontweight="bold", color="#92400E")
        lines = [f"{k_all} / {n_eff} ({pct(k_all):.1f}%) validated effectors shared by all {n_conditions} conditions",
                 f"{k_one} / {n_eff} ({pct(k_one):.1f}%) unique to a single antibiotic"]
        if partial:
            lines.append(f"{len(partial)} / {n_eff} ({pct(len(partial)):.1f}%) partially shared across "
                         f"{min(partial.values())}-{max(partial.values())} conditions"
                         + (", led by" if lead else ""))
            if lead:
                lines += textwrap.wrap(", ".join(f"{fmt_gene(g)} ({partial[g]}/{n_conditions})"
                                                 for g in lead), 68)
        yy = 11.6
        for ln in lines:
            put(5.6, yy, ln, 10.2, color="#7C2D12")
            yy -= 2.1

    # ------------------------------------------------------------- write ---
    out_png = out_dir / "Fig2_poster.png"
    out_svg = out_dir / "Fig2_poster.svg"
    fig.savefig(out_png, dpi=POSTER_DPI, facecolor=POSTER_BG)
    fig.savefig(out_svg, facecolor=POSTER_BG)
    plt.close(fig)

    with open(out_dir / "Fig2_Outer_Labels.csv", "w", newline="", encoding="utf-8") as fh:
        w_ = csv.writer(fh)
        w_.writerow(["gene_id", "label", "functional_class", "recurrence",
                     "conditions", "label_side"])
        w_.writerows(sorted(label_rows, key=lambda r: (CLASS_ORDER.index(r[2]), -r[3], r[1])))
    print(f"poster: {len(label_rows)} outside labels (recurrence >= {LABEL_MIN_REC})")
    print(f"wrote {out_png}")
    print(f"wrote {out_svg}")


if __name__ == "__main__":
    main()
