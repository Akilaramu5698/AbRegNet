"""
s1_fetch_strains.py -- download proteomes and ordered CDS annotation for the
six A. baumannii benchmark strains.


Outputs (data_strains/)
  <strain>.faa        protein FASTA, header rewritten to >locus_tag strain=...
  <strain>_cds.csv    ordered CDS table: locus_tag, gene, product, start, end,
                      strand, cds_index (genome order, per replicon)
  fetch_summary.csv   what was downloaded, with CDS counts


Run:  python s1_fetch_strains.py
"""
from __future__ import annotations

import os
import re
import time
import urllib.parse
import urllib.request

import pandas as pd

EMAIL = ""            # <- put your address here (NCBI etiquette)
TOOL = "grn-ortholog-pipeline"
OUT = "data_strains"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# strain -> (nuccore accessions for all replicons, locus prefix, note)
STRAINS = {
    "ATCC_19606": (["CP045110.1"], "FQU82_", "type strain, ST52"),
    "ATCC_17978": (["CP000521.1"], "A1S_",   "lab reference, ST437"),
    "AB5075_UW":  (["CP008706.1"], "ABUW_",  "IC-I ST1 XDR"),
    "AYE":        (["CU459141.1"], "ABAYE",  "IC-I ST1 XDR"),
    "ACICU":      (["CP000863.1"], "ACICU_", "IC-II ST2 CR-XDR"),
    "MDR_ZJ06":   (["CP001937.2"], "ABZJ_",  "IC-II ST2 CR-XDR"),
}

HEADER_FIELD = re.compile(r"\[(\w+)=([^\]]*)\]")
LOCATION = re.compile(r"(complement\()?(?:join\()?<?(\d+)\.\.>?(\d+)")


def efetch(acc: str) -> str:
    params = {"db": "nuccore", "id": acc, "rettype": "fasta_cds_aa",
              "retmode": "text", "tool": TOOL}
    if EMAIL:
        params["email"] = EMAIL
    url = f"{EUTILS}?{urllib.parse.urlencode(params)}"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if attempt == 3:
                raise
            time.sleep(3 * (attempt + 1))


def parse_cds_fasta(text: str, strain: str, replicon: str):
    """Return (list of (header_id, sequence), DataFrame of CDS annotation)."""
    records, rows, seq, cur = [], [], [], None
    for line in text.splitlines():
        if line.startswith(">"):
            if cur:
                records.append((cur["locus_tag"], "".join(seq)))
            seq = []
            fields = dict(HEADER_FIELD.findall(line))
            loc = fields.get("location", "")
            m = LOCATION.search(loc)
            start, end = (int(m.group(2)), int(m.group(3))) if m else (0, 0)
            strand = "-" if (m and m.group(1)) or loc.startswith("complement") else "+"
            lt = fields.get("locus_tag") or fields.get("protein_id") or f"{replicon}_{len(rows)}"
            cur = {"strain": strain, "replicon": replicon, "locus_tag": lt,
                   "gene": fields.get("gene", ""), "product": fields.get("protein", ""),
                   "protein_id": fields.get("protein_id", ""),
                   "start": start, "end": end, "strand": strand,
                   "pseudo": "pseudo" in fields}
            rows.append(cur)
        elif cur is not None:
            seq.append(line.strip())
    if cur:
        records.append((cur["locus_tag"], "".join(seq)))
    df = pd.DataFrame(rows)
    if len(df):
        df = df.sort_values(["replicon", "start"]).reset_index(drop=True)
        df["cds_index"] = df.groupby("replicon").cumcount()
    return records, df


def main():
    os.makedirs(OUT, exist_ok=True)
    if not EMAIL:
        print("  note: EMAIL is empty - set it in the script; NCBI throttles "
              "anonymous requests")
    summary = []
    for strain, (accs, prefix, note) in STRAINS.items():
        all_recs, frames = [], []
        for acc in accs:
            print(f"[fetch] {strain} {acc} ...", flush=True)
            txt = efetch(acc)
            if not txt.startswith(">"):
                raise SystemExit(f"{strain}/{acc}: unexpected response "
                                 f"(first 200 chars): {txt[:200]}")
            recs, df = parse_cds_fasta(txt, strain, acc)
            all_recs += recs
            frames.append(df)
            time.sleep(0.4)                       # be polite to the API
        cds = pd.concat(frames, ignore_index=True)
        faa = os.path.join(OUT, f"{strain}.faa")
        with open(faa, "w", encoding="utf-8") as fh:
            for lt, s in all_recs:
                if s:
                    fh.write(f">{lt} strain={strain}\n")
                    for i in range(0, len(s), 60):
                        fh.write(s[i:i + 60] + "\n")
        cds.to_csv(os.path.join(OUT, f"{strain}_cds.csv"), index=False)
        observed = cds["locus_tag"].str.startswith(prefix).mean() if len(cds) else 0
        summary.append({"strain": strain, "accessions": ";".join(accs),
                        "expected_locus_prefix": prefix,
                        "frac_loci_matching_prefix": round(float(observed), 3),
                        "n_cds": len(cds), "n_proteins_written": len(all_recs),
                        "note": note})
        print(f"  {strain}: {len(cds)} CDS, prefix match "
              f"{observed:.1%}")
    s = pd.DataFrame(summary)
    s.to_csv(os.path.join(OUT, "fetch_summary.csv"), index=False)
    print(s.to_string(index=False))
    bad = s[s["frac_loci_matching_prefix"] < 0.5]
    if len(bad):
        print("\nWARNING: locus prefix mismatch for: "
              + ", ".join(bad['strain']) + " - check the accession is the "
              "strain you intended before trusting downstream results.")


if __name__ == "__main__":
    main()
