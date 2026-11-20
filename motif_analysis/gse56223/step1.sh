#!/usr/bin/env bash
# STEP 1 — Map symbols to locus_tags via reference_genes.csv
# FIX: when a symbol maps to multiple locus_tags (paralogues),
#      ALL copies are included as separate edges
echo "============================================"
echo " STEP 1 — Mapping symbols to locus_tags"
echo "============================================"

TARGETS="ami.tsv"
REF_GENES="reference_genes.csv"

for f in "$TARGETS" "$REF_GENES"; do
    [[ -f "$f" ]] || { echo "  ERROR: $f not found"; exit 1; }
done

sed -i 's/\r//' "$TARGETS"
sed -i 's/\r//' "$REF_GENES"
echo "  ✓ Cleaned line endings"

mkdir -p beds fastas meme_out fimo_out results

python3 - <<'PYMAP'
from collections import defaultdict

# ── Load reference_genes.csv ─────────────────────────────────────────
# sym_to_loci: symbol → LIST of all locus_tags (handles paralogues)
sym_to_loci = defaultdict(list)
locus_to_sym = {}

with open("reference_genes.csv") as f:
    next(f)
    for line in f:
        parts = line.strip().split(",")
        if len(parts) < 2: continue
        locus = parts[0].strip()
        sym   = parts[1].strip().lower()
        locus_to_sym[locus] = sym
        sym_to_loci[sym].append(locus)   # ← store ALL copies, not just first

# Report duplicates found
multi = {s:l for s,l in sym_to_loci.items() if len(l)>1}
print(f"  ✓ reference_genes.csv: {len(sym_to_loci)} unique symbols loaded")
print(f"  ℹ  Symbols with multiple locus_tags (paralogues): {len(multi)}")

# ── Convert ami.tsv edges — expand duplicates ─────────────────────────
edges   = []
missing = []
skipped_dup = []

with open("ami.tsv") as f:
    next(f)
    for line in f:
        parts = line.strip().split("\t")
        if len(parts) < 2: continue
        reg_sym = parts[0].strip().lower()
        tgt_sym = parts[1].strip().lower()

        reg_loci = sym_to_loci.get(reg_sym, [])
        tgt_loci = sym_to_loci.get(tgt_sym, [])

        if not reg_loci:
            missing.append(f"  SKIP regulator '{reg_sym}' — not in reference_genes.csv")
            continue
        if not tgt_loci:
            missing.append(f"  SKIP target '{tgt_sym}' — not in reference_genes.csv")
            continue

        # Use FIRST locus_tag for regulator (regulator is one specific TF)
        reg_locus = reg_loci[0]

        # Expand ALL copies for target (include every paralogue)
        for tgt_locus in tgt_loci:
            edges.append((reg_sym, tgt_sym, reg_locus, tgt_locus))

        if len(tgt_loci) > 1:
            skipped_dup.append(
                f"  EXPANDED '{tgt_sym}' → {len(tgt_loci)} copies: {tgt_loci}"
            )

# ── Write ami_locus.tsv ───────────────────────────────────────────────
with open("ami_locus.tsv", "w") as f:
    f.write("Regulator_sym\tTarget_sym\tRegulator_locus\tTarget_locus\n")
    for r_sym, t_sym, r_locus, t_locus in edges:
        f.write(f"{r_sym}\t{t_sym}\t{r_locus}\t{t_locus}\n")

# ── Write locus→symbol lookup ─────────────────────────────────────────
with open("locus_to_symbol.tsv", "w") as f:
    f.write("Locus_tag\tSymbol\n")
    for locus, sym in sorted(locus_to_sym.items()):
        f.write(f"{locus}\t{sym}\n")

print(f"  ✓ ami_locus.tsv: {len(edges)} edges (original 42 + paralogue expansions)")
if missing:
    for m in missing: print(m)
else:
    print(f"  ✓ All symbols found in reference_genes.csv")

if skipped_dup:
    print(f"\n  Paralogue expansions (all copies included):")
    for s in skipped_dup: print(s)
PYMAP

echo ""
echo "  Regulators: $(awk 'NR>1{print $1}' ami_locus.tsv | sort -u | tr '\n' ' ')"
echo "  Total edges (including paralogue copies): $(awk 'NR>1' ami_locus.tsv | wc -l)"
echo ""
echo "  STEP 1 DONE ✓"
