# *Acinetobacter baumannii* ATCC 19606 — Antibiotic Resistance GRN Pipeline

End-to-end pipeline for constructing and analysing condition-specific gene regulatory networks (GRNs) in *A. baumannii* ATCC 19606 under seven antibiotic conditions: **amikacin, ciprofloxacin, doripenem, meropenem, polymyxin, tigecycline**, and a **colistin-resistant** strain.

> **🗄️ Database:** The consensus GRN networks produced by this pipeline are publicly available in **AbRegNet** — the *Acinetobacter baumannii* GRN Database.
> Integrating transcriptomics data from 5 inference methods to decode antibiotic resistance, regulatory mechanisms, and adaptive responses.
> 🔗 https://bic.sastra.edu/grnnet/

---

## Pipeline overview

```
Raw FASTQs
    │
    ├─ 1. Quality control         fastqc.nf
    ├─ 2. Adapter trimming        run_fasp.sh  (fastp, Q20)
    ├─ 3. Alignment & counting    index_to_counts.nf
    │                              Bowtie2 → SAMtools → featureCounts
    │
    ├─ 4. Count merging           1_merge_count.ipynb
    │      Merge *_counts.txt → merged_counts.csv
    │      Filter strain-zero genes
    │      Split per GEO ID → GEO_grouped/<GEO>_counts.csv
    │
    ├─ 5. Normalisation           upper_quartile_normalisation.R
    │      Upper-quartile norm (edgeR) → log2 CPM per GEO dataset
    │      Merge all GEO datasets → all_GEO_merged_logCPM.csv
    │
    ├─ 6. Batch correction        imodulon_combat.R
    │      ComBat (SVA) across 7 GEO batches → global_ComBat_logCPM.csv
    │      Final merged matrix   → all96_merged_logCPM_ComBat.csv
    │      Per-antibiotic slices → <Drug>_for_GRN.csv  (input to comhub)
    │
    ├─ 7.  GRN construction (ComHub  ·  iModulonMiner  (external))
    │      → per-antibiotic GRN edge CSVs (source, target)
    │
    ├─ 8. High-confidence filter  high_threshold_across_networks.py
    ├─ 9. Cross-method consensus  consensus_across_antibiotics.py
    │                              + RegulonDB validation
    │                              → networks stored in AbRegNet database
    │                                https://bic.sastra.edu/grnnet/
    │
    ├─ 10. GRN analysis
    │       step1_1_summary_table.py           topology, FFLs, bi-fans, hub chains
    │       step1_2_centrality_analysis.ipynb  degree / betweenness / closeness
    │       step2_tiers_seesaw.py              3-tier hub stratification + seesaw test
    │       step3_hub_motifs.py                hub-driven FFLs and bi-fans
    │       step4_motif_enrichment.py          Maslov–Sneppen null-model (n=1,000)
    │       step5_coverage_leaves.py           hub motif leaf-gene coverage
    │       step6_regulatory_network_tiers.py  pan-network tier figure + Cytoscape export
    │       step7_pan_grn_analysis.py          Fig. 2 pan-GRN XGMML
    │       regulator_analysis.ipynb           TF family classification
    │
    ├─ 11. Cross-strain motifs    cross_strain_compare.py
    │
    └─ 12. Ortholog & synteny
            s2_hub_symbol_orthologs.py   RBH ortholog propagation (BLAST+)
            s3_microsynteny.py           chromosomal S_syn scoring (±2 CDS)
            s4_validate_synteny.py       background / control / window validation
```

---

## Steps

### 1 · Quality control — `fastqc.nf`

Runs FastQC on every FASTQ in a directory via Nextflow.

```bash
nextflow run fastqc.nf \
  --fastq     'raw_reads/*.fastq' \
  --qc_report results/fastqc_report
```

Outputs HTML + ZIP per sample to `--qc_report`.

---

### 2 · Adapter trimming — `run_fasp.sh`

Loops over every `*.fastq` in the current directory and runs fastp (Q20, 30 % unqualified-base cap, 4 threads).

```bash
bash run_fasp.sh
# outputs: fastp_output/<sample>_trimmed.fastq
#          fastp_output/<sample>_fastp.html / .json
```

---

### 3 · Alignment and counting — `index_to_counts.nf`

Bowtie2 index → alignment → SAM→BAM → sort → flagstat → featureCounts.

| Parameter | Description | Default |
|---|---|---|
| `--ref` | Reference genome FASTA | *(required)* |
| `--reads` | Glob for input FASTQs | *(required)* |
| `--annotation` | GFF / GTF annotation | *(required)* |
| `--indexdir` | Bowtie2 index output | `./index` |
| `--alignedir` | BAM output | `./alignment` |
| `--countdir` | Count output | `./counts` |
| `--single_end` | `true` for single-end reads | `false` |

```bash
# Paired-end
nextflow run index_to_counts.nf \
  --ref        reference/ATCC19606.fasta \
  --reads      'fastp_output/*_{1,2}.fastq' \
  --annotation reference/ATCC19606.gff \
  --alignedir  results/alignment \
  --countdir   results/counts

# Single-end — append: --single_end true
```

featureCounts flags: `-t gene -g gene_id`; paired-end adds `-p --countReadPairs`. BAMs are name-sorted for single-end and coordinate-sorted for paired-end.

---

### 4 · Count merging — `1_merge_count.ipynb`

Three steps in one notebook:

1. Reads all `*_counts.txt` from featureCounts, keeps `Geneid` + count column per sample, outer-joins into `merged_counts.csv`.
2. Filters genes that are zero across **all** samples of any one strain (strain-zero filter), writing `merged_counts_strain_filtered.csv`. Requires `metadata.final.csv` with columns `Run`, `Strain`.
3. Groups samples by GEO dataset ID and saves one file per GEO dataset to `GEO_grouped/<GEO>_counts.csv`. Requires `metadata.final.csv` with column `GEO ID`.

```bash
jupyter nbconvert --to notebook --execute 1_merge_count.ipynb
```

---

### 5 · Normalisation — `upper_quartile_normalisation.R`

Applies upper-quartile normalisation (edgeR) followed by log₂ CPM conversion to each GEO dataset separately, then merges all datasets into one matrix.

```r
# Runs automatically over all *_counts.csv in the working directory
Rscript upper_quartile_normalisation.R
```

**What it does, per GEO dataset:**
- Reads `<GEO>_counts.csv` from `GEO_grouped/`
- Creates a DGEList and calls `calcNormFactors(method = "upperquartile")`
- Converts to log₂ CPM with `prior.count = 1`
- Saves `<GEO>_upperquartile_logCPM.csv`

**Then merges all datasets:**
- Checks gene row identity across files
- Column-binds into `all_GEO_merged_logCPM.csv`

Requires `metadata.final.csv` (columns: `Run`, `GEO ID`).

---

### 6 · Batch correction — `imodulon_combat.R`

Removes GEO-dataset batch effects using ComBat (SVA package) while preserving the biological antibiotic signal.

```r
Rscript imodulon_combat.R
```

**What it does:**

| Step | Detail |
|---|---|
| Input | `all_GEO_merged_logCPM.csv` + `metadata_global_dataset_wise.csv` |
| Batch variable | `GEO ID` (7 batches, ~90 samples) |
| Biological covariate | `~ Antibiotic` (preserved in `mod`) |
| Global ComBat output | `global_ComBat_logCPM.csv` |
| Tigecycline | Single batch — no correction; saved as `tigecycline_logCPM_no_ComBat.csv` |
| Final merged matrix | `all96_merged_logCPM_ComBat.csv` (96 samples) |
| Per-antibiotic slices | `Amikacin_for_GRN.csv`, `Ciprofloxacin_for_GRN.csv`, `Meropenem_for_GRN.csv`, `Colistin-Polymyxin_for_GRN.csv` |

The per-antibiotic slices are the direct input to ComHub / iModulonMiner (step 7).

---

### 7 · ICA decomposition — ComHub & iModulonMiner *(external)*

Run ComHub and iModulonMiner on the per-antibiotic expression matrices from step 6 to produce GRN edge CSVs. See their documentation:
- ComHub: https://gitlab.com/Gustafsson-lab/comhub
- iModulonMiner: https://github.com/SBRG/iModulonMiner

The GRN analysis scripts expect seven edge files with columns `source, target`, placed beside the scripts, in `data/`, or at the path set by `GRN_DATA_DIR`:

| Condition | File |
|---|---|
| Amikacin | `ami.grn.csv` |
| Ciprofloxacin | `cip.grn.csv` |
| Doripenem | `dori.grn.csv` |
| Meropenem | `mero.grn.csv` |
| Polymyxin | `poly.grn.csv` |
| Tigecycline | `tige.grn.csv` |
| Colistin-resistant | `colis.resis.grn.csv` |

Gene symbols are lowercased automatically (`FadR` → `fadr`).

---

### 8–9 · Consensus analysis

```bash
# Apply literature-recommended thresholds per method
# Input: 2.networks/<antibiotic>/<method>.csv  (GeneA_Symbol, GeneB_Symbol, Score)
python high_threshold_across_networks.py

# Merge methods, test 2+/3+/4+/all-5 support levels, validate vs RegulonDB
# Input: high_confidence_output/  +  regulondb.csv
python consensus_across_antibiotics.py
```

| Method | Threshold |
|---|---|
| GENIE3 | Top 1 % |
| ARACNE | Top 1 % |
| PCC | Score ≥ 0.8 |
| CLR (R) | Score ≥ 2.5 |
| ElasticNet bootstrap | Score ≥ 0.6 |

---

### 10 · GRN analysis

All scripts import `grn_common.py`. Run from the directory containing the scripts and GRN edge files.

```bash
python step1_1_summary_table.py                # topology, FFLs, bi-fans, hub chains
jupyter nbconvert --execute step1_2_centrality_analysis.ipynb
python step2_tiers_seesaw.py                   # 3-tier stratification 
python step3_hub_motifs.py                     # hub-driven FFLs and bi-fans
python step4_motif_enrichment.py               # null-model enrichment (swap, default)
python step5_coverage_leaves.py                # leaf-gene coverage
python step6_regulatory_network_tiers.py       # pan-network tier figure
python step7_pan_grn_analysis.py               # pan-GRN XGMML
jupyter nbconvert --execute regulator_analysis.ipynb
```

**11 Master Hub regulators** (second-stage significance filter in `step4`):

> *arsR2 · betI · dmlR · fadR · lrp · paaX · pgrR · rcsC · rhaS · rpoH · yybR*

**Per-condition hub panels** are defined in `grn_common.py → HUBS_BY_CONDITION`.

---

### 11 · Cross-strain motif comparison — `cross_strain_compare.py`

```bash
cd MOTIF          # folder with antibiotic subfolders containing strain results
python3 cross_strain_compare.py
# outputs: cross_strain_results/Table1–5 CSVs
```

---

### 12 · Ortholog and synteny analysis

Requires BLAST+ on `PATH`, strain FASTAs in `data_strains/<strain>.faa`, `hub_ids.csv`, and `hub_proteins.faa`.

```bash
python s2_hub_symbol_orthologs.py   # RBH ortholog propagation (E≤1e-10, cov≥85%, id≥85%)
python s3_microsynteny.py           # S_syn chromosomal neighbourhood scoring (window ±2 CDS)
python s4_validate_synteny.py       # 4-check validation: background, controls, windows, self-consistency
```

---

## Requirements

| Tool | Version | Used in |
|---|---|---|
| Nextflow | ≥ 23 | steps 1, 3 |
| FastQC | ≥ 0.11 | step 1 |
| fastp | ≥ 0.23 | step 2 |
| Bowtie2 | ≥ 2.5 | step 3 |
| SAMtools | ≥ 1.17 | step 3 |
| featureCounts (Subread) | ≥ 2.0 | step 3 |
| Python | ≥ 3.9 | steps 4, 8–12 |
| pandas, numpy | — | steps 4, 8–12 |
| networkx | ≥ 3.0 | steps 8–11 |
| scipy, statsmodels, matplotlib, seaborn | — | steps 8–12 |
| R | ≥ 4.3 | steps 5, 6 |
| edgeR | ≥ 3.40 | step 5 |
| sva (ComBat) | ≥ 3.46 | step 6 |
| BLAST+ | ≥ 2.14 | step 12 |
| ComHub, pymodulon | — | step 7 |

```bash
# Python
pip install networkx pandas numpy scipy statsmodels matplotlib seaborn comhub pymodulon

# R
install.packages("BiocManager")
BiocManager::install(c("edgeR", "sva"))
```

---

## AbRegNet — GRN Database

The consensus GRN networks produced by this pipeline are deposited in **AbRegNet**, a publicly accessible web database hosted at SASTRA University:

🔗 **https://bic.sastra.edu/grnnet/**

AbRegNet integrates transcriptomics data from five GRN inference methods (CLR, PCC, ARACNE, GENIE3, ElasticNet) to decode antibiotic resistance, regulatory mechanisms, and adaptive responses in *A. baumannii*. Users can browse, search, and download:

- Condition-specific GRN edge lists (per antibiotic)
- Consensus networks at multiple method-support levels


---

## iModulon analysis

Coming soon — functional GO/COG/KEGG annotation, differential iModulon activity between conditions.

---


