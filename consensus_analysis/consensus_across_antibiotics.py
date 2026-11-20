
import pandas as pd
import numpy as np
import os


# INPUT / OUTPUT CONFIGURATION


# Root folder containing one subfolder per antibiotic, each with the five
# threshold-filtered method CSVs inside it:
#   INPUT_ROOT/<antibiotic>/clr_R_network.csv
#   INPUT_ROOT/<antibiotic>/pcc_network.csv
#   INPUT_ROOT/<antibiotic>/aracne_network.csv
#   INPUT_ROOT/<antibiotic>/genie3_network.csv
#   INPUT_ROOT/<antibiotic>/elasticnet_bootstrap_network.csv
INPUT_ROOT = "high_confidence_output"

NETWORK_FILES = {
    "CLR": "clr_R_network.csv",
    "PCC": "pcc_network.csv",
    "ARACNE": "aracne_network.csv",
    "GENIE3": "genie3_network.csv",
    "ElasticNet": "elasticnet_bootstrap_network.csv"
}

REGULONDB_FILE = "regulondb.csv"

OUTPUT_DIR = "consensus_analysis"
os.makedirs(OUTPUT_DIR, exist_ok=True)

TOTAL_METHODS = len(NETWORK_FILES)
SUPPORT_LEVELS = [2, 3, 4, TOTAL_METHODS]


# LOAD REGULONDB (once - shared across all antibiotics)


print("\nLoading RegulonDB...")

regdb = pd.read_csv(REGULONDB_FILE)

regdb["Regulator"] = (
    regdb["Regulator"]
    .astype(str)
    .str.lower()
    .str.strip()
)

regdb["Regulated_gene"] = (
    regdb["Regulated_gene"]
    .astype(str)
    .str.lower()
    .str.strip()
)

# create undirected validated pairs
validated_pairs = set()

for _, row in regdb.iterrows():

    g1 = min(row["Regulator"], row["Regulated_gene"])
    g2 = max(row["Regulator"], row["Regulated_gene"])

    validated_pairs.add((g1, g2))

print(f"Loaded RegulonDB interactions: {len(validated_pairs)}")



# PER-ANTIBIOTIC CONSENSUS ANALYSIS


def run_consensus_for_antibiotic(antibiotic, antibiotic_input_dir, antibiotic_output_dir):
    """
    Runs the full combine -> consensus -> support-level -> RegulonDB
    validation workflow for a single antibiotic folder. Returns a list of
    summary dict rows (one per support level) for this antibiotic.
    """
    os.makedirs(antibiotic_output_dir, exist_ok=True)

    print("\n" + "#" * 60)
    print(f"ANTIBIOTIC: {antibiotic}")
    print("#" * 60)


    # LOAD NETWORKS FOR THIS ANTIBIOTIC


    print("\nLoading GRN networks...")

    all_edges = []

    for method, file in NETWORK_FILES.items():

        print("\n")
        print(f"METHOD: {method}")
       

        file_path = os.path.join(antibiotic_input_dir, file)

        if not os.path.exists(file_path):
            print(f"File not found: {file_path}")
            continue

        df = pd.read_csv(file_path)

        print(f"Original edges: {len(df)}")

    
        # CLEAN SYMBOLS
    

        df["GeneA_Symbol"] = (
            df["GeneA_Symbol"]
            .astype(str)
            .str.lower()
            .str.strip()
        )

        df["GeneB_Symbol"] = (
            df["GeneB_Symbol"]
            .astype(str)
            .str.lower()
            .str.strip()
        )

        # remove nan and empty symbols
        df = df[
            (df["GeneA_Symbol"] != "nan") &
            (df["GeneB_Symbol"] != "nan") &
            (df["GeneA_Symbol"] != "") &
            (df["GeneB_Symbol"] != "")
        ]

        print(f"After removing missing symbols: {len(df)}")

    
        # REMOVE SELF LOOPS
    

        df = df[
            df["GeneA_Symbol"] != df["GeneB_Symbol"]
        ]

        print(f"After removing self-loops: {len(df)}")

    
        # CREATE UNDIRECTED PAIRS
    

        df["g1"] = df.apply(
            lambda x: min(x["GeneA_Symbol"], x["GeneB_Symbol"]),
            axis=1
        )

        df["g2"] = df.apply(
            lambda x: max(x["GeneA_Symbol"], x["GeneB_Symbol"]),
            axis=1
        )

    
        # KEEP STRONGEST EDGE
    

        df = (
            df.sort_values("Score", ascending=False)
              .drop_duplicates(["g1", "g2"])
        )

        print(f"Unique edges retained: {len(df)}")

    
        # ADD METHOD NAME
    

        df["method"] = method

        # keep needed columns
        df = df[[
            "g1",
            "g2",
            "Score",
            "method"
        ]]

        all_edges.append(df)

    if not all_edges:
        print(f"\nNo network files found for antibiotic '{antibiotic}' - skipping.")
        return []

    
    # COMBINE ALL METHODS (this antibiotic only)
    

    print(f"COMBINING NETWORKS - {antibiotic}")
  

    combined = pd.concat(all_edges, ignore_index=True)

    print(f"Total combined edges: {len(combined)}")

    combined.to_csv(
        os.path.join(antibiotic_output_dir, "all_combined_edges.csv"),
        index=False
    )

    
    # BUILD CONSENSUS TABLE
    

    print("\nBuilding consensus network...")

    consensus = (
        combined
        .groupby(["g1", "g2"])
        .agg({
            "method": lambda x: sorted(set(x)),
            "Score": ["mean", "max", "count"]
        })
    )

    # flatten columns
    consensus.columns = [
        "methods",
        "mean_score",
        "max_score",
        "support"
    ]

    consensus = consensus.reset_index()

    # convert methods to string
    consensus["method_string"] = consensus["methods"].apply(
        lambda x: ",".join(x)
    )

    print(f"Unique consensus edges: {len(consensus)}")

    consensus.to_csv(
        os.path.join(antibiotic_output_dir, "consensus_all_edges.csv"),
        index=False
    )

    
    # TEST CONSENSUS LEVELS
    

  
    print(f"CONSENSUS BENCHMARKING - {antibiotic}")
    

    antibiotic_summary = []

    for level in SUPPORT_LEVELS:

        print(f"CONSENSUS LEVEL")
      
        # SELECT SUPPORT LEVEL
    

        if level == TOTAL_METHODS:

            subset = consensus[
                consensus["support"] == TOTAL_METHODS
            ].copy()

            label = f"all_{TOTAL_METHODS}_methods"

        else:

            subset = consensus[
                consensus["support"] >= level
            ].copy()

            label = f"{level}_plus_methods"

        print(f"Group: {label}")

    
        # VALIDATE AGAINST REGULONDB
    

        subset["validated"] = subset.apply(
            lambda row: (row["g1"], row["g2"]) in validated_pairs,
            axis=1
        )

        validated = subset[
            subset["validated"]
        ].copy()

        novel = subset[
            ~subset["validated"]
        ].copy()

    
        # STATISTICS
    

        total_edges = len(subset)
        validated_edges = len(validated)
        novel_edges = len(novel)

        validation_rate = (
            validated_edges / total_edges * 100
            if total_edges > 0 else 0
        )

    
        # SAVE FILES
    

        subset.to_csv(
            os.path.join(
                antibiotic_output_dir,
                f"consensus_{label}.csv"
            ),
            index=False
        )

        validated.to_csv(
            os.path.join(
                antibiotic_output_dir,
                f"validated_{label}.csv"
            ),
            index=False
        )

        novel.to_csv(
            os.path.join(
                antibiotic_output_dir,
                f"novel_{label}.csv"
            ),
            index=False
        )

    
        # PRINT RESULTS
    

        print(f"Total edges      : {total_edges}")
        print(f"Validated edges  : {validated_edges}")
        print(f"Novel edges      : {novel_edges}")
        print(f"Validation rate  : {validation_rate:.2f}%")

    
        # STORE SUMMARY
    

        antibiotic_summary.append({
            "antibiotic": antibiotic,
            "support_level": label,
            "total_edges": total_edges,
            "validated_edges": validated_edges,
            "novel_edges": novel_edges,
            "validation_rate": round(validation_rate, 2)
        })

    # per-antibiotic summary file (same as before, just scoped to this folder)
    antibiotic_summary_df = pd.DataFrame(antibiotic_summary)
    antibiotic_summary_df.to_csv(
        os.path.join(antibiotic_output_dir, "consensus_summary.csv"),
        index=False
    )

    return antibiotic_summary



# MAIN: LOOP OVER EVERY ANTIBIOTIC FOLDER


if not os.path.exists(INPUT_ROOT):
    raise SystemExit(f"Error: input root '{INPUT_ROOT}' not found!")

all_summaries = []

antibiotics = sorted(
    d for d in os.listdir(INPUT_ROOT)
    if os.path.isdir(os.path.join(INPUT_ROOT, d))
)

if not antibiotics:
    raise SystemExit(f"Error: no antibiotic subfolders found under '{INPUT_ROOT}'!")

for antibiotic in antibiotics:
    antibiotic_input_dir = os.path.join(INPUT_ROOT, antibiotic)
    antibiotic_output_dir = os.path.join(OUTPUT_DIR, antibiotic)

    summary_rows = run_consensus_for_antibiotic(
        antibiotic, antibiotic_input_dir, antibiotic_output_dir
    )
    all_summaries.extend(summary_rows)


# COMBINED SUMMARY ACROSS ALL ANTIBIOTICS


combined_summary_df = pd.DataFrame(all_summaries)

combined_summary_file = os.path.join(OUTPUT_DIR, "combined_consensus_summary.csv")
combined_summary_df.to_csv(combined_summary_file, index=False)


print("FINAL COMBINED SUMMARY (all antibiotics)")


print(combined_summary_df)

print(f"\nCombined summary saved: {combined_summary_file}")
print(f"Per-antibiotic outputs saved under: {OUTPUT_DIR}/<antibiotic>/")

print("\nDONE.")