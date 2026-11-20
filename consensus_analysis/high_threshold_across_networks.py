# Applying literature-recommended default thresholds to extract the HIGH-CONFIDENCE final network per method

#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
warnings.filterwarnings('ignore')


# CONFIGURATION


network_root = "2.networks"

# Literature-recommended default thresholds (NO optimization) - used to
# extract the high-confidence final network per method.
LITERATURE_THRESHOLDS = {
    "genie3_network": {"type": "percent", "value": 1},
    "aracne_network": {"type": "percent", "value": 1},
    "pcc_network": {"type": "absolute", "value": 0.8},
    "clr_R_network": {"type": "absolute", "value": 2.5},
    "elasticnet_bootstrap_network": {"type": "absolute", "value": 0.6},
}
LITERATURE_THRESHOLDS_LOWER = {k.lower(): v for k, v in LITERATURE_THRESHOLDS.items()}

# Output root - mirrors the input `networks/<antibiotic>/<method>.csv`
# structure: high_confidence_output/<antibiotic>/<method>.csv, one file per
# method per antibiotic (instead of one combined CSV across everything).
output_root = "high_confidence_output"



# APPLY THRESHOLD (literature default)


def apply_literature_threshold(mapped_sorted, method, total_edges, config):
    if config["type"] == "percent":
        cutoff_n = int(total_edges * (config["value"] / 100.0))
        cutoff_n = max(1, cutoff_n)
        subset = mapped_sorted.head(cutoff_n).copy()
        label = f"Top {config['value']}%"
        threshold_value = config["value"]
    else:
        subset = mapped_sorted[mapped_sorted["Score"] >= config["value"]].copy()
        label = f">= {config['value']}"
        threshold_value = config["value"]

    return subset, label, threshold_value


def _load_and_clean_network_file(antibiotic_path, file):
    """Shared loader used by the main filtering loop."""
    mapped = pd.read_csv(os.path.join(antibiotic_path, file))

    required_cols = ["GeneA_Symbol", "GeneB_Symbol", "Score"]
    if not all(col in mapped.columns for col in required_cols):
        return None

    mapped = mapped.dropna(subset=['GeneA_Symbol', 'GeneB_Symbol', 'Score'])
    mapped = mapped[mapped["GeneA_Symbol"] != mapped["GeneB_Symbol"]].copy()
    if len(mapped) == 0:
        return None

    mapped['GeneA_Symbol'] = mapped['GeneA_Symbol'].astype(str).str.lower().str.strip()
    mapped['GeneB_Symbol'] = mapped['GeneB_Symbol'].astype(str).str.lower().str.strip()
    mapped['Score'] = pd.to_numeric(mapped['Score'], errors='coerce')
    mapped = mapped.dropna(subset=['Score'])
    if len(mapped) == 0:
        return None

    return mapped.sort_values("Score", ascending=False)


# PLOTTING FUNCTIONS


def plot_network_summary(network_summary):
    fig, ax = plt.subplots(figsize=(10, 6))
    pivot_data = network_summary.pivot(index='Method', columns='Antibiotic', values='Predicted_Edges')
    pivot_data = pivot_data.fillna(0)
    pivot_data.plot(kind='bar', ax=ax, width=0.8, edgecolor='black')

    ax.set_xlabel('Method', fontsize=12, fontweight='bold')
    ax.set_ylabel('Number of Edges in Final Network', fontsize=12, fontweight='bold')
    ax.set_title('Final High-Confidence Network Sizes After Literature-Recommended Thresholds',
                 fontsize=14, fontweight='bold')
    ax.legend(title='Antibiotic', bbox_to_anchor=(1.05, 1), loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    plt.savefig('network_sizes.png', dpi=300, bbox_inches='tight')
    plt.show()



# MAIN FILTERING LOOP (literature default thresholds only)


def run_filtering(network_root, LITERATURE_THRESHOLDS_LOWER, output_root):
    """
    Applies the literature threshold once per Antibiotic x Method to build
    the high-confidence edge subset. No gold-standard evaluation.

    Saves output mirroring the input layout:
        high_confidence_output/<antibiotic>/<method>.csv
    i.e. one file per method inside each antibiotic's folder, same as the
    input networks/<antibiotic>/<method>.csv structure - five method files
    per antibiotic folder.
    """
    network_summary_rows = []
    high_confidence_edges = []

    for antibiotic in sorted(os.listdir(network_root)):
        antibiotic_path = os.path.join(network_root, antibiotic)
        if not os.path.isdir(antibiotic_path):
            continue

        print(f"\n{'='*60}\nProcessing Antibiotic: {antibiotic}\n{'='*60}")

        for file in sorted(os.listdir(antibiotic_path)):
            if not file.endswith(".csv"):
                continue

            method_raw = file.replace(".csv", "")
            method_key = method_raw.lower()
            if method_key not in LITERATURE_THRESHOLDS_LOWER:
                continue

            method = method_raw.upper()
            config = LITERATURE_THRESHOLDS_LOWER[method_key]
            print(f"\n  Method: {method}")

            try:
                mapped_sorted = _load_and_clean_network_file(antibiotic_path, file)
            except Exception as e:
                print(f"     Error reading file: {e}")
                continue
            if mapped_sorted is None:
                print(f"     Skipping: missing columns or no edges after cleaning")
                continue

            total_edges = len(mapped_sorted)
            print(f"    Total edges: {total_edges:,}")

            # Apply the literature threshold per Antibiotic x Method to
            # build the high-confidence edge list to export.
            subset, threshold_label, threshold_val = apply_literature_threshold(
                mapped_sorted, method, total_edges, config
            )
            print(f"    After {threshold_label}: {len(subset):,} edges")

            hc_subset = subset.copy()
            hc_subset['Antibiotic'] = antibiotic
            hc_subset['Method'] = method
            hc_subset['Threshold_Applied'] = threshold_label
            hc_subset = hc_subset[['Antibiotic', 'Method', 'GeneA_Symbol', 'GeneB_Symbol',
                                    'Score', 'Threshold_Applied']]
            high_confidence_edges.append(hc_subset)

            # Save per-antibiotic/per-method file, mirroring input structure:
            # high_confidence_output/<antibiotic>/<method_raw>.csv
            out_dir = os.path.join(output_root, antibiotic)
            os.makedirs(out_dir, exist_ok=True)
            out_path = os.path.join(out_dir, f"{method_raw}.csv")
            hc_subset.drop_duplicates(
                subset=['GeneA_Symbol', 'GeneB_Symbol']
            ).to_csv(out_path, index=False)
            print(f"    Saved: {out_path}")

            network_summary_rows.append({
                'Antibiotic': antibiotic,
                'Method': method,
                'Threshold_Applied': threshold_label,
                'Threshold_Value': threshold_val,
                'Threshold_Type': config['type'],
                'Total_Edges': total_edges,
                'Predicted_Edges': len(subset),
            })

    hc_df = pd.concat(high_confidence_edges, ignore_index=True) if high_confidence_edges else pd.DataFrame()
    network_summary_df = pd.DataFrame(network_summary_rows)
    return network_summary_df, hc_df



# MAIN EXECUTION


if __name__ == "__main__":

    print("="*80)
    print("GRN VALIDATION PIPELINE - PUBLICATION READY (v5)")
    print("="*80)
    print("\nKey Features:")
    print("  - Literature-recommended default thresholds -> high-confidence")
    print("    edge filtering (no cutoff optimization / sensitivity sweep,")
    print("    no gold-standard evaluation)")
    print("="*80)

    if not os.path.exists(network_root):
        print(f"\nError: Network root directory '{network_root}' not found!")
        raise SystemExit(1)

    #  FILTERING (literature default thresholds only)
    print(f"\nRunning threshold filtering...")
    network_summary, high_confidence_df = run_filtering(
        network_root, LITERATURE_THRESHOLDS_LOWER, output_root
    )

    if len(network_summary) == 0:
        print("\nNo results generated! Check your input files.")
        raise SystemExit(1)

    network_summary.to_csv("final_networks_summary.csv", index=False)
    print(f"\nSummary saved to: final_networks_summary.csv")

    # HIGH-CONFIDENCE EDGES (literature thresholds only) 
    # Per-antibiotic/per-method files are already saved inside run_filtering()
    # under high_confidence_output/<antibiotic>/<method>.csv (mirrors input
    # layout). A single combined CSV is also kept here as a convenience
    # rollup across all antibiotics/methods.
    if len(high_confidence_df) > 0:
        high_confidence_df = high_confidence_df.drop_duplicates(
            subset=['Antibiotic', 'Method', 'GeneA_Symbol', 'GeneB_Symbol']
        )
        high_confidence_df.to_csv("high_confidence_edges_combined.csv", index=False)
        print(f"\nCombined high-confidence edges (all antibiotics/methods) saved to: "
              f"high_confidence_edges_combined.csv ({len(high_confidence_df):,} edges)")
        print(f"Per-antibiotic/per-method files saved under: {output_root}/<antibiotic>/<method>.csv")
    else:
        print("\nNo high-confidence edges generated - check your input files.")

    # ========== PLOTS ==========
    print(f"\nGenerating figures...")
    try:
        plot_network_summary(network_summary)
        print(f"   Saved: network_sizes.png")
    except Exception as e:
        print(f"   Error: {e}")

    print("\n" + "="*80)
    print("ANALYSIS COMPLETE")
    print("="*80)
    print("\nOUTPUT FILES:")
    print("   1. final_networks_summary.csv                 - edge counts per method/antibiotic after literature thresholds")
    print(f"   2. {output_root}/<antibiotic>/<method>.csv       - filtered high-confidence edges, mirrors input structure (5 method files per antibiotic)")
    print("   3. high_confidence_edges_combined.csv         - convenience rollup of all antibiotics/methods in one file")
    print("   4. network_sizes.png                          - figure")