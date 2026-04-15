#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["axes.labelsize"] = 16

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
CSV_REVIEWED = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/new_filters/afdb_mhc_annotations_reviewed.csv"
CSV_LENGTH = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/new_filters/afdb_mhc_annotations_length_175_185.csv"

# ------------------------------------------------------------
# Load and clean function
# ------------------------------------------------------------
def load_and_clean(csv_path):
    df = pd.read_csv(csv_path)
    df = df[["avg_pLDDT", "mapped_class", "mapped_superclass"]]
    df = df.dropna(subset=["avg_pLDDT"])

    df["avg_pLDDT"] = (
        df["avg_pLDDT"]
        .astype(str)
        .str.replace(",", ".", regex=False)
    )

    df["avg_pLDDT"] = pd.to_numeric(df["avg_pLDDT"], errors="coerce")
    df = df.dropna(subset=["avg_pLDDT"])

    return df

df_reviewed = load_and_clean(CSV_REVIEWED)
df_length = load_and_clean(CSV_LENGTH)

# Add a dataset column to distinguish
df_reviewed["Dataset"] = "Reviewed"
df_length["Dataset"] = "Length Filter"

# Combine for seaborn plotting
df_combined = pd.concat([df_reviewed, df_length], ignore_index=True)

# ------------------------------------------------------------
# Seaborn plot function (side-by-side violins with median lines)
# ------------------------------------------------------------
def plot_side_by_side_violin(df):
    plt.figure(figsize=(14, 6))
    sns.set(style="white")

    # Fixed class order (exactly as requested)
    ordered_classes = [
        "HLA Class Ia",
        "HLA Class IA SC",
        "HLA Class Ib",
        "MR1",
        "MR1 SC",
        "CD1",
        "CD1 SC",
        "CD1 chimera",
        "EPCR",
        "ZAG",
        "HFE",
        "FcRn",
        "MIC",
        "ULBP",
        "UL18",
        "OMPC",
        "2L"
    ]

    # Keep only classes present in dataframe (preserve order)
    ordered_classes = [c for c in ordered_classes if c in df["mapped_class"].unique()]

    # Violin plot
    ax = sns.violinplot(
        x="mapped_class",
        y="avg_pLDDT",
        hue="Dataset",
        data=df,
        order=ordered_classes,
        inner="quartile",
        scale="width",
        width=0.8,
        palette={"Reviewed": "#ff9999", "Length Filter": "#99ccff"},
        linewidth=1.5
    )

    # Overlay jittered points
    sns.stripplot(
        x="mapped_class",
        y="avg_pLDDT",
        hue="Dataset",
        data=df,
        dodge=True,
        order=ordered_classes,
        palette={"Reviewed": "#8b0000", "Length Filter": "#00008b"},
        size=4,
        alpha=0.5,
        jitter=True,
        ax=ax,
        linewidth=0
    )

    # Draw medians manually in gray
    classes = ordered_classes
    for i, cls in enumerate(classes):
        for dataset in ["Reviewed", "Length Filter"]:
            vals = df[(df["mapped_class"] == cls) & (df["Dataset"] == dataset)]["avg_pLDDT"]
            if not vals.empty:
                median_val = np.median(vals)
                offset = -0.2 if dataset == "Reviewed" else 0.2
                plt.hlines(
                    median_val,
                    i + offset - 0.2,
                    i + offset + 0.2,
                    colors="gray",
                    linewidth=1.5
                )

    # X-ticks
    ax.set_xticks(range(len(classes)))
    ax.set_xticklabels(classes, rotation=45, ha="right")

    # Remove duplicate legends
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys(), title="", loc="upper right")

    # Remove grid lines but keep ticks
    ax.grid(False)

    plt.ylabel("Average pLDDT")
    plt.xlabel("Class")

    # Title removed as requested

    plt.tight_layout()

    plt.savefig(
        "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/plddt/Average_pLDDT_distribution_between_AFDB_groups.png",
        dpi=300
    )
    plt.close()


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plot_side_by_side_violin(df_combined)
