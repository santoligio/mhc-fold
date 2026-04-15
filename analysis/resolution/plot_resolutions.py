#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams["axes.labelsize"] = 16

# ------------------------------------------------------------
# Input
# ------------------------------------------------------------
CSV_FILE = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/filtered/pdb_mhc_annotations_filtered.csv"

# ------------------------------------------------------------
# Load and clean data
# ------------------------------------------------------------
df = pd.read_csv(CSV_FILE)

# Keep only necessary columns
df = df[["resolution", "mapped_class", "mapped_superclass"]]

# Remove missing resolution values
df = df.dropna(subset=["resolution"])

# Ensure numeric
df["resolution"] = df["resolution"].astype(str).str.replace(",", ".", regex=False)
df["resolution"] = pd.to_numeric(df["resolution"], errors="coerce")
df = df.dropna(subset=["resolution"])

# ------------------------------------------------------------
# Plotting function using seaborn
# ------------------------------------------------------------
def plot_violin_resolution(df):
    
    # Fixed class order
    ordered_classes = [
        "HLA Class Ia",
        "HLA Class Ia SC",
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
        "OMCP",
        "2L"
    ]

    # Keep only classes present in dataframe (preserves order)
    ordered_classes = [c for c in ordered_classes if c in df["mapped_class"].unique()]

    plt.figure(figsize=(14, 6))
    sns.set(style="white")

    # Violin plot
    ax = sns.violinplot(
        x="mapped_class",
        y="resolution",
        data=df,
        order=ordered_classes,
        inner="quartile",
        scale="width",
        palette="Pastel1"
    )

    # Overlay jittered points
    sns.stripplot(
        x="mapped_class",
        y="resolution",
        data=df,
        order=ordered_classes,
        color="k",
        size=3,
        alpha=0.5,
        jitter=True
    )

    ax.grid(False)

    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")
    plt.ylabel("Resolution (Å)")
    plt.xlabel("Class")
    
    # No title (removed as requested)

    plt.tight_layout()

    # Save figure
    filename = "pdb_res_distribution.png"
    plt.savefig(f"/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/resolution/{filename}", dpi=300)
    plt.close()


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plot_violin_resolution(df)
