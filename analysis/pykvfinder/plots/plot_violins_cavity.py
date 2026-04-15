#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Input files (Cavities)
# ------------------------------------------------------------
CSV_PDB_CAV  = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/pykvfinder/pdb/pdb_cavities/filter_cavities_out.csv"
CSV_AFDB_CAV = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/pykvfinder/afdb/afdb_cavities/filter_cavities_out.csv"

CSV_PDB_ANNOT  = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/filtered/pdb_mhc_annotations_filtered.csv"
CSV_AFDB_ANNOT = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/new_filters/afdb_mhc_annotations_reviewed.csv"

OUTDIR = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/pykvfinder"

# ------------------------------------------------------------
# Ordered classes
# ------------------------------------------------------------
ordered_classes = [
    "HLA Class Ia", "HLA Class Ia SC", "HLA Class Ib",
    "MR1", "MR1 SC",
    "CD1", "CD1 SC", "CD1 chimera",
    "EPCR", "ZAG", "HFE", "FcRn",
    "MIC", "ULBP",
    "UL18", "OMCP", "2L"
]

# ------------------------------------------------------------
# Manual class colors
# ------------------------------------------------------------
CLASS_COLORS = {
    "HLA Class Ia": "#3182bd",
    "HLA Class Ia SC": "#9ecae1",
    "HLA Class Ib": "#756bb1",
    "MR1": "#fd8d3c",
    "MR1 SC": "#fdae6b",
    "CD1": "#fd8d3c",
    "CD1 SC": "#fdae6b",
    "CD1 chimera": "#fdd0a2",
    "EPCR": "#fd8d3c",
    "ZAG": "#fd8d3c",
    "HFE": "#fd8d3c",
    "FcRn": "#fd8d3c",
    "MIC": "#fd8d3c",
    "ULBP": "#fd8d3c",
    "UL18": "#74c476",
    "OMCP": "#74c476",
    "2L": "#74c476"
}

# ------------------------------------------------------------
# Load cavity data
# ------------------------------------------------------------
df_pdb = pd.read_csv(CSV_PDB_CAV)
df_afdb = pd.read_csv(CSV_AFDB_CAV)

# ------------------------------------------------------------
# Collapse duplicated structures
# ------------------------------------------------------------
df_pdb = (
    df_pdb
    .groupby("structure", as_index=False)
    .agg(volume=("volume", "sum"),
         avg_depth=("avg_depth", "mean"))
)

df_afdb = (
    df_afdb
    .groupby("structure", as_index=False)
    .agg(volume=("volume", "sum"),
         avg_depth=("avg_depth", "mean"))
)

# ------------------------------------------------------------
# Extract IDs
# ------------------------------------------------------------
df_pdb["pdb_id"] = df_pdb["structure"].str[:4].str.upper()

def extract_uniprot(name):
    prefix = str(name).split("_")[0]
    return prefix.split("-")[1].upper()

df_afdb["uniprot_id"] = df_afdb["structure"].apply(extract_uniprot)

# ------------------------------------------------------------
# Load annotations
# ------------------------------------------------------------
df_pdb_annot  = pd.read_csv(CSV_PDB_ANNOT)
df_afdb_annot = pd.read_csv(CSV_AFDB_ANNOT)

df_pdb_annot["pdb_id"] = df_pdb_annot["pdb_id"].astype(str).str.upper()
df_afdb_annot["uniprot_id"] = df_afdb_annot["uniprot_id"].astype(str).str.upper()

# ------------------------------------------------------------
# Merge annotations
# ------------------------------------------------------------
df_pdb = pd.merge(
    df_pdb,
    df_pdb_annot[["pdb_id", "mapped_class"]],
    on="pdb_id",
    how="left"
)
df_pdb["source"] = "PDB"

df_afdb = pd.merge(
    df_afdb,
    df_afdb_annot[["uniprot_id", "mapped_class"]],
    on="uniprot_id",
    how="left"
)
df_afdb["source"] = "AFDB"

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------
df = pd.concat([df_pdb, df_afdb], ignore_index=True)
df = df.dropna(subset=["volume", "avg_depth", "mapped_class"])
df = df[df["mapped_class"].isin(ordered_classes)].copy()

# Keep only present classes
present_classes = [
    cls for cls in ordered_classes
    if cls in df["mapped_class"].unique()
]

df["mapped_class"] = pd.Categorical(
    df["mapped_class"],
    categories=present_classes,
    ordered=True
)

palette_present = [CLASS_COLORS[c] for c in present_classes]

sns.set_style("white")

# ------------------------------------------------------------
# Violin 1: Volume
# ------------------------------------------------------------
plt.figure(figsize=(14, 6))

sns.violinplot(
    data=df,
    x="mapped_class",
    y="volume",
    order=present_classes,
    palette=palette_present,
    inner="box",
    cut=0,
    linewidth=1
)

sns.stripplot(
    data=df[df["source"] == "PDB"],
    x="mapped_class",
    y="volume",
    order=present_classes,
    color="black",
    marker="o",
    size=3,
    jitter=0.25,
    alpha=0.5
)

sns.stripplot(
    data=df[df["source"] == "AFDB"],
    x="mapped_class",
    y="volume",
    order=present_classes,
    color="black",
    marker="^",
    size=3,
    jitter=0.25,
    alpha=0.5
)

plt.xticks(rotation=60, ha="right")
plt.xlabel("Class")
plt.ylabel("Cavity Volume (Å³)")
sns.despine()
plt.tight_layout()

plt.savefig(f"{OUTDIR}/violin_cavity_volume_PDB_AFDB.png",
            dpi=300, bbox_inches="tight")
plt.close()

# ------------------------------------------------------------
# Violin 2: Average Depth
# ------------------------------------------------------------
plt.figure(figsize=(14, 6))

sns.violinplot(
    data=df,
    x="mapped_class",
    y="avg_depth",
    order=present_classes,
    palette=palette_present,
    inner="box",
    cut=0,
    linewidth=1
)

sns.stripplot(
    data=df[df["source"] == "PDB"],
    x="mapped_class",
    y="avg_depth",
    order=present_classes,
    color="black",
    marker="o",
    size=3,
    jitter=0.25,
    alpha=0.5
)

sns.stripplot(
    data=df[df["source"] == "AFDB"],
    x="mapped_class",
    y="avg_depth",
    order=present_classes,
    color="black",
    marker="^",
    size=3,
    jitter=0.25,
    alpha=0.5
)

plt.xticks(rotation=60, ha="right")
plt.xlabel("Class")
plt.ylabel("Average Cavity Depth (Å)")
sns.despine()
plt.tight_layout()

plt.savefig(f"{OUTDIR}/violin_cavity_avg_depth_PDB_AFDB.png",
            dpi=300, bbox_inches="tight")
plt.close()

print("Cavity violin plots (Volume and Avg Depth) saved successfully.")
