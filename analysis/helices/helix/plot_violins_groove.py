#!/usr/bin/env python3

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
CSV_PDB_MAIN   = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/pdb_info_helices.csv"
CSV_AFDB_MAIN  = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/afdb_info_helices.csv"

CSV_PDB_ANNOT  = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/filtered/pdb_mhc_annotations_filtered.csv"
CSV_AFDB_ANNOT = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/new_filters/afdb_mhc_annotations_reviewed.csv"

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
# Load data
# ------------------------------------------------------------
df_pdb_main   = pd.read_csv(CSV_PDB_MAIN)
df_afdb_main  = pd.read_csv(CSV_AFDB_MAIN)
df_pdb_annot  = pd.read_csv(CSV_PDB_ANNOT)
df_afdb_annot = pd.read_csv(CSV_AFDB_ANNOT)

# ---- PDB MERGE ----
df_pdb_main["pdb_id"] = df_pdb_main["file"].str[:4].str.upper()
df_pdb_annot["pdb_id"] = df_pdb_annot["pdb_id"].astype(str).str.upper()

df_pdb = pd.merge(
    df_pdb_main,
    df_pdb_annot[["pdb_id", "mapped_class"]],
    on="pdb_id",
    how="left"
)
df_pdb["source"] = "PDB"

# ---- AFDB MERGE ----
def extract_uniprot(filename):
    prefix = str(filename).split("_")[0]
    return prefix.split("-")[1].upper()

df_afdb_main["uniprot_id"] = df_afdb_main["file"].apply(extract_uniprot)
df_afdb_annot["uniprot_id"] = df_afdb_annot["uniprot_id"].astype(str).str.upper()

df_afdb = pd.merge(
    df_afdb_main,
    df_afdb_annot[["uniprot_id", "mapped_class"]],
    on="uniprot_id",
    how="left"
)
df_afdb["source"] = "AFDB"

# ------------------------------------------------------------
# Combine both sources
# ------------------------------------------------------------
df = pd.concat([df_pdb, df_afdb], ignore_index=True)
df = df.dropna(subset=["helix_pct_total", "helix_com_dist", "mapped_class"])
df = df[df["mapped_class"].isin(ordered_classes)].copy()

# ---- Keep only classes that are present ----
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

# ------------------------------------------------------------
# Plot style
# ------------------------------------------------------------
sns.set_style("white")

# ------------------------------------------------------------
# Violin 1: Helix COM Distance
# ------------------------------------------------------------
plt.figure(figsize=(14, 6))

sns.violinplot(
    data=df,
    x="mapped_class",
    y="helix_com_dist",
    order=present_classes,
    palette=palette_present,
    inner="box",
    cut=0,
    linewidth=1
)

# PDB circles
sns.stripplot(
    data=df[df["source"] == "PDB"],
    x="mapped_class",
    y="helix_com_dist",
    order=present_classes,
    color="black",
    marker="o",
    size=3,
    jitter=0.25,
    alpha=0.5
)

# AFDB triangles
sns.stripplot(
    data=df[df["source"] == "AFDB"],
    x="mapped_class",
    y="helix_com_dist",
    order=present_classes,
    color="black",
    marker="^",
    size=3,
    jitter=0.25,
    alpha=0.5
)

plt.xticks(rotation=60, ha="right")
plt.xlabel("Class")
plt.ylabel("Helix COM Distance (Å)")
sns.despine()
plt.tight_layout()

plt.savefig(
    "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/violin_helix_com_dist_PDB_AFDB.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# Violin 2: Helix Percentage Total
# ------------------------------------------------------------
plt.figure(figsize=(14, 6))

sns.violinplot(
    data=df,
    x="mapped_class",
    y="helix_pct_total",
    order=present_classes,
    palette=palette_present,
    inner="box",
    cut=0,
    linewidth=1
)

sns.stripplot(
    data=df[df["source"] == "PDB"],
    x="mapped_class",
    y="helix_pct_total",
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
    y="helix_pct_total",
    order=present_classes,
    color="black",
    marker="^",
    size=3,
    jitter=0.25,
    alpha=0.5
)

plt.xticks(rotation=60, ha="right")
plt.xlabel("Class")
plt.ylabel("Helical Residues (%)")
sns.despine()
plt.tight_layout()

plt.savefig(
    "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/violin_helix_pct_total_PDB_AFDB.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

print("Unified violin plots with PDB (circles) and AFDB (triangles) saved successfully.")