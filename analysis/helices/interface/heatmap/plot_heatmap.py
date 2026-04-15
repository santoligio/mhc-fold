#!/usr/bin/env python3

import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# ------------------------------------------------------------
# Input file
# ------------------------------------------------------------
TXT_FILE = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/helix_statistics_afdb_class.txt"

# ------------------------------------------------------------
# All possible amino acids
# ------------------------------------------------------------
AA_FULLNAME = {
    'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine',
    'ASP': 'Aspartic acid', 'CYS': 'Cysteine', 'GLN': 'Glutamine',
    'GLU': 'Glutamic acid', 'GLY': 'Glycine', 'HIS': 'Histidine',
    'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine',
    'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline',
    'SER': 'Serine', 'THR': 'Threonine', 'TRP': 'Tryptophan',
    'TYR': 'Tyrosine', 'VAL': 'Valine'
}

# ------------------------------------------------------------
# Dayhoff order
# ------------------------------------------------------------
AA_ORDER = [
    # Sulfur polymerization
    'CYS',

    # Small
    'GLY', 'SER', 'THR', 'ALA', 'PRO',

    # Acid and amide
    'ASP', 'GLU', 'ASN', 'GLN',

    # Basic
    'ARG', 'HIS', 'LYS',

    # Hydrophobic
    'LEU', 'VAL', 'MET', 'ILE',

    # Aromatic
    'TYR', 'PHE', 'TRP'
]

# ------------------------------------------------------------
# Parse file
# ------------------------------------------------------------
data = {}
current_group = None

with open(TXT_FILE, "r") as f:
    for line in f:
        line = line.strip()

        # Detect group
        if line.startswith("Group:"):
            current_group = line.replace("Group:", "").strip()
            data[current_group] = {aa: 0.0 for aa in AA_ORDER}

        # Match amino acid lines
        match = re.match(r".+\((\w{3})\)\s*:\s*([\d\.]+)%", line)
        if match and current_group:
            aa = match.group(1)
            value = float(match.group(2))
            if aa in AA_ORDER:
                data[current_group][aa] = value

# ------------------------------------------------------------
# Create DataFrame
# ------------------------------------------------------------
df = pd.DataFrame.from_dict(data, orient="index")
df = df[AA_ORDER]  # ensure consistent column order

# ------------------------------------------------------------
# Fix class order (row order)
# ------------------------------------------------------------
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
    "OMPC",
    "2L"
]

# Keep only classes present in dataframe (preserve order)
ordered_classes = [c for c in ordered_classes if c in df.index]

df = df.loc[ordered_classes]

# ------------------------------------------------------------
# Plot heatmap
# ------------------------------------------------------------
plt.figure(figsize=(14, 6))

sns.heatmap(
    df,
    cmap="viridis",
    linewidths=0.5,
    linecolor="white",
    cbar_kws={"label": "Percentage (%)"}
)

# ------------------------------------------------------------
# Color amino acid labels by Dayhoff group
# ------------------------------------------------------------
group_colors = {
    'CYS': '#C75C1A',

    'GLY': '#1B5E20', 'SER': '#1B5E20', 'THR': '#1B5E20',
    'ALA': '#1B5E20', 'PRO': '#1B5E20',

    'ASP': '#8B1E3F', 'GLU': '#8B1E3F', 'ASN': '#8B1E3F', 'GLN': '#8B1E3F',

    'ARG': '#005F73', 'HIS': '#005F73', 'LYS': '#005F73',

    'LEU': '#2F2F2F', 'VAL': '#2F2F2F', 'MET': '#2F2F2F', 'ILE': '#2F2F2F',

    'TYR': '#8C6D1F', 'PHE': '#8C6D1F', 'TRP': '#8C6D1F'
}

ax = plt.gca()
for tick_label in ax.get_xticklabels():
    aa = tick_label.get_text()
    if aa in group_colors:
        tick_label.set_color(group_colors[aa])

plt.title("Amino Acid Composition per Class (AFDB)")

plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)

plt.tight_layout()
plt.savefig(
    "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/heatmap_afdb.png",
    dpi=300,
    bbox_inches="tight"   # important so legend is not cut
)

plt.close()
