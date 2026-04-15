#!/usr/bin/env python3

# ---------------------
# PCA for TM-align all-vs-all (chunked output) — AFDB
# ---------------------

import pandas as pd
import numpy as np
import os
import glob
from sklearn.decomposition import PCA

# ------------------------------------------------------------
# Inputs / Outputs
# ------------------------------------------------------------
TMALIGN_DIR = "/mnt/4TB/giovanna/foldseek/version_02/analysis/tmalign/reviewed+length"
TMALIGN_PATTERN = "afdb_vs_afdb_part*.csv"

OUTPUT_PCA = "/mnt/4TB/giovanna/foldseek/version_02/analysis/pca/reviewed+length/afdb2_pca_result_mhc.csv"

# ------------------------------------------------------------
# Helper: extract UniProt ID from AFDB filename
# (MUST match TM-align AFDB logic)
# ------------------------------------------------------------
def extract_uniprot_id(filename):
    """
    Example:
    AF-A0A060CZ20_trimmed_mhc.pdb -> A0A060CZ20
    """
    if pd.isna(filename):
        return None

    name = os.path.splitext(filename)[0]
    name = name.split("_")[0]          # AF-A0A060CZ20
    parts = name.split("-")

    if len(parts) < 2:
        return None

    return parts[1].upper()


# ------------------------------------------------------------
# Load all TM-align chunks
# ------------------------------------------------------------
csv_files = sorted(
    glob.glob(os.path.join(TMALIGN_DIR, TMALIGN_PATTERN))
)

if not csv_files:
    raise RuntimeError("No TM-align chunk CSV files found")

print(f"Loading {len(csv_files)} TM-align chunks")

df_list = []
for f in csv_files:
    df = pd.read_csv(f)
    df_list.append(df)

read_data = pd.concat(df_list, ignore_index=True)

print(f"Loaded {len(read_data):,} pairwise comparisons")

# ------------------------------------------------------------
# Extract UniProt IDs (AFDB)
# ------------------------------------------------------------
read_data["id1"] = read_data["file1"].apply(extract_uniprot_id)
read_data["id2"] = read_data["file2"].apply(extract_uniprot_id)

# Drop invalid rows
read_data = read_data.dropna(subset=["id1", "id2", "tm_score_avg"])

# ------------------------------------------------------------
# Collect unique AFDB entries
# ------------------------------------------------------------
all_ids = pd.unique(
    read_data[["id1", "id2"]].values.ravel("K")
)

all_ids = sorted(all_ids)
n = len(all_ids)

print(f"Found {n} unique AFDB structures")

id_index = {uid: i for i, uid in enumerate(all_ids)}

# ------------------------------------------------------------
# Build distance matrix (1 - TM-score)
# ------------------------------------------------------------
feat_matrix = np.ones((n, n), dtype=np.float32)
np.fill_diagonal(feat_matrix, 0.0)

for _, row in read_data.iterrows():
    i = id_index[row["id1"]]
    j = id_index[row["id2"]]
    dist = 1.0 - row["tm_score_avg"]

    feat_matrix[i, j] = dist
    feat_matrix[j, i] = dist

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------
pca = PCA(n_components=2)
pca_result = pca.fit_transform(feat_matrix)

pca_df = pd.DataFrame({
    "UniProt_ID": all_ids,
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
})

pca_df.to_csv(OUTPUT_PCA, index=False)

print("PCA finished")
print(f"Explained variance: {pca.explained_variance_ratio_}")
print(f"Saved PCA results to {OUTPUT_PCA}")
