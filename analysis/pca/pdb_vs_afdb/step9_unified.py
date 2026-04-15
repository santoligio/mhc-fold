#!/usr/bin/env python3

# ---------------------
# PCA for TM-align all-vs-all (chunked output) — AFDB + PDB
# ---------------------

import pandas as pd
import numpy as np
import os
import glob
from sklearn.decomposition import PCA

# ------------------------------------------------------------
# Inputs / Outputs
# ------------------------------------------------------------
TMALIGN_DIR = "/mnt/4TB/giovanna/foldseek/version_02/analysis/tmalign/pdb_vs_afdb"
TMALIGN_PATTERN = "pdb_vs_afdb_part*.csv"

OUTPUT_PCA = "/mnt/4TB/giovanna/foldseek/version_02/analysis/pca/pdb_vs_afdb/pdb_vs_afdb_pca_result_mhc.csv"

# ------------------------------------------------------------
# Helper: extract ID from filename (AFDB or PDB)
# ------------------------------------------------------------
def extract_structure_id(filename):
    """
    Handles both:

    AFDB:
        AF-A0A060CZ20_trimmed_mhc.pdb  -> A0A060CZ20

    PDB:
        3mre_trimmed.pdb              -> 3MRE
        3mre_A.pdb                    -> 3MRE
    """
    if pd.isna(filename):
        return None

    name = os.path.splitext(os.path.basename(filename))[0]

    # ---------------- AFDB ----------------
    if name.startswith("AF-"):
        name = name.split("_")[0]      # AF-A0A060CZ20
        parts = name.split("-")
        if len(parts) >= 2:
            return parts[1].upper()
        return None

    # ---------------- PDB ----------------
    # Take first 4 characters as PDB ID
    pdb_id = name[:4]

    if len(pdb_id) == 4:
        return pdb_id.upper()

    return None


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
# Extract structure IDs (AFDB + PDB)
# ------------------------------------------------------------
read_data["id1"] = read_data["file1"].apply(extract_structure_id)
read_data["id2"] = read_data["file2"].apply(extract_structure_id)

# Drop invalid rows
read_data = read_data.dropna(subset=["id1", "id2", "tm_score_avg"])

# ------------------------------------------------------------
# Collect unique structures
# ------------------------------------------------------------
all_ids = pd.unique(
    read_data[["id1", "id2"]].values.ravel("K")
)

all_ids = sorted(all_ids)
n = len(all_ids)

print(f"Found {n} unique structures")

if n == 0:
    raise RuntimeError("No valid structures found after ID extraction.")

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
    "UniProt_ID": all_ids,   # kept same column name for compatibility
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
})

pca_df.to_csv(OUTPUT_PCA, index=False)

print("PCA finished")
print(f"Explained variance: {pca.explained_variance_ratio_}")
print(f"Saved PCA results to {OUTPUT_PCA}")
