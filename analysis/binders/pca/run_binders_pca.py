# ---------------------
# PCA for TM-align all-vs-all (chunked output)
# ---------------------

import pandas as pd
import numpy as np
import os
import glob
from sklearn.decomposition import PCA

# ------------------------------------------------------------
# Inputs / Outputs
# ------------------------------------------------------------
TMALIGN_DIR = "/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/tmalign/"
TMALIGN_PATTERN = "binders_all-vs-all_part*.csv"

OUTPUT_PCA = "/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/pca/pca_result_binders.csv"

# ------------------------------------------------------------
# Helper: extract PDB + chain (binder-specific)
# ------------------------------------------------------------
def extract_pdb_chain(filename):
    """
    Example:
    1lp9_binder_chainD.pdb → 1LP9_D
    """
    if pd.isna(filename):
        return None

    name = os.path.splitext(filename)[0]
    parts = name.split("_")

    if len(parts) < 3:
        return None

    pdb_id = parts[0].upper()

    chain_part = parts[-1]  # chainD
    if not chain_part.lower().startswith("chain"):
        return None

    chain = chain_part.replace("chain", "").replace("CHAIN", "")

    return f"{pdb_id}_{chain}"


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
# Extract binder identifiers (PDB_CHAIN)
# ------------------------------------------------------------
read_data["pdb1"] = read_data["file1"].apply(extract_pdb_chain)
read_data["pdb2"] = read_data["file2"].apply(extract_pdb_chain)

# Drop invalid rows
read_data = read_data.dropna(subset=["pdb1", "pdb2", "tm_score_avg"])

# ------------------------------------------------------------
# Collect unique binder structures
# ------------------------------------------------------------
all_structures = pd.unique(
    read_data[["pdb1", "pdb2"]].values.ravel("K")
)

all_structures = sorted(all_structures)
n = len(all_structures)

print(f"Found {n} unique binder structures")

structure_index = {pdb: i for i, pdb in enumerate(all_structures)}

# ------------------------------------------------------------
# Build distance matrix (1 - TM-score)
# ------------------------------------------------------------
feat_matrix = np.ones((n, n), dtype=np.float32)
np.fill_diagonal(feat_matrix, 0.0)

for _, row in read_data.iterrows():
    i = structure_index[row["pdb1"]]
    j = structure_index[row["pdb2"]]
    dist = 1.0 - row["tm_score_avg"]

    feat_matrix[i, j] = dist
    feat_matrix[j, i] = dist

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------
pca = PCA(n_components=2)
pca_result = pca.fit_transform(feat_matrix)

pca_df = pd.DataFrame({
    "Binder": all_structures,
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
})

pca_df.to_csv(OUTPUT_PCA, index=False)

print("PCA finished")
print(f"Explained variance: {pca.explained_variance_ratio_}")
print(f"Saved PCA results to {OUTPUT_PCA}")
