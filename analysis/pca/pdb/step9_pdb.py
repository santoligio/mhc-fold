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
TMALIGN_DIR = "/mnt/4TB/giovanna/foldseek/version_02/analysis/tmalign/filtered_pdb_avg"
TMALIGN_PATTERN = "mhc_all-vs-all_part*.csv"

OUTPUT_PCA = "/mnt/4TB/giovanna/foldseek/version_02/analysis/pca/teste_avg/pdb_pca_result_mhc.csv"

# ------------------------------------------------------------
# Helper: normalize PDB ID (MUST match TM-align logic)
# ------------------------------------------------------------
def extract_base_pdb(filename):
    """
    Examples:
    1abc-assembly1.pdb -> 1ABC
    1abc_A.pdb         -> 1ABC
    """
    if pd.isna(filename):
        return None
    return str(filename).split("-")[0][:4].upper()

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
# Extract PDB IDs (consistent with TM-align)
# ------------------------------------------------------------
read_data["pdb1"] = read_data["file1"].apply(extract_base_pdb)
read_data["pdb2"] = read_data["file2"].apply(extract_base_pdb)

# Drop invalid rows
read_data = read_data.dropna(subset=["pdb1", "pdb2", "tm_score_avg"])

# ------------------------------------------------------------
# Collect unique PDBs
# ------------------------------------------------------------
all_pdbs = pd.unique(
    read_data[["pdb1", "pdb2"]].values.ravel("K")
)

all_pdbs = sorted(all_pdbs)
n = len(all_pdbs)

print(f"Found {n} unique PDB structures")

pdb_index = {pdb: i for i, pdb in enumerate(all_pdbs)}

# ------------------------------------------------------------
# Build distance matrix (1 - TM-score)
# ------------------------------------------------------------
feat_matrix = np.ones((n, n), dtype=np.float32)
np.fill_diagonal(feat_matrix, 0.0)

for _, row in read_data.iterrows():
    i = pdb_index[row["pdb1"]]
    j = pdb_index[row["pdb2"]]
    dist = 1.0 - row["tm_score_avg"]

    feat_matrix[i, j] = dist
    feat_matrix[j, i] = dist

# ------------------------------------------------------------
# PCA
# ------------------------------------------------------------
pca = PCA(n_components=2)
pca_result = pca.fit_transform(feat_matrix)

pca_df = pd.DataFrame({
    "PDB": all_pdbs,
    "PC1": pca_result[:, 0],
    "PC2": pca_result[:, 1],
})

pca_df.to_csv(OUTPUT_PCA, index=False)

print("PCA finished")
print(f"Explained variance: {pca.explained_variance_ratio_}")
print(f"Saved PCA results to {OUTPUT_PCA}")
