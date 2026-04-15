#!/usr/bin/env python3

import pandas as pd

# =========================
# INPUT FILES
# =========================
ANNOTATIONS_CSV = "binders_annotations_merged.csv"
GENE_MAPPING_CSV = "binders_mapping_gene_name.csv"
GO_MAPPING_CSV = "binders_mapping_go_annotations.csv"
CHAINS_ANALYSIS_CSV = "binders_chains_analysis_merged.csv"

# =========================
# OUTPUT FILES
# =========================
OUTPUT_ANNOTATIONS_CSV = "binders_annotations_filtered.csv"
OUTPUT_CHAINS_CSV = "binders_chains_analysis_filtered.csv"

# =========================
# LOAD FILES
# =========================
df_annotations = pd.read_csv(ANNOTATIONS_CSV)
df_gene_map = pd.read_csv(GENE_MAPPING_CSV)
df_go_map = pd.read_csv(GO_MAPPING_CSV)
df_chains = pd.read_csv(CHAINS_ANALYSIS_CSV)

# =========================
# Prepare mapping dictionaries
# =========================
gene_map_dict = df_gene_map.set_index("original_gene").to_dict("index")
go_map_dict = df_go_map.set_index("binder_name_raw").to_dict("index")

# =========================
# Process rows
# =========================
output_rows = []

for idx, row in df_annotations.iterrows():

    gene_name = str(row["gene_name"]).strip()

    # Exclude B2M
    if gene_name == "B2M":
        continue

    classification = str(row["classification"]).strip()

    mapped_class = None
    mapped_superclass = None

    # -------------------------
    # CASE 1: gene_name present
    # -------------------------
    if gene_name != "Missing entry":

        if classification in gene_map_dict:
            mapped_class = gene_map_dict[classification]["mapped_class"]
            mapped_superclass = gene_map_dict[classification]["mapped_superclass"]
        else:
            print(f"WARNING: No gene mapping found for classification '{classification}'"
                  f"(PDB: {row['pdb_id']}, chain: {row['new_chain']})")

    # -------------------------
    # CASE 2: gene_name missing
    # -------------------------
    else:

        if classification in go_map_dict:
            regex_mapping = go_map_dict[classification]["regex_mapping"]
            manual_mapping = go_map_dict[classification]["manual_mapping"]

            if str(regex_mapping).lower() != "unknown":
                mapped_class = regex_mapping
            else:
                mapped_class = manual_mapping

            mapped_superclass = go_map_dict[classification]["mapped_superclass"]
        else:
            print(f"WARNING: No GO mapping found for classification '{classification}'"
                  f"(PDB: {row['pdb_id']}, chain: {row['new_chain']})")

    # Append row only if mapping was found
    if mapped_class is not None and mapped_superclass is not None:
        output_rows.append({
            "pdb_id": row["pdb_id"],
            "new_chain": row["new_chain"],
            "original_chain": row["original_chain"],
            "uniprot_id": row["uniprot_id"],
            "organism": row["organism"],
            "classification": classification,
            "mapped_class": mapped_class,
            "mapped_superclass": mapped_superclass
        })

# =========================
# Save filtered annotations
# =========================
df_output = pd.DataFrame(output_rows)
df_output.to_csv(OUTPUT_ANNOTATIONS_CSV, index=False)

print("\nFiltering complete (annotations).")
print(f"Output written to: {OUTPUT_ANNOTATIONS_CSV}")
print(f"Total entries kept: {len(df_output)}")

# =========================
# Filter chains analysis
# =========================

# Get surviving pdb_ids
valid_pdb_ids = set(df_output["pdb_id"])

# Keep only matching pdb_ids
df_chains_filtered = df_chains.merge(
    df_output[["pdb_id", "new_chain"]],
    on=["pdb_id", "new_chain"],
    how="inner"
)

df_chains_filtered.to_csv(OUTPUT_CHAINS_CSV, index=False)

print("\nChains analysis filtering complete.")
print(f"Output written to: {OUTPUT_CHAINS_CSV}")
print(f"Total chains entries kept: {len(df_chains_filtered)}")
