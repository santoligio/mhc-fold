#!/usr/bin/env python3

import pandas as pd

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
PDB_ANN_FILE = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/pdb/mhc/pdb_mhc_annotations_edited.csv"
AFDB_ANN_FILE = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/afdb/afdb_mhc_annotations_edited.csv"

PDB_ASSEMBLIES_FILE = "/mnt/4TB/giovanna/foldseek/version_02/analysis/pdb_assemblies_analysis.csv"
AFDB_MODELS_FILE = "/mnt/4TB/giovanna/foldseek/version_02/analysis/afdb_models_analysis.csv"

AFDB_REPRESENTATIVES_FILE = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/filtered/representatives_afdb.csv"
GENE_MAPPING_FILE = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/filtered/gene_mapping.csv"
AFDB_UNIPROT_REMOVE_FILE = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/filtered/uniprot_to_remove/uniprot_proteins_to_remove.txt"

# ------------------------------------------------------------
# Output files
# ------------------------------------------------------------
PDB_ANN_OUT = "pdb_mhc_annotations_filtered.csv"
AFDB_ANN_OUT = "afdb_mhc_annotations_filtered.csv"

PDB_ASSEMBLIES_OUT = "pdb_assemblies_analysis_filtered.csv"
AFDB_MODELS_OUT = "afdb_models_analysis_filtered.csv"

REMOVED_LOG_FILE = "removed_ids.log"

# ------------------------------------------------------------
# Global filters
# ------------------------------------------------------------
ALLOWED_SPECIES = {
    "Homo sapiens",
    "Human cytomegalovirus (strain AD169)",
    "Cowpox virus (strain Brighton Red)",
    "Yaba-like disease virus",
}

EXCLUDED_GENE_NAMES = {
    "HLA",
    "HLA locus",
    "MHC",
    "HLA-DRB1",
    "HLA-DPB1",
    "DKFZp686N10220",
    "DKFZp686P19218",
    "FLJ45422",
    "LOC554223",
    "DKFZp762B162",
}

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
pdb_ann = pd.read_csv(PDB_ANN_FILE)
afdb_ann = pd.read_csv(AFDB_ANN_FILE)

pdb_assemblies = pd.read_csv(PDB_ASSEMBLIES_FILE)
afdb_models = pd.read_csv(AFDB_MODELS_FILE)

afdb_reps = pd.read_csv(AFDB_REPRESENTATIVES_FILE)
gene_mapping = pd.read_csv(GENE_MAPPING_FILE)

removed_log = []

def log_removal(pdb, uniprot_id, database, reason):
    removed_log.append(
        {
            "pdb": pdb,
            "uniprot_id": uniprot_id,
            "database": database,
            "reason": reason,
        }
    )

# ------------------------------------------------------------
# ------------------ PDB-SPECIFIC RULES ----------------------
# ------------------------------------------------------------

# 1. Apply author-assigned overrides
for idx, row in pdb_ann.iterrows():
    if pd.notna(row.get("gene_name_assigned")):
        pdb_ann.at[idx, "gene_name"] = row["gene_name_assigned"]

    if pd.notna(row.get("grouped_gene_assigned")):
        pdb_ann.at[idx, "gene_name"] = row["grouped_gene_assigned"]

    if pd.notna(row.get("organism_assigned")):
        pdb_ann.at[idx, "organism"] = row["organism_assigned"]

# 2. Resolve possibly chimeric mappings
if "possibly_chimeric" in pdb_ann.columns:
    chimeric = pdb_ann[pdb_ann["possibly_chimeric"] == "yes"]

    to_drop = []

    for (pdb_id, chain), group in chimeric.groupby(["pdb_id", "chain"]):
        if len(group) > 1:
            preferred = group[group["grouped_gene_assigned"].notna()]
            if not preferred.empty:
                dropped = group.index.difference(preferred.index)
                for i in dropped:
                    log_removal(
                        pdb=pdb_ann.at[i, "pdb_id"],
                        uniprot_id=pdb_ann.at[i, "uniprot_id"]
                        if "uniprot_id" in pdb_ann.columns
                        else None,
                        database="pdb",
                        reason="chimeric_mapping_resolution",
                    )
                to_drop.extend(dropped)

    pdb_ann = pdb_ann.drop(index=to_drop)

# ------------------------------------------------------------
# ------------------ AFDB-SPECIFIC RULES ----------------------
# ------------------------------------------------------------

#hla_reps = afdb_reps[afdb_reps["gene"].isin({"HLA-A", "HLA-B", "HLA-C"})]
#allowed_models = set(hla_reps["model"])
#
#mask_hla = afdb_ann["gene_name"].isin({"HLA-A", "HLA-B", "HLA-C"}) & ~afdb_ann[
#    "pdb_id"
#].isin(allowed_models)
#
#for _, row in afdb_ann[mask_hla].iterrows():
#    log_removal(
#        pdb=f"AF-{row['uniprot_id']}-model",
#        uniprot_id=row["uniprot_id"],
#        database="afdb",
#        reason="non_representative_hla_model",
#    )
#
#afdb_ann = afdb_ann[~mask_hla]

# ------------------------------------------------------------
# -------------------- GLOBAL FILTERS -------------------------
# ------------------------------------------------------------

def apply_global_filters(df, database):
    mask_species = ~df["organism"].isin(ALLOWED_SPECIES)
    for _, row in df[mask_species].iterrows():
        log_removal(
            pdb=row["pdb_id"],
            uniprot_id=row["uniprot_id"] if "uniprot_id" in df.columns else None,
            database=database,
            reason="excluded_species",
        )
    df = df[~mask_species]

    mask_gene = df["gene_name"].isin(EXCLUDED_GENE_NAMES)
    for _, row in df[mask_gene].iterrows():
        log_removal(
            pdb=row["pdb_id"],
            uniprot_id=row["uniprot_id"] if "uniprot_id" in df.columns else None,
            database=database,
            reason="excluded_gene_name",
        )
    df = df[~mask_gene]

    return df

pdb_ann = apply_global_filters(pdb_ann, "pdb")
afdb_ann = apply_global_filters(afdb_ann, "afdb")

# ------------------------------------------------------------
# -------- GLOBAL GENE NAME CORRECTION (DB-SPECIFIC) ---------
# ------------------------------------------------------------

pdb_gene_map = dict(
    zip(
        gene_mapping[gene_mapping["database"] == "pdb"]["gene_name"],
        gene_mapping[gene_mapping["database"] == "pdb"]["mapped_gene_name"],
    )
)

afdb_gene_map = dict(
    zip(
        gene_mapping[gene_mapping["database"] == "afdb"]["gene_name"],
        gene_mapping[gene_mapping["database"] == "afdb"]["mapped_gene_name"],
    )
)

pdb_ann["gene_name"] = pdb_ann["gene_name"].replace(pdb_gene_map)
afdb_ann["gene_name"] = afdb_ann["gene_name"].replace(afdb_gene_map)

# ------------------------------------------------------------
# -------- APPEND CLASS / SUPERCLASS ANNOTATION --------------
# ------------------------------------------------------------

gene_class_map = (
    gene_mapping[
        ["mapped_gene_name", "mapped_class", "mapped_superclass"]
    ]
    .drop_duplicates(subset="mapped_gene_name")
    .set_index("mapped_gene_name")
)

pdb_ann = pdb_ann.join(gene_class_map, on="gene_name")
afdb_ann = afdb_ann.join(gene_class_map, on="gene_name")

# ------------------------------------------------------------
# -------- REMOVE ENTRIES WITH MISSING GENE NAME -------------
# ------------------------------------------------------------

mask_missing_pdb = pdb_ann["gene_name"] == "Missing entry"
for _, row in pdb_ann[mask_missing_pdb].iterrows():
    log_removal(
        pdb=row["pdb_id"],
        uniprot_id=row["uniprot_id"] if "uniprot_id" in pdb_ann.columns else None,
        database="pdb",
        reason="missing_gene_name",
    )
pdb_ann = pdb_ann[~mask_missing_pdb]

mask_missing_afdb = afdb_ann["gene_name"] == "Missing entry"
for _, row in afdb_ann[mask_missing_afdb].iterrows():
    log_removal(
        pdb=f"AF-{row['uniprot_id']}-model",
        uniprot_id=row["uniprot_id"],
        database="afdb",
        reason="missing_gene_name",
    )
afdb_ann = afdb_ann[~mask_missing_afdb]

# ------------------------------------------------------------
# -------- REMOVE COLUMNS (GLOBAL / AFDB-SPECIFIC) -----------
# ------------------------------------------------------------

for df in (pdb_ann, afdb_ann):
    if "uniprot_length" in df.columns:
        df.drop(columns="uniprot_length", inplace=True)

if "chain" in afdb_ann.columns:
    afdb_ann.drop(columns="chain", inplace=True)

# ------------------------------------------------------------
# ------------------ TEMPORARY FILTERS -----------------------
# ------------------------------------------------------------

if "target_length" in afdb_ann.columns:
    mask_len = afdb_ann["target_length"] < 140
    for _, row in afdb_ann[mask_len].iterrows():
        log_removal(
            pdb=f"AF-{row['uniprot_id']}-model",
            uniprot_id=row["uniprot_id"],
            database="afdb",
            reason="target_length_lt_140",
        )
    afdb_ann = afdb_ann[~mask_len]

with open(AFDB_UNIPROT_REMOVE_FILE) as fh:
    uniprot_to_remove = {
        line.strip()
        for line in fh
        if line.strip() and not line.startswith("#")
    }

if "uniprot_id" in afdb_ann.columns:
    mask_rm = afdb_ann["uniprot_id"].isin(uniprot_to_remove)
    for _, row in afdb_ann[mask_rm].iterrows():
        log_removal(
            pdb=f"AF-{row['uniprot_id']}-model",
            uniprot_id=row["uniprot_id"],
            database="afdb",
            reason="explicit_uniprot_removal",
        )
    afdb_ann = afdb_ann[~mask_rm]

# ------------------------------------------------------------
# -------- REMOVE AUTHOR-ASSIGNED COLUMNS (PDB ONLY) ---------
# ------------------------------------------------------------

PDB_AUTHOR_COLUMNS = [
    "possibly_chimeric",
    "grouped_gene_assigned",
    "gene_name_assigned",
    "organism_assigned",
    "comments",
]

pdb_ann = pdb_ann.drop(
    columns=[c for c in PDB_AUTHOR_COLUMNS if c in pdb_ann.columns]
)

# ------------------------------------------------------------
# -------- PRUNE ASSEMBLIES / MODELS FILES -------------------
# ------------------------------------------------------------

pdb_assemblies["_pdb_base"] = (
    pdb_assemblies["pdb"]
    .str.split("-", n=1)
    .str[0]
    .str.upper()
)

valid_pdb_ids = set(pdb_ann["pdb_id"].str.upper())

pdb_assemblies = pdb_assemblies[
    pdb_assemblies["_pdb_base"].isin(valid_pdb_ids)
].drop(columns="_pdb_base")


valid_afdb_ids = set(afdb_ann["pdb_id"])
afdb_models = afdb_models[
    afdb_models["pdb"].isin(valid_afdb_ids)
]

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------
pdb_ann.to_csv(PDB_ANN_OUT, index=False)
afdb_ann.to_csv(AFDB_ANN_OUT, index=False)

pdb_assemblies.to_csv(PDB_ASSEMBLIES_OUT, index=False)
afdb_models.to_csv(AFDB_MODELS_OUT, index=False)

pd.DataFrame(removed_log).to_csv(REMOVED_LOG_FILE, index=False)
