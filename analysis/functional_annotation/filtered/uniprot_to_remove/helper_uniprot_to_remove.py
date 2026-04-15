#!/usr/bin/env python3

import requests
import pandas as pd

# ------------------------------------------------------------
# Inputs
# ------------------------------------------------------------
AFDB_MODELS_FILE = (
    "/mnt/4TB/giovanna/foldseek/version_02/filter/step1/afdb/afdb_models.csv"
)

MHC_ANNOTATIONS_FILE = (
    "/mnt/4TB/giovanna/foldseek/version_02/analysis/"
    "functional_annotation/pdb/mhc/pdb_mhc_annotations.csv"
)

UNIPROT_REMOVE_URL = (
    "https://ftp.ebi.ac.uk/pub/contrib/UniProt/proteomes/"
    "proteins_to_remove_from_UniProtKB.txt"
)

# ------------------------------------------------------------
# Output
# ------------------------------------------------------------
OUTPUT_FILE = "uniprot_proteins_to_remove.txt"

# ------------------------------------------------------------
# Load AFDB models and extract UniProt IDs
# ------------------------------------------------------------
afdb_models = pd.read_csv(AFDB_MODELS_FILE)

if "pdb" not in afdb_models.columns:
    raise ValueError("afdb_models.csv must contain a 'pdb' column")

def extract_uniprot_from_afdb(model_id):
    """
    Example:
    AF-Q2YHQ6-F1-model_v6 -> Q2YHQ6
    """
    try:
        return model_id.split("-")[1]
    except Exception:
        return None

afdb_uniprot_ids = set(
    afdb_models["pdb"]
    .dropna()
    .astype(str)
    .apply(extract_uniprot_from_afdb)
    .dropna()
)

print(f"Extracted {len(afdb_uniprot_ids)} UniProt IDs from AFDB models")

# ------------------------------------------------------------
# Load UniProt IDs from MHC annotations
# ------------------------------------------------------------
mhc_ann = pd.read_csv(MHC_ANNOTATIONS_FILE)

if "uniprot_id" not in mhc_ann.columns:
    raise ValueError("mhc_annotations.csv must contain a 'uniprot_id' column")

mhc_uniprot_ids = set(
    mhc_ann["uniprot_id"]
    .dropna()
    .astype(str)
    .str.strip()
)

print(f"Extracted {len(mhc_uniprot_ids)} UniProt IDs from MHC annotations")

# ------------------------------------------------------------
# Union of all UniProt IDs to check
# ------------------------------------------------------------
all_uniprot_ids = afdb_uniprot_ids | mhc_uniprot_ids

print(f"Total unique UniProt IDs to check: {len(all_uniprot_ids)}")

# ------------------------------------------------------------
# Stream UniProt removal list and match
# ------------------------------------------------------------
matched_ids = set()

with requests.get(UNIPROT_REMOVE_URL, stream=True) as r:
    r.raise_for_status()

    for line in r.iter_lines(decode_unicode=True):
        if not line or line.startswith("#"):
            continue

        uniprot_id = line.strip()

        if uniprot_id in all_uniprot_ids:
            matched_ids.add(uniprot_id)

print(f"Found {len(matched_ids)} UniProt IDs marked for removal")

# ------------------------------------------------------------
# Write output
# ------------------------------------------------------------
with open(OUTPUT_FILE, "w") as out:
    for uid in sorted(matched_ids):
        out.write(f"{uid}\n")

print(f"Wrote results to {OUTPUT_FILE}")
