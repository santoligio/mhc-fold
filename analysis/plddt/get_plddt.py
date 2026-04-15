#!/usr/bin/env python3
# ------------------------------------------------------------
# AFDB trimmed MHC models:
#  - Average pLDDT (from B-factor, residue-based)
#  - Extract amino acid sequence (1-letter)
#  - Write AFDB CSV
#  - Append avg_pLDDT + seq to mhc_annotations.csv
#   (with AF-ID normalization)
# ------------------------------------------------------------

import os
import csv
import re
from tqdm import tqdm
from Bio.PDB import PDBParser, is_aa
from Bio.SeqUtils import seq1

# === PATHS ===
INPUT_DIR  = "/mnt/4TB/giovanna/foldseek/version_02/filter/step4/afdb/trimmed_mhc"
OUTPUT_DIR = "/mnt/4TB/giovanna/foldseek/version_02/analysis/plddt"

AFDB_CSV = os.path.join(OUTPUT_DIR, "afdb_trimmed_plddt.csv")
MHC_CSV  = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/afdb/mhc_annotations.csv"
MHC_OUT  = os.path.join(OUTPUT_DIR, "mhc_annotations_with_plddt.csv")

os.makedirs(OUTPUT_DIR, exist_ok=True)

parser = PDBParser(QUIET=True)

# ------------------------------------------------------------
# Helper: normalize AF identifiers
# ------------------------------------------------------------
def normalize_af_id(name):
    """
    Extract canonical AF-<UNIPROT> identifier
    Examples:
      AF-G3SUA4_trimmed_mhc      -> AF-G3SUA4
      AF-I7B294-F1-model_v6      -> AF-I7B294
    """
    m = re.search(r"(AF-[A-Z0-9]+)", name)
    return m.group(1) if m else None


results = []
plddt_map = {}   # AF-ID -> (avg_pLDDT, seq)

pdb_files = sorted(
    f for f in os.listdir(INPUT_DIR)
    if f.endswith("_trimmed_mhc.pdb")
)

# ------------------------------------------------------------
# Process AFDB PDB files
# ------------------------------------------------------------
for fname in tqdm(pdb_files, desc="Processing AFDB models"):
    pdb_path = os.path.join(INPUT_DIR, fname)
    model_id = fname.replace(".pdb", "")

    structure = parser.get_structure(model_id, pdb_path)

    residue_plddt = []
    sequence = []

    for model in structure:
        for chain in model:
            for res in chain:
                if not is_aa(res, standard=True):
                    continue

                bvals = [atom.get_bfactor() for atom in res.get_atoms()]
                if not bvals:
                    continue

                residue_plddt.append(sum(bvals) / len(bvals))

                try:
                    sequence.append(seq1(res.get_resname()))
                except KeyError:
                    sequence.append("X")

    if not residue_plddt:
        continue

    avg_plddt = round(sum(residue_plddt) / len(residue_plddt), 2)
    seq = "".join(sequence)

    results.append([model_id, avg_plddt, seq])

    af_id = normalize_af_id(model_id)
    if af_id:
        plddt_map[af_id] = (avg_plddt, seq)

# ------------------------------------------------------------
# Write AFDB CSV (unchanged behavior)
# ------------------------------------------------------------
with open(AFDB_CSV, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["AF_model", "avg_pLDDT", "seq"])
    writer.writerows(results)

print(f"Wrote {len(results)} AFDB entries to:")
print(AFDB_CSV)

# ------------------------------------------------------------
# Merge into mhc_annotations.csv
# ------------------------------------------------------------
with open(MHC_CSV) as fin, open(MHC_OUT, "w", newline="") as fout:
    reader = csv.DictReader(fin)
    fieldnames = reader.fieldnames + ["avg_pLDDT", "seq"]

    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        # clean whitespace just in case
        row = {k: v.strip() if isinstance(v, str) else v for k, v in row.items()}

        pdb_id = row["pdb_id"]
        af_id = normalize_af_id(pdb_id)

        if af_id and af_id in plddt_map:
            row["avg_pLDDT"], row["seq"] = plddt_map[af_id]
        else:
            row["avg_pLDDT"] = ""
            row["seq"] = ""

        writer.writerow(row)

print("Updated MHC annotations written to:")
print(MHC_OUT)
