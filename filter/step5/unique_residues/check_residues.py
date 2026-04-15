#!/usr/bin/env python3
import os
import pandas as pd
from Bio.PDB import PDBParser, is_aa

# === CONFIGURATION ===
PDB_DIR = "/mnt/4TB/giovanna/foldseek/version_02/filter/step4/pdb/remapped_cifs"
RES_CSV = "/mnt/4TB/giovanna/foldseek/version_02/filter/step3/aux_scripts/unique_residues.csv"
CHAIN_CSV = "/mnt/4TB/giovanna/foldseek/version_02/filter/step5/pdb/mhc_contacts.csv"
OUTPUT_CSV = "/mnt/4TB/giovanna/foldseek/version_02/filter/step5/unique_residues/unique_residues2.csv"

# === LOAD AND NORMALIZE CSV FILES ===
res_df = pd.read_csv(RES_CSV)
res_df['resname'] = res_df['resname'].astype(str).str.strip().str.upper()

chain_df = pd.read_csv(CHAIN_CSV)
chain_df['pdb_id'] = chain_df['pdb_id'].astype(str).str.strip()
chain_df['Chain'] = chain_df['Chain'].astype(str).str.strip().str.upper()
chain_df['type'] = chain_df['type'].astype(str).str.strip().str.lower()

# Prepare output storage
output_rows = []
seen = set()  # track (resname, chain_type) to avoid repeats

# Initialize PDB parser
parser = PDBParser(QUIET=True)

# Iterate over PDB files
for pdb_file in os.listdir(PDB_DIR):
    if not pdb_file.endswith(".pdb"):
        continue

    # Extract base pdb_id from filename (before _filtered_complex)
    pdb_id = os.path.basename(pdb_file).split('_')[0]

    pdb_path = os.path.join(PDB_DIR, pdb_file)

    try:
        structure = parser.get_structure(pdb_id, pdb_path)
    except Exception as e:
        print(f"ERROR: Could not parse {pdb_file}: {e}")
        continue

    # Filter chain types for this PDB
    pdb_chain_types = chain_df[chain_df['pdb_id'] == pdb_id]

    # Iterate over models, chains, residues
    for model in structure:
        for chain in model:
            chain_id = chain.get_id().strip().upper()  # normalize PDB chain ID

            # Determine chain type
            chain_type_row = pdb_chain_types[pdb_chain_types['Chain'] == chain_id]
            chain_type = chain_type_row['type'].values[0] if not chain_type_row.empty else "unknown"

            for residue in chain:
                # Skip standard amino acids
                if is_aa(residue, standard=True):
                    continue

                # Normalize residue name
                resname = residue.get_resname().strip().upper()

                # Use combination of residue + chain_type to avoid duplicates
                key = (resname, chain_type)
                if key in seen:
                    continue

                match = res_df[res_df['resname'] == resname]
                if not match.empty:
                    for _, row in match.iterrows():
                        row_dict = row.to_dict()
                        row_dict['pdb_id'] = pdb_id
                        row_dict['chain'] = chain_id
                        row_dict['chain_type'] = chain_type
                        output_rows.append(row_dict)
                        seen.add(key)
                else:
                    print(f"WARNING: Residue {resname} in {pdb_file} chain {chain_id} not found in reference CSV")

# Save results
if output_rows:
    out_df = pd.DataFrame(output_rows)
    out_df.to_csv(OUTPUT_CSV, index=False)
    print(f"Saved output CSV with {len(output_rows)} entries to {OUTPUT_CSV}")
else:
    print("No matching residues found.")
