#!/usr/bin/env python3

import os
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select, is_aa

# ============================================================
# CONFIGURATION
# ============================================================
DATABASE = "pdb"

PDB_INPUT_DIR = (
    f"/mnt/4TB/giovanna/foldseek/version_02/filter/step5/"
    f"{DATABASE}/trimmed_mhc_complexes/"
)

CONTACTS_CSV = (
    f"/mnt/4TB/giovanna/foldseek/version_02/filter/step5/"
    f"{DATABASE}/mhc_contacts.csv"
)

OUTPUT_BASE = (
    f"/mnt/4TB/giovanna/foldseek/version_02/filter/step6/"
    f"{DATABASE}"
)

OUTPUT_MHC_ONLY     = os.path.join(OUTPUT_BASE, "mhc_only")
OUTPUT_BINDERS      = os.path.join(OUTPUT_BASE, "binders")
OUTPUT_MHC_ALL      = os.path.join(OUTPUT_BASE, "mhc_all")
OUTPUT_MHC_LIGANDS  = os.path.join(OUTPUT_BASE, "mhc_ligands")

FLAGGED_CSV = os.path.join(
    OUTPUT_BASE,
    "nonstandard_residues.csv"
)

for d in [
    OUTPUT_MHC_ONLY,
    OUTPUT_BINDERS,
    OUTPUT_MHC_ALL,
    OUTPUT_MHC_LIGANDS
]:
    os.makedirs(d, exist_ok=True)

# ============================================================
# SELECTOR
# ============================================================
class ChainSetSelect(Select):
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return chain.id in self.chain_ids

# ============================================================
# UTIL
# ============================================================
def chain_has_nonstandard_residues(chain):
    for res in chain:
        if not is_aa(res, standard=True):
            return True
    return False

# ============================================================
# MAIN
# ============================================================
def main():

    contacts_df = pd.read_csv(CONTACTS_CSV)

    binders = {}
    ligands = {}

    for pdb_id, group in contacts_df.groupby("pdb_id"):
        binders[pdb_id] = list(
            group[group["type"] == "binder"]["Chain"]
        )
        ligands[pdb_id] = list(
            group[group["type"] == "ligand"]["Chain"]
        )

    parser = PDBParser(QUIET=True)
    flagged = []

    for fname in os.listdir(PDB_INPUT_DIR):

        if not fname.endswith(".pdb"):
            continue

        pdb_id = fname.replace("_filtered_complex.pdb", "")
        pdb_file = os.path.join(PDB_INPUT_DIR, fname)

        structure = parser.get_structure(pdb_id, pdb_file)
        model = structure[0]

        if "A" not in model:
            print(f"[WARNING] {pdb_id}: MHC chain A not found")
            continue

        io = PDBIO()
        io.set_structure(structure)

        # ====================================================
        # MHC ONLY
        # ====================================================
        if chain_has_nonstandard_residues(model["A"]):
            flagged.append((pdb_id, "A"))

        io.save(
            os.path.join(
                OUTPUT_MHC_ONLY,
                f"{pdb_id}_MHC_groove.pdb"
            ),
            ChainSetSelect(["A"])
        )

        # ====================================================
        # BINDERS — ONLY binder chain
        # ====================================================
        for c in binders.get(pdb_id, []):

            if c not in model:
                print(f"[WARNING] {pdb_id}: binder chain {c} not found")
                continue

            if chain_has_nonstandard_residues(model[c]):
                flagged.append((pdb_id, c))

            io.save(
                os.path.join(
                    OUTPUT_BINDERS,
                    f"{pdb_id}_binder_chain{c}.pdb"
                ),
                ChainSetSelect([c])
            )

        # ====================================================
        # MHC + LIGANDS
        # ====================================================
        ligand_chains = [
            c for c in ligands.get(pdb_id, [])
            if c in model
        ]

        mhc_all_chains = ["A"] + ligand_chains

        # ---- mhc_all: ALWAYS written
        io.save(
            os.path.join(
                OUTPUT_MHC_ALL,
                f"{pdb_id}_MHC_all.pdb"
            ),
            ChainSetSelect(mhc_all_chains)
        )

        # ---- mhc_ligands: ONLY if ligands exist
        if ligand_chains:
            for c in ligand_chains:
                if chain_has_nonstandard_residues(model[c]):
                    flagged.append((pdb_id, c))

            io.save(
                os.path.join(
                    OUTPUT_MHC_LIGANDS,
                    f"{pdb_id}_MHC_ligands.pdb"
                ),
                ChainSetSelect(mhc_all_chains)
            )

    # ====================================================
    # REPORT
    # ====================================================
    if flagged:
        pd.DataFrame(
            sorted(set(flagged)),
            columns=["pdb_id", "Chain"]
        ).to_csv(FLAGGED_CSV, index=False)

        print(f"[INFO] Non-standard residues written to {FLAGGED_CSV}")

    print("[DONE] Step6 finished")

if __name__ == "__main__":
    main()
