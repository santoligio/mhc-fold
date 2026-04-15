#!/usr/bin/env python3

import os
import string
import pandas as pd

from itertools import product
from Bio.PDB import MMCIFParser, MMCIFIO, NeighborSearch


# =========================
# Configuration
# =========================

CIF_DIR   = "/mnt/4TB/giovanna/foldseek/version_02/filter/step2/pdb/1_assemblies"
CSV_FILE = "/mnt/4TB/giovanna/foldseek/version_02/filter/step1/pdb/pdb_assemblies.csv"

BASE_OUT_DIR     = "/mnt/4TB/giovanna/foldseek/version_02/filter/step3/pdb-teste"
REMAPPED_CIF_DIR = os.path.join(BASE_OUT_DIR, "remapped_cifs")

GAPS_CSV = os.path.join(BASE_OUT_DIR, "mhc_chain_gaps.csv")

CUTOFF_DISTANCE = 10.0

os.makedirs(BASE_OUT_DIR, exist_ok=True)
os.makedirs(REMAPPED_CIF_DIR, exist_ok=True)


# =========================
# Utility: load MHC chain mapping
# =========================
def load_mhc_chain_map(csv_file):

    df = pd.read_csv(csv_file)
    df_primary = df[df["status"].str.lower() == "primary"]

    mhc_map = {}
    for _, row in df_primary.iterrows():
        pdb_id = row["pdb"].split("-")[0].upper()
        mhc_map[pdb_id] = row["chain"]

    return df, df_primary, mhc_map


# =========================
# Renumber ONLY the MHC chain
# =========================
def renumber_mhc_chain(chain):

    residues = list(chain.get_residues())
    for res in residues:
        chain.detach_child(res.id)

    new_resseq = 1
    for res in residues:
        hetflag, _, _ = res.id
        res.id = (hetflag, new_resseq, " ")
        chain.add(res)
        new_resseq += 1


# =========================
# Cutoff filter
# =========================
def find_far_chains(structure, mhc_chain_id, cutoff=CUTOFF_DISTANCE):

    model = structure[0]
    if mhc_chain_id not in model:
        return set()

    mhc_chain = model[mhc_chain_id]
    mhc_atoms = [a for a in mhc_chain.get_atoms() if a.id == "CA"]
    if not mhc_atoms:
        return set()

    ns = NeighborSearch(list(structure.get_atoms()))
    near_chain_ids = set()

    for atom in mhc_atoms:
        for nb in ns.search(atom.coord, cutoff, level="A"):
            near_chain_ids.add(nb.get_parent().get_parent().id)

    far = set(chain.id for chain in model) - near_chain_ids
    far.discard(mhc_chain_id)
    return far


# =========================
# Chain remapping utilities
# =========================
def rename_chain_tmp(structure):

    rename_needed = {}
    for model in structure:
        for chain in model:
            old_id = chain.id
            tmp_id = f"_TMP_{old_id}"
            chain.id = tmp_id
            rename_needed[tmp_id] = old_id
    return rename_needed


def chain_id_generator():

    letters = string.ascii_uppercase
    n = 1
    while True:
        for combo in product(letters, repeat=n):
            yield "".join(combo)
        n += 1


def remap_chain_ids(structure, rename_needed, mhc_chain_id):

    original_ids = sorted(rename_needed.values())
    binder_ids = [cid for cid in original_ids if cid != mhc_chain_id]

    gen = chain_id_generator()
    target_ids = {mhc_chain_id: next(gen)}

    for cid in binder_ids:
        target_ids[cid] = next(gen)

    for model in structure:
        for chain in model:
            if chain.id in rename_needed:
                chain.id = target_ids[rename_needed[chain.id]]

    return target_ids


# =========================
# GAP CALCULATION (NEW)
# =========================
def calculate_mhc_gaps(chain, tstart, tend):
    """
    tstart/tend are in RENumbered coordinates.
    Gaps are counted in AUTH coordinates.
    """

    # Collect residues in original auth order
    auth_resseqs = []
    for res in chain.get_residues():
        hetflag, auth_resseq, _ = res.id
        if hetflag.strip():
            continue
        auth_resseqs.append(auth_resseq)

    if not auth_resseqs:
        return 0

    # Map renumbered index → auth resseq
    renum_to_auth = {
        idx + 1: auth
        for idx, auth in enumerate(auth_resseqs)
    }

    # Extract auth residues corresponding to tstart–tend
    region_auth = [
        renum_to_auth[i]
        for i in range(tstart, tend + 1)
        if i in renum_to_auth
    ]

    if len(region_auth) < 2:
        return 0

    # Count gaps in auth numbering
    total_gaps = 0
    for prev, curr in zip(region_auth, region_auth[1:]):
        gap = curr - prev - 1
        if gap > 0:
            total_gaps += gap

    return total_gaps


# =========================
# Main per-file processing
# =========================
def process_cif(cif_path, mhc_chain_id, tstart, tend, out_dir):

    pdb_id = os.path.basename(cif_path).split("-")[0].upper()

    parser = MMCIFParser(QUIET=True, auth_chains=True)
    structure = parser.get_structure(pdb_id, cif_path)
    model = structure[0]

    if mhc_chain_id not in model:
        print(f"[SKIP] {pdb_id}: MHC chain {mhc_chain_id} not found")
        return None, None

    mhc_chain = model[mhc_chain_id]

    # --- GAP CALCULATION BEFORE RENNUMBERING ---
    total_gaps = calculate_mhc_gaps(mhc_chain, tstart, tend)

    # --- RENNUMBER MHC ---
    renumber_mhc_chain(mhc_chain)

    far_chains = find_far_chains(structure, mhc_chain_id, CUTOFF_DISTANCE)
    for cid in far_chains:
        if cid in model:
            model.detach_child(cid)

    if far_chains:
        print(f"[FILTER] {pdb_id}: removed far chains {sorted(far_chains)}")

    rename_needed = rename_chain_tmp(structure)
    target_ids = remap_chain_ids(structure, rename_needed, mhc_chain_id)

    out_cif = os.path.join(out_dir, f"{pdb_id.lower()}_remapped.cif")
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(out_cif)

    print(f"[OK] {pdb_id} → {out_cif}")

    return target_ids, total_gaps


# =========================
# Entry point
# =========================
def main():

    df_csv, df_primary, mhc_map = load_mhc_chain_map(CSV_FILE)

    mapping_rows = []
    updated_rows = []
    gap_rows = []

    for _, row in df_primary.iterrows():

        pdb_id = row["pdb"].split("-")[0].upper()
        cif_name = f"{pdb_id.lower()}-assembly1.cif"
        cif_path = os.path.join(CIF_DIR, cif_name)

        if not os.path.exists(cif_path):
            continue

        chain_map, gaps = process_cif(
            cif_path,
            row["chain"],
            int(row["tstart"]),
            int(row["tend"]),
            REMAPPED_CIF_DIR
        )

        if chain_map is None:
            continue

        gap_rows.append({
            "pdb_id": pdb_id,
            "total_gaps": gaps
        })

        for old_chain, new_chain in chain_map.items():
            mapping_rows.append({
                "pdb": pdb_id,
                "old_chain": old_chain,
                "new_chain": new_chain
            })

    # --- write chain mapping CSV ---
    if mapping_rows:
        mapping_df = pd.DataFrame(mapping_rows)
        mapping_csv = os.path.join(BASE_OUT_DIR, "chain_id_mapping.csv")
        mapping_df.to_csv(mapping_csv, index=False)
        print(f"[MAPPING CSV WRITTEN] {mapping_csv}")

    # --- write gaps CSV ---
    if gap_rows:
        gaps_df = pd.DataFrame(gap_rows)
        gaps_df.to_csv(GAPS_CSV, index=False)
        print(f"[GAPS CSV WRITTEN] {GAPS_CSV}")

    # --- write remapped pdb_assemblies CSV ---
    for _, row in df_csv.iterrows():

        pdb_id = row["pdb"].split("-")[0].upper()
        old_chain = row["chain"]

        matches = [
            m for m in mapping_rows
            if m["pdb"] == pdb_id and m["old_chain"] == old_chain
        ]

        if not matches:
            continue

        row = row.copy()
        row["chain"] = matches[0]["new_chain"]
        updated_rows.append(row)

    if updated_rows:
        out_df = pd.DataFrame(updated_rows)
        out_csv = os.path.join(BASE_OUT_DIR, "pdb_assemblies_remapped.csv")
        out_df.to_csv(out_csv, index=False)
        print(f"[ASSEMBLIES CSV WRITTEN] {out_csv}")


if __name__ == "__main__":
    main()
