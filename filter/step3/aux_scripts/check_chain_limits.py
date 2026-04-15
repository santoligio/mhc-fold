#!/usr/bin/env python3

import os
from Bio.PDB import MMCIFParser


# =========================
# Configuration
# =========================

CIF_DIR = "/mnt/4TB/giovanna/foldseek/version_02/filter/step3/remapped_cifs"


# =========================
# Core logic
# =========================

def count_chains(cif_path):
    """
    Parse CIF and count chains in model 0.

    Returns:
        int: number of chains
    """
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("struct", cif_path)

    model = structure[0]
    return len(list(model.get_chains()))


def scan_cif_directory(cif_dir):
    """
    Scan directory and report CIFs with >26 chains.
    """
    too_many = []
    total = 0

    for fname in sorted(os.listdir(cif_dir)):
        if not fname.endswith(".cif"):
            continue

        total += 1
        cif_path = os.path.join(cif_dir, fname)

        try:
            n_chains = count_chains(cif_path)
        except Exception as e:
            print(f"[ERROR] {fname}: {e}")
            continue

        if n_chains > 26:
            pdb_id = fname.split("-")[0].upper()
            too_many.append((pdb_id, fname, n_chains))
            print(f"[PDB-UNSAFE] {fname} → {n_chains} chains")

    print("\n===== SUMMARY =====")
    print(f"Total CIFs scanned: {total}")
    print(f"CIFs with >26 chains: {len(too_many)}")

    return too_many


# =========================
# Entry point
# =========================

if __name__ == "__main__":
    scan_cif_directory(CIF_DIR)
