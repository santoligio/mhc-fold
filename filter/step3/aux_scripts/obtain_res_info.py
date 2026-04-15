#!/usr/bin/env python3

import os
import csv
from collections import defaultdict

CIF_FOLDER = "/mnt/4TB/giovanna/foldseek/version_02/filter/step2/pdb/1_assemblies"
OUTPUT_CSV = "unique_residues.csv"


def extract_resnames_from_atoms(cif_path):
    """
    Extract all RESNAMEs actually used in the structure
    from _atom_site.label_comp_id
    """
    resnames = set()
    in_loop = False
    headers = []
    comp_idx = None

    with open(cif_path) as f:
        for line in f:
            line = line.strip()

            if line.startswith("loop_"):
                in_loop = True
                headers = []
                comp_idx = None
                continue

            if in_loop and line.startswith("_atom_site."):
                headers.append(line)
                if line == "_atom_site.label_comp_id":
                    comp_idx = len(headers) - 1
                continue

            if in_loop and headers and not line.startswith("_"):
                if comp_idx is None:
                    continue

                fields = line.split()
                if len(fields) > comp_idx:
                    resnames.add(fields[comp_idx])
                continue

            if in_loop and not line:
                in_loop = False

    return resnames


def extract_chem_comp_table(cif_path):
    """
    Extract _chem_comp.id → (type, name)
    """
    chem_comp = {}
    in_loop = False
    headers = []

    idx_id = idx_type = idx_name = None

    with open(cif_path) as f:
        for line in f:
            line = line.strip()

            if line.startswith("loop_"):
                in_loop = True
                headers = []
                continue

            if in_loop and line.startswith("_chem_comp."):
                headers.append(line)
                if line == "_chem_comp.id":
                    idx_id = len(headers) - 1
                elif line == "_chem_comp.type":
                    idx_type = len(headers) - 1
                elif line == "_chem_comp.name":
                    idx_name = len(headers) - 1
                continue

            if in_loop and headers and not line.startswith("_"):
                if None in (idx_id, idx_type, idx_name):
                    continue

                fields = line.split(None, len(headers) - 1)
                if len(fields) <= max(idx_id, idx_type, idx_name):
                    continue

                cid = fields[idx_id]
                ctype = fields[idx_type].strip('"')
                cname = fields[idx_name].strip('"')

                chem_comp[cid] = (ctype, cname)
                continue

            if in_loop and not line:
                in_loop = False

    return chem_comp


def collect_residue_information(cif_folder, output_csv):
    results = {}

    for fname in os.listdir(cif_folder):
        if not fname.endswith(".cif"):
            continue

        pdb_id = fname.split("-")[0].lower()
        cif_path = os.path.join(cif_folder, fname)

        resnames = extract_resnames_from_atoms(cif_path)
        chem_comp = extract_chem_comp_table(cif_path)

        if not chem_comp:
            print(f"[WARNING] {pdb_id}: no _chem_comp information found")

        for res in resnames:
            if res not in results:
                if res in chem_comp:
                    ctype, cname = chem_comp[res]
                    results[res] = (ctype, cname, "")
                else:
                    results[res] = ("NA", "NA", pdb_id)

    with open(output_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["resname", "chem_comp_type", "chem_comp_name", "missing_in_pdb"])

        for resname, (ctype, cname, missing) in sorted(results.items()):
            writer.writerow([resname, ctype, cname, missing])


if __name__ == "__main__":
    collect_residue_information(CIF_FOLDER, OUTPUT_CSV)
    print(f"[DONE] Residue summary written to {OUTPUT_CSV}")
