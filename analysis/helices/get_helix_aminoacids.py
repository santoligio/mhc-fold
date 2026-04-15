#!/usr/bin/env python3

import os
import argparse
import tempfile
from pathlib import Path
from collections import defaultdict

import pandas as pd
from Bio.PDB import PDBParser, MMCIFIO
from Bio.PDB.DSSP import DSSP
from Bio.Data.IUPACData import protein_letters_1to3_extended


AA_FULLNAME = {
    'ALA': 'Alanine', 'ARG': 'Arginine', 'ASN': 'Asparagine',
    'ASP': 'Aspartic acid', 'CYS': 'Cysteine', 'GLN': 'Glutamine',
    'GLU': 'Glutamic acid', 'GLY': 'Glycine', 'HIS': 'Histidine',
    'ILE': 'Isoleucine', 'LEU': 'Leucine', 'LYS': 'Lysine',
    'MET': 'Methionine', 'PHE': 'Phenylalanine', 'PRO': 'Proline',
    'SER': 'Serine', 'THR': 'Threonine', 'TRP': 'Tryptophan',
    'TYR': 'Tyrosine', 'VAL': 'Valine'
}

def parse_resnum(value):
    if pd.isna(value):
        return None
    return int(value)


def write_temp_cif(structure):
    model = structure[0]
    io = MMCIFIO()
    io.set_structure(model)
    tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".cif")
    tmp_path = tmp_file.name
    tmp_file.close()
    io.save(tmp_path)
    return tmp_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb_dir", required=True)
    parser.add_argument("--info_csv", required=True)
    parser.add_argument("--dssp_exec", default="mkdssp")
    parser.add_argument("--mode", choices=["all", "class", "superclass"], default="all")
    parser.add_argument("--annotations_csv", required=False)

    args = parser.parse_args()

    pdb_dir = Path(args.pdb_dir)

    helix_info = pd.read_csv(args.info_csv)
    helix_info = helix_info[helix_info["status"] == "OK"]

    helix_cols = ['helix1_start', 'helix1_end', 'helix2_start', 'helix2_end']
    for col in helix_cols:
        helix_info[col] = pd.to_numeric(helix_info[col], errors='coerce').astype('Int64')

    for col in ['file', 'chain']:
        helix_info[col] = helix_info[col].astype(str).str.strip()

    annotation_map = {}

    if args.mode in ["class", "superclass"]:
        if args.annotations_csv is None:
            raise ValueError("annotations_csv is required when mode is class or superclass")

        annotations = pd.read_csv(args.annotations_csv)

        # Normalize PDB IDs
        if "pdb_id" in annotations.columns:
            annotations["pdb_id"] = annotations["pdb_id"].astype(str).str.lower().str[:4]

        # Normalize UniProt IDs
        if "uniprot_id" in annotations.columns:
            annotations["uniprot_id"] = annotations["uniprot_id"].astype(str).str.upper()

        for _, row in annotations.iterrows():
            if "pdb_id" in row and pd.notna(row.get("pdb_id", None)):
                annotation_map[row["pdb_id"]] = {
                    "mapped_class": row["mapped_class"],
                    "mapped_superclass": row["mapped_superclass"]
                }
            if "uniprot_id" in row and pd.notna(row.get("uniprot_id", None)):
                annotation_map[row["uniprot_id"]] = {
                    "mapped_class": row["mapped_class"],
                    "mapped_superclass": row["mapped_superclass"]
                }

    parser_pdb = PDBParser(QUIET=True)

    group_counts = defaultdict(lambda: defaultdict(int))
    group_totals = defaultdict(int)

    for pdb_file in pdb_dir.glob("*.pdb"):

        stem = pdb_file.stem

        # ---------- Unified ID logic ----------
        if stem.upper().startswith("AF-"):
            base = stem.split("_")[0]
            uniprot_id = base.split("-")[1].upper()
            id_key = uniprot_id
        else:
            id_key = stem[:4].lower()
        # --------------------------------------

        print(f"Reading {stem} ... ", end="", flush=True)

        if pdb_file.name not in helix_info["file"].values:
            print("SKIPPED (not in helix CSV)")
            continue

        matching_rows = helix_info[helix_info["file"] == pdb_file.name]

        try:
            structure = parser_pdb.get_structure(stem, pdb_file)
        except Exception as e:
            print(f"ERROR: {type(e).__name__}")
            continue

        try:
            tmp_cif_path = write_temp_cif(structure)
            model = structure[0]
            dssp = DSSP(model, tmp_cif_path, dssp=args.dssp_exec)
        except Exception as e:
            print(f"ERROR: {type(e).__name__}")
            if os.path.exists(tmp_cif_path):
                os.unlink(tmp_cif_path)
            continue

        if args.mode == "all":
            group_key = "all"
        else:
            group_key = annotation_map.get(id_key, None)
            if group_key is None:
                print("SKIPPED (no annotation match)")
                os.unlink(tmp_cif_path)
                continue
            group_key = group_key[f"mapped_{args.mode}"]

        for _, row in matching_rows.iterrows():
            chain_id = row["chain"]
            if chain_id.upper() == "NOCHAININFO":
                chain_id = list(model.child_dict.keys())[0]

            h1_start = parse_resnum(row["helix1_start"])
            h1_end   = parse_resnum(row["helix1_end"])
            h2_start = parse_resnum(row["helix2_start"])
            h2_end   = parse_resnum(row["helix2_end"])

            if None in [h1_start, h1_end, h2_start, h2_end]:
                continue

            for key in dssp.keys():
                chain, res_id = key
                hetfield, resseq, icode = res_id

                if chain != chain_id or hetfield != " ":
                    continue

                in_helix = (h1_start <= resseq <= h1_end) or (h2_start <= resseq <= h2_end)
                if not in_helix:
                    continue

                aa1 = dssp[key][1]
                aa3 = protein_letters_1to3_extended.get(aa1.upper(), aa1).upper()
                rsa = dssp[key][3]

                if rsa is None:
                    continue

                if rsa > 0.2:
                    group_counts[group_key][aa3] += 1
                    group_totals[group_key] += 1

        os.unlink(tmp_cif_path)
        print("OK")

    output_file = f"helix_statistics_{args.mode}.txt"

    with open(output_file, "w") as out:
        out.write("=== Results ===\n\n")

        for group_key in sorted(group_counts.keys()):
            total_selected = group_totals[group_key]

            if total_selected == 0:
                out.write(f"{group_key}: No residues matched criteria.\n\n")
                continue

            out.write(f"Group: {group_key}\n")
            for resname, count in sorted(group_counts[group_key].items()):
                percent = (count / total_selected) * 100
                full_name = AA_FULLNAME.get(resname, resname)
                out.write(f"{full_name:20s} ({resname}) : {percent:6.2f}%\n")

            out.write(f"Total residues considered: {total_selected}\n\n")

    print(f"\nResults written to {output_file}")


if __name__ == "__main__":
    main()
