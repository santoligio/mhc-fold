#!/usr/bin/env python3

import csv
from Bio.PDB import Residue, is_aa


INPUT_CSV  = "unique_residues.csv"
OUTPUT_CSV = "unique_residues.csv"   # overwrite same file


def biopython_recognizes_resname(resname):
    """
    Create a dummy Residue and test with Biopython is_aa
    """
    dummy_residue = Residue.Residue(
        (" ", 1, " "),
        resname,
        ""
    )
    return is_aa(dummy_residue, standard=False)


def annotate_csv_with_biopython(input_csv, output_csv):
    rows = []

    with open(input_csv, newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames + ["biopython_included"]

        for row in reader:
            resname = row["resname"]
            row["biopython_included"] = (
                "yes" if biopython_recognizes_resname(resname) else "no"
            )
            rows.append(row)

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    annotate_csv_with_biopython(INPUT_CSV, OUTPUT_CSV)
    print(f"[DONE] Biopython annotation added to {OUTPUT_CSV}")
