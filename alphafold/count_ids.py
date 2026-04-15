#!/usr/bin/env python3

import sys

def extract_af_uniprot(col2):
    """
    Extract UniProt accession from AlphaFold ID.
    Example:
    AF-Q95387-F1-model_v6 -> Q95387
    """
    if col2.startswith("AF-"):
        parts = col2.split("-")
        if len(parts) >= 2:
            return parts[1]
    return None


def count_unique_af_ids(filepath):
    unique_ids = set()

    with open(filepath, "r") as f:
        for line in f:
            if not line.strip():
                continue

            cols = line.strip().split()
            if len(cols) < 2:
                continue

            col2 = cols[1]
            uniprot_id = extract_af_uniprot(col2)

            if uniprot_id:
                unique_ids.add(uniprot_id)

    return unique_ids


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_unique_af_ids.py input_file.tsv")
        sys.exit(1)

    filepath = sys.argv[1]
    unique_ids = count_unique_af_ids(filepath)

    print(f"Number of unique AlphaFold UniProt IDs in column 2: {len(unique_ids)}")
    print("\nUnique IDs:")
    for uid in sorted(unique_ids):
        print(uid)
