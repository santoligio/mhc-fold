#!/usr/bin/env python3

import sys

def count_unique_column2(filepath):
    unique_ids = set()

    with open(filepath, "r") as f:
        for line in f:
            if not line.strip():
                continue  # skip empty lines

            cols = line.strip().split("\t")

            if len(cols) < 2:
                continue  # skip malformed lines

            col2 = cols[1]

            # Extract ID before first "-"
            pdb_id = col2.split("-")[0].lower()

            unique_ids.add(pdb_id)

    return unique_ids


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python count_unique_ids.py input_file.tsv")
        sys.exit(1)

    filepath = sys.argv[1]
    unique_ids = count_unique_column2(filepath)

    print(f"Number of unique IDs in column 2: {len(unique_ids)}")
    print("\nUnique IDs:")
    for uid in sorted(unique_ids):
        print(uid)
