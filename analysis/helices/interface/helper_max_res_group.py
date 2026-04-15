#!/usr/bin/env python3

import re
from collections import defaultdict

# ------------------------------------------------------------
# Residue grouping
# ------------------------------------------------------------

AA_GROUPS = {
    "Sulfur polymerization": ['CYS'],
    "Small": ['GLY', 'SER', 'THR', 'ALA', 'PRO'],
    "Acid and amide": ['ASP', 'GLU', 'ASN', 'GLN'],
    "Basic": ['ARG', 'HIS', 'LYS'],
    "Hydrophobic": ['LEU', 'VAL', 'MET', 'ILE'],
    "Aromatic": ['TYR', 'PHE', 'TRP'],
}

# Reverse lookup: residue -> group
AA_TO_GROUP = {
    residue: group
    for group, residues in AA_GROUPS.items()
    for residue in residues
}

# ------------------------------------------------------------
# Parser
# ------------------------------------------------------------

def parse_file(filename):
    data = defaultdict(lambda: defaultdict(float))
    current_class = None

    with open(filename) as f:
        for line in f:
            line = line.strip()

            if line.startswith("Group:"):
                current_class = line.split(":")[1].strip()

            elif "(" in line and "%" in line:
                match = re.search(r"\((\w{3})\)\s*:\s*([\d.]+)%", line)
                if match:
                    residue = match.group(1)
                    percent = float(match.group(2))

                    if residue in AA_TO_GROUP:
                        group_type = AA_TO_GROUP[residue]
                        data[current_class][group_type] += percent

    return data


# ------------------------------------------------------------
# Ranking
# ------------------------------------------------------------

def print_top_three(data):
    print("\n=== Top 3 classes per residue group ===\n")

    for group_type in AA_GROUPS.keys():

        # Build ranking list
        ranking = []
        for class_name, values in data.items():
            value = values.get(group_type, 0.0)
            ranking.append((class_name, value))

        # Sort descending
        ranking.sort(key=lambda x: x[1], reverse=True)

        print(f"\n{group_type}")
        print("-" * len(group_type))

        for i, (class_name, value) in enumerate(ranking[:3], start=1):
            print(f"{i}. {class_name:10s} {value:6.2f}%")

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

if __name__ == "__main__":
    filename = "helix_statistics_pdb_class.txt"  # change to your file
    data = parse_file(filename)
    print_top_three(data)
