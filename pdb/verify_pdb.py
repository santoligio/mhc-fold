# file: check_pdb_groove.py

# Replace these with your actual filenames
index_file = "index.txt"       # File with PDB codes, one per line
search_file = "dbf_pdb_aln"     # File to search for {PDB}_groove_aligned

# Read the entire search file into memory
with open(search_file, "r") as f:
    search_lines = f.readlines()

# Read PDB codes from the index file
with open(index_file, "r") as f:
    pdb_codes = [line.strip() for line in f if line.strip()]

# Check each PDB code
for pdb in pdb_codes:
    target = f"{pdb}_groove_aligned"
    found = False
    
    for line in search_lines:
        if target in line:
            found = True
            break
    
    if found:
        print(f"{pdb}: ✅ OK")
    else:
        print(f"{pdb}: ❌ Not found")
