#!/usr/bin/env python3
# ------------------------------------------------------------
# Fetch PDB resolution from RCSB
# ------------------------------------------------------------
# Inputs:
#   - pdb_assemblies_analysis.csv (must contain pdb_id)
#   - mhc_annotations.csv
#
# Outputs:
#   1) pdb_resolution.csv  -> pdb_id,resolution (Brazilian format)
#   2) mhc_annotations_with_resolution.csv
# ------------------------------------------------------------

import csv
import time
import requests
import os
from tqdm import tqdm

# === PATHS ===
BASE_DIR = "/mnt/4TB/giovanna/foldseek/version_02/analysis/resolution"

ASSEMBLIES_CSV = "/mnt/4TB/giovanna/foldseek/version_02/analysis/pdb_assemblies_analysis.csv"
MHC_CSV        = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/pdb/mhc/pdb_mhc_annotations_edited.csv"

RESO_CSV_OUT   = os.path.join(BASE_DIR, "pdb_resolution.csv")
MHC_OUT        = os.path.join(BASE_DIR, "mhc_annotations_with_resolution.csv")

# === RCSB API ===
RCSB_ENTRY_API = "https://data.rcsb.org/rest/v1/core/entry/{}"
REQUEST_DELAY = 0.2  # seconds (≈5 requests/sec, RCSB-safe)

# ------------------------------------------------------------
# Helper: Brazilian resolution formatting
# ------------------------------------------------------------
def format_resolution_br(res):
    """
    Convert resolution to Brazilian format:
    2.01 -> '2,01'
    """
    if res == "" or res is None:
        return ""
    try:
        return f"{float(res):.2f}".replace(".", ",")
    except ValueError:
        return ""

# ------------------------------------------------------------
# Load PDB IDs
# ------------------------------------------------------------
pdb_ids = set()

with open(ASSEMBLIES_CSV) as f:
    reader = csv.DictReader(f)
    for row in reader:
        pdb = row["pdb"].split("-")[0].lower()
        if pdb:
            pdb_ids.add(pdb)

pdb_ids = sorted(pdb_ids)

print(f"Found {len(pdb_ids)} unique PDB IDs")

# ------------------------------------------------------------
# Fetch resolution data (with progress bar)
# ------------------------------------------------------------
resolution_map = {}  # pdb_id -> resolution (float or "")

for pdb_id in tqdm(pdb_ids, desc="Fetching RCSB resolutions"):
    url = RCSB_ENTRY_API.format(pdb_id)

    try:
        r = requests.get(url, timeout=10)

        if r.status_code != 200:
            resolution_map[pdb_id] = ""
            continue

        data = r.json()

        res = data.get("rcsb_entry_info", {}) \
                  .get("resolution_combined", [])

        resolution_map[pdb_id] = res[0] if res else ""

    except Exception:
        resolution_map[pdb_id] = ""

    time.sleep(REQUEST_DELAY)

print(f"Resolution retrieved for {len(resolution_map)} entries")

# ------------------------------------------------------------
# Write resolution-only CSV (Brazilian format)
# ------------------------------------------------------------
with open(RESO_CSV_OUT, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["pdb_id", "resolution"])
    for pdb_id in pdb_ids:
        res = resolution_map.get(pdb_id, "")
        writer.writerow([pdb_id.upper(), format_resolution_br(res)])

print("Wrote resolution table to:")
print(RESO_CSV_OUT)

# ------------------------------------------------------------
# Append resolution to mhc_annotations.csv (Brazilian format)
# ------------------------------------------------------------
with open(MHC_CSV) as fin, open(MHC_OUT, "w", newline="") as fout:
    reader = csv.DictReader(fin)
    fieldnames = reader.fieldnames + ["resolution"]

    writer = csv.DictWriter(fout, fieldnames=fieldnames)
    writer.writeheader()

    for row in reader:
        pdb_id = row["pdb_id"].upper()
        res = resolution_map.get(pdb_id, "")
        row["resolution"] = format_resolution_br(res)
        writer.writerow(row)

print("Updated MHC annotations written to:")
print(MHC_OUT)
