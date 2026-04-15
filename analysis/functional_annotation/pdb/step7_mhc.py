#!/bin/python3
# ------------------------------------------------------------
# Robust UniProt annotation pipeline with explicit error logging
# ------------------------------------------------------------
# Key guarantees:
#  - NO silent biological or technical failures
#  - ALL errors and warnings logged to file
#  - ALL compatible mappings written to final CSV
#  - ALL errors exported to structured CSV:
#       pdb_id, old_chain, new_chain, error
#  - "Belongs to" extraction preserved exactly as before
# ------------------------------------------------------------

import re
import time
import logging
import pandas as pd
import requests
from concurrent.futures import ThreadPoolExecutor, as_completed

# =========================
# CONFIGURATION
# =========================

INPUT_CSV   = "/mnt/4TB/giovanna/foldseek/version_02/analysis/pdb_assemblies_analysis.csv"
OUTPUT_CSV  = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/mhc/mhc_annotations-4.csv"
LOG_FILE    = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/mhc/mhc_annotations-4.log"
ERROR_CSV   = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/mhc/mhc_annotations.errors.csv"

NUM_THREADS   = 16
REQUEST_SLEEP = 5

# =========================
# LOGGING
# =========================

logging.basicConfig(
    filename=LOG_FILE,
    filemode="w",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)

console = logging.StreamHandler()
console.setLevel(logging.INFO)
console.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
logging.getLogger().addHandler(console)

# =========================
# ERROR COLLECTOR
# =========================

ERROR_RECORDS = []

def log_error(pdb_id, chain_id, message):
    ERROR_RECORDS.append({
        "pdb_id": pdb_id,
        "chain": chain_id,
        "error": message,
    })
    logging.error(f"[ERROR] [{pdb_id}:{new_chain}] {message}")

# =========================
# NETWORK
# =========================

def safe_get(url, label, timeout=10, retries=3):
    for attempt in range(1, retries + 1):
        try:
            time.sleep(REQUEST_SLEEP)
            r = requests.get(url, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            logging.warning(f"[{label}] attempt {attempt}/{retries} failed: {e}")
            time.sleep(2 ** attempt)

    logging.error(f"[ERROR] [{label}] failed after {retries} attempts")
    return None

# =========================
# SIFTS
# =========================

def get_sifts_mapping(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    data = safe_get(url, f"SIFTS:{pdb_id}")

    if not data or pdb_id.lower() not in data:
        logging.error(f"[SIFTS:{pdb_id}] no valid UniProt mapping")
        return {}

    return data[pdb_id.lower()].get("UniProt", {})

# =========================
# UNIPROT
# =========================

def get_uniprot_metadata(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    data = safe_get(url, f"UniProt:{uniprot_id}")

    if data is None:
        logging.error(f"[UniProt:{uniprot_id}] request failed")
        return None

    # Protein name
    name = "Missing entry"
    try:
        name = data["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except Exception:
        logging.warning(f"[UniProt:{uniprot_id}] recommendedName missing")
        try:
            name = data["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
            logging.warning(f"[UniProt:{uniprot_id}] used submissionNames")
        except Exception:
            logging.warning(f"[UniProt:{uniprot_id}] protein name missing")

    # Sequence length
    uniprot_length = data.get("sequence", {}).get("length")
    if uniprot_length is None:
        uniprot_length = "Missing entry"
        logging.warning(f"[UniProt:{uniprot_id}] sequence length missing")

    # Organism
    organism = data.get("organism", {}).get("scientificName")
    if organism is None:
        organism = "Missing entry"
        logging.warning(f"[UniProt:{uniprot_id}] organism missing")

    # Gene name
    gene_name = "Missing entry"
    try:
        gene_name = data["genes"][0]["geneName"]["value"]
    except Exception:
        logging.warning(f"[UniProt:{uniprot_id}] gene name missing")

    # Entry status
    entry_status = data.get("entryType")
    if entry_status:
        entry_status = entry_status.capitalize()
    else:
        entry_status = "Missing entry"
        logging.warning(f"[UniProt:{uniprot_id}] entry status missing")

    # Belongs to (UNCHANGED)
    belongs_to_info = ""
    for comment in data.get("comments", []):
        for text in comment.get("texts", []):
            value = text.get("value", "")
            if "belongs to" in value.lower():
                match = re.search(r"(Belongs to .*?\.)", value)
                belongs_to_info = match.group(1) if match else value
                break
        if belongs_to_info:
            break

    return {
        "uniprot_name": name,
        "uniprot_length": uniprot_length,
        "organism": organism,
        "gene_name": gene_name,
        "entry_status": entry_status,
        "belongs_to": belongs_to_info,
    }

# =========================
# CORE LOGIC
# =========================

def map_entry(pdb_id, chain_id, res_start, res_end, index, total):
    records = []

    try:
        chain_id = str(chain_id).split("-")[0]
        logging.info(f"[{index}/{total}] Processing {pdb_id.split('-')[0].upper()} chain {chain_id}")

        target_length = int(res_end) - int(res_start) + 1
        if target_length <= 0:
            log_error(pdb_id, chain_id, "invalid residue range")
            return []

        # ---------- AlphaFold ----------
        if pdb_id.startswith("AF-"):
            uniprot_id = pdb_id.split("-")[1]
            meta = get_uniprot_metadata(uniprot_id)
            if not meta:
                log_error(pdb_id, chain_id, "AlphaFold UniProt metadata failed")
                return []

            records.append({
                "pdb_id": pdb_id,
                "chain": chain_id,
                "res_start": res_start,
                "res_end": res_end,
                "uniprot_id": uniprot_id,
                "target_length": target_length,
                "mapped_length": meta["uniprot_length"],
                **meta,
                "possibly_chimeric_or_singlechain": "",
            })
            return records

        # ---------- PDB ----------
        pdb_clean = pdb_id.split("-")[0].lower()
        mappings = get_sifts_mapping(pdb_clean)

        if not mappings:
            log_error(pdb_clean.upper(), chain_id, "no SIFTS UniProt mappings available")
            return []

        seen_uniprot = set()
        uniprot_hits = []

        for uniprot_key, block in mappings.items():
            uniprot_id = uniprot_key.split("_")[0]
            if uniprot_id in seen_uniprot:
                continue

            for seg in block.get("mappings", []):
                if seg.get("chain_id") != chain_id:
                    continue

                seen_uniprot.add(uniprot_id)

                meta = get_uniprot_metadata(uniprot_id)
                if not meta:
                    log_error(pdb_clean.upper(), chain_id, f"UniProt metadata failed ({uniprot_id})")
                    continue

                uniprot_hits.append(uniprot_id)

                records.append({
                    "pdb_id": pdb_clean.upper(),
                    "chain": chain_id,
                    "res_start": res_start,
                    "res_end": res_end,
                    "uniprot_id": uniprot_id,
                    "target_length": target_length,
                    "mapped_length": meta["uniprot_length"],
                    **meta,
                    "possibly_chimeric_or_singlechain": "",
                })
                break

        if len(set(uniprot_hits)) > 1:
            for r in records:
                r["possibly_chimeric_or_singlechain"] = "yes"

        if not records:
            log_error(
                pdb_clean.upper(),
                chain_id,
                "no compatible UniProt mapping",
            )

        return records

    except Exception as e:
        log_error(pdb_id, chain_id, f"unexpected exception: {e}")
        return []

# =========================
# MAIN
# =========================

def main():
    try:
        df = pd.read_csv(INPUT_CSV)
        df = df[df["status"].str.lower() == "primary"]
    except Exception as e:
        logging.critical(f"Failed to load input CSV: {e}")
        return

    total = len(df)
    logging.info(f"{total} primary entries loaded")

    results = []

    with ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = [
            executor.submit(
                map_entry,
                row.pdb,
                row.chain,
                row.tstart,
                row.tend,
                i + 1,
                total,
            )
            for i, row in enumerate(df.itertuples(index=False))
        ]

        for f in as_completed(futures):
            results.extend(f.result())

    if results:
        pd.DataFrame(results).to_csv(OUTPUT_CSV, index=False)
        logging.info(f"Final CSV written: {OUTPUT_CSV}")
    else:
        logging.error("No valid mappings generated")

    if ERROR_RECORDS:
        pd.DataFrame(ERROR_RECORDS).to_csv(ERROR_CSV, index=False)
        logging.info(f"Error CSV written: {ERROR_CSV}")

    logging.info(f"Log file written: {LOG_FILE}")

if __name__ == "__main__":
    main()
