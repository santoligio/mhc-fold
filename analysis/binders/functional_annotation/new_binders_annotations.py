#!/bin/python3
# ------------------------------------------------------------
# Robust binder UniProt + GO fallback annotation pipeline
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

INPUT_CSV      = "/mnt/4TB/giovanna/foldseek/version_02/filter/step5/pdb/mhc_contacts.csv"
CHAIN_MAP_CSV  = "/mnt/4TB/giovanna/foldseek/version_02/filter/step3/pdb/chain_id_mapping.csv"
FILTER_CSV     = "/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/filtered/pdb_assemblies_analysis_filtered.csv"

OUTPUT_CSV = "/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/functional_annotation/binders_annotations_merged.csv"
ERROR_CSV  = "/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/functional_annotation/binders_annotations_errors.csv"
MERGED_CHAIN_CSV = "/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/functional_annotation/binders_chains_analysis_merged.csv"

LOG_FILE   = "/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/functional_annotation/binders_annotations.log"

NUM_THREADS   = 16
REQUEST_SLEEP = 5

DATA_API = "https://data.rcsb.org/graphql"


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
# GLOBAL CONTAINERS
# =========================

ERROR_RECORDS = []
CHAIN_SUMMARY = []
VALID_PDB_SET = set()

ALLOWED_ORGANISMS = {
    "Homo sapiens",
    "Human adenovirus E serotype 4",
    "Human adenovirus C serotype 6",
    "Human cytomegalovirus (strain AD169)",
    "Echovirus E6",
    "Echovirus E11",
    "Echovirus E30",
    "synthetic construct",
    "Plasmodium falciparum",
    "Plasmodium falciparum HB3",
    "Echovirus E18",
}

def log_error(pdb_id, old_chain, new_chain, message):
    ERROR_RECORDS.append({
        "pdb_id": pdb_id,
        "original_chain": old_chain,
        "new_chain": new_chain,
        "error": message,
    })
    logging.error(f"[{pdb_id}:{new_chain}] {message}")


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
        except requests.RequestException:
            logging.warning(f"{label} attempt {attempt}/{retries} failed")
            time.sleep(2 ** attempt)
    logging.error(f"{label} failed after {retries} attempts")
    return None


# =========================
# GO FALLBACK
# =========================

QUERY = """
query($pdb: String!) {
  entry(entry_id: $pdb) {
    polymer_entities {
      rcsb_polymer_entity_container_identifiers {
        entity_id
        auth_asym_ids
      }
      rcsb_polymer_entity {
        pdbx_description
      }
      rcsb_entity_source_organism {
        ncbi_scientific_name
      }
    }
  }
}
"""

def go_fallback(pdb_id, old_chain):
    try:
        r = requests.post(
            DATA_API,
            json={"query": QUERY, "variables": {"pdb": pdb_id}},
            timeout=10
        )
        r.raise_for_status()

        entry = r.json().get("data", {}).get("entry")
        if not entry:
            return None, None

        for ent in entry.get("polymer_entities", []):
            auth_ids = (
                ent.get("rcsb_polymer_entity_container_identifiers", {})
                .get("auth_asym_ids") or []
            )

            if old_chain in auth_ids:
                polymer = (
                    ent.get("rcsb_polymer_entity", {})
                    .get("pdbx_description", "Missing entry")
                )

                organism_list = ent.get("rcsb_entity_source_organism") or []
                organism = (
                    organism_list[0].get("ncbi_scientific_name")
                    if organism_list else "Missing entry"
                )

                return polymer, organism

        return None, None

    except Exception as e:
        logging.warning(f"GO fallback failed for {pdb_id}:{old_chain} ({e})")
        return None, None


# =========================
# SIFTS + UNIPROT
# =========================

def get_sifts_mapping(pdb_id):
    url = f"https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/{pdb_id.lower()}"
    data = safe_get(url, f"SIFTS:{pdb_id}")
    if not data or pdb_id.lower() not in data:
        return {}
    return data[pdb_id.lower()].get("UniProt", {})


def get_uniprot_metadata(uniprot_id):
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    data = safe_get(url, f"UniProt:{uniprot_id}")
    if not data:
        return None

    name = data.get("proteinDescription", {}) \
               .get("recommendedName", {}) \
               .get("fullName", {}) \
               .get("value", "Missing entry")

    gene_name = "Missing entry"
    try:
        gene_name = data["genes"][0]["geneName"]["value"]
    except Exception:
        pass

    organism = data.get("organism", {}).get("scientificName", "Missing entry")

    return {
        "uniprot_name": name,
        "gene_name": gene_name,
        "organism": organism,
    }


# =========================
# CORE LOGIC
# =========================

def map_entry(pdb_id, new_chain, chain_map):

    pdb_clean = pdb_id.split("-")[0].upper()

    old_chain_raw = chain_map.get(pdb_clean, {}).get(new_chain)
    if not old_chain_raw:
        log_error(pdb_clean, None, new_chain, "chain mapping missing")
        return []

    old_chain = str(old_chain_raw).split("-")[0]

    CHAIN_SUMMARY.append({
        "pdb_id": pdb_clean[:4],
        "original_chain": old_chain,
        "new_chain": new_chain
    })

    mappings = get_sifts_mapping(pdb_clean)

    seen = []
    for uniprot_key, block in mappings.items():
        uniprot_id = uniprot_key.split("_")[0]
        for seg in block.get("mappings", []):
            if seg.get("chain_id") == old_chain:
                meta = get_uniprot_metadata(uniprot_id)
                if meta:
                    seen.append((uniprot_id, meta))
                break

    # 🚨 MULTIPLE MAPPINGS → FORCE GO
    if len(seen) > 1:
        seen = []

    use_go = False

    if not seen:
        use_go = True
    else:
        meta = seen[0][1]
        gene = meta["gene_name"]
        name = meta["uniprot_name"]
        org  = meta["organism"]

        gene_missing = gene in [None, "", "Missing entry"]
        name_missing = name in [None, "", "Missing entry"]
        org_missing  = org in [None, "", "Missing entry"]

        if (gene_missing and name_missing) or org_missing or gene=="Genome polyprotein":

    if use_go:
        polymer, organism = go_fallback(pdb_clean, old_chain)
        if polymer:
            return [{
                "pdb_id": pdb_clean,
                "new_chain": new_chain,
                "original_chain": old_chain_raw,
                "uniprot_id": "GO_FALLBACK",
                "gene_name": "Missing entry",
                "uniprot_name": polymer,
                "organism": organism if organism else "Missing entry",
                "classification": polymer
            }]
        else:
            log_error(pdb_clean, old_chain_raw, new_chain, "GO fallback failed")
            return []

    # NORMAL CASE
    uniprot_id, meta = seen[0]
    gene = meta["gene_name"]
    name = meta["uniprot_name"]

    classification = gene if gene not in [None, "", "Missing entry"] else name

    return [{
        "pdb_id": pdb_clean,
        "new_chain": new_chain,
        "original_chain": old_chain_raw,
        "uniprot_id": uniprot_id,
        "gene_name": gene,
        "uniprot_name": name,
        "organism": meta["organism"],
        "classification": classification
    }]


# =========================
# MAIN
# =========================

def main():

    global VALID_PDB_SET

    # ---------------------------------------
    # LOAD PDB FILTER FIRST (HARD FILTER)
    # ---------------------------------------
    filter_df = pd.read_csv(FILTER_CSV)

    VALID_PDB_SET = set(
        filter_df["pdb"].astype(str).str[:4].str.upper()
    )

    # ---------------------------------------
    # LOAD INPUT AND FILTER EARLY
    # ---------------------------------------
    df = pd.read_csv(INPUT_CSV)
    df = df[df["type"].str.lower() == "binder"]

    # 🚨 HARD SKIP BEFORE ANY THREAD OR API CALL
    df["pdb_clean"] = df["pdb_id"].astype(str).str.split("-").str[0].str.upper()
    df = df[df["pdb_clean"].str[:4].isin(VALID_PDB_SET)]

    if df.empty:
        logging.warning("No PDB entries passed pdb_assemblies_analysis_filtered.csv filter")
        return

    # ---------------------------------------
    # LOAD CHAIN MAP
    # ---------------------------------------
    df_chain_map = pd.read_csv(CHAIN_MAP_CSV)
    chain_map = {}
    for _, row in df_chain_map.iterrows():
        chain_map.setdefault(row["pdb"].upper(), {})[row["new_chain"]] = row["old_chain"]

    results = []

    # ---------------------------------------
    # THREAD EXECUTION
    # ---------------------------------------
    with ThreadPoolExecutor(max_workers=NUM_THREADS) as executor:
        futures = [
            executor.submit(map_entry, row.pdb_id, row.Chain, chain_map)
            for row in df.itertuples(index=False)
        ]
        for f in as_completed(futures):
            results.extend(f.result())

    # ---------------------------------------
    # APPLY ORGANISM FILTER FIRST
    # ---------------------------------------
    if results:
        results_df = pd.DataFrame(results)

        results_df = results_df[
            results_df["organism"].isin(ALLOWED_ORGANISMS)
        ]
        
        if not results_df.empty:
            results_df.to_csv(OUTPUT_CSV, index=False)

            # ---------------------------------------
            # NOW BUILD CHAIN SUMMARY FROM FILTERED DATA
            # ---------------------------------------
            chain_summary_df = (
                results_df[["pdb_id", "original_chain", "new_chain"]]
                .drop_duplicates()
            )
        
            chain_summary_df.to_csv(MERGED_CHAIN_CSV, index=False)

    # ---------------------------------------
    # ERRORS
    # ---------------------------------------
    if ERROR_RECORDS:
        pd.DataFrame(ERROR_RECORDS).to_csv(ERROR_CSV, index=False)

if __name__ == "__main__":
    main()
