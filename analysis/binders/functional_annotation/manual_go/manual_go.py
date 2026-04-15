#!/usr/bin/env python3

import requests
import argparse

DATA_API = "https://data.rcsb.org/graphql"

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

def query_pdb(pdb_id):
    r = requests.post(
        DATA_API,
        json={"query": QUERY, "variables": {"pdb": pdb_id}},
        timeout=10
    )
    r.raise_for_status()
    return r.json()["data"]["entry"]

def find_chain(entry, chains_to_try):
    for chain in chains_to_try:
        if chain is None:
            continue

        for ent in entry.get("polymer_entities", []):
            ids = ent.get("rcsb_polymer_entity_container_identifiers", {})
            auth_ids = ids.get("auth_asym_ids") or []

            if chain in auth_ids:
                desc = ent.get("rcsb_polymer_entity", {}).get(
                    "pdbx_description", "Missing entry"
                )

                org_list = ent.get("rcsb_entity_source_organism") or []
                organism = (
                    org_list[0].get("ncbi_scientific_name")
                    if org_list else "Missing entry"
                )

                return chain, desc, organism

    return None, None, None

def main():
    parser = argparse.ArgumentParser(
        description="Query PDB entry and retrieve chain description and organism."
    )

    parser.add_argument(
        "pdb_id",
        type=str,
        help="PDB ID (e.g., 6LA6)"
    )

    parser.add_argument(
        "chain1",
        type=str,
        help="Primary chain to search (e.g., EEE or E)"
    )

    parser.add_argument(
        "--chain2",
        type=str,
        default=None,
        help="Optional secondary chain fallback"
    )

    args = parser.parse_args()

    pdb_id = args.pdb_id.upper()
    chains_to_try = [args.chain1]

    if args.chain2:
        chains_to_try.append(args.chain2)

    entry = query_pdb(pdb_id)
    chain, desc, organism = find_chain(entry, chains_to_try)

    if chain:
        print(f"Found chain: {chain}")
        print(f"Description: {desc}")
        print(f"Organism: {organism}")
    else:
        print(f"No matching chain found ({', '.join(filter(None, chains_to_try))})")

if __name__ == "__main__":
    main()