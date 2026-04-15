#!/bin/python3
# ---------------------
# DSSP-safe version (Biopython >=1.85)
# Uses PDBs directly, no temp files
# ---------------------

import os
import pandas as pd
from Bio.PDB import PDBParser, PDBIO, Select, NeighborSearch, is_aa
from Bio.PDB.DSSP import DSSP

# === CONFIGURABLE PATHS ===
DATABASE        = "pdb-teste"
INPUT_CSV       = f"/mnt/4TB/giovanna/foldseek/version_02/filter/step3/pdb/pdb_assemblies_remapped.csv"
PDB_DIR         = f"/mnt/4TB/giovanna/foldseek/version_02/filter/step4/pdb/trimmed_mhc"
OUTPUT_PDB_DIR  = f"/mnt/4TB/giovanna/foldseek/version_02/filter/step5/{DATABASE}/trimmed_mhc_complexes"

CONTACTS_CSV    = f"/mnt/4TB/giovanna/foldseek/version_02/filter/step5/{DATABASE}/mhc_contacts.csv"
LOG_CSV         = f"/mnt/4TB/giovanna/foldseek/version_02/filter/step5/{DATABASE}/mhc_log.csv"
BLACKLIST_CSV   = f"/mnt/4TB/giovanna/foldseek/version_02/filter/step5/blacklist.csv"

DISTANCE_CUTOFF = 4.0
MIN_BINDER_AA   = 30
HELIX_CODES     = {"H", "G", "I"}

# ============================================================
# === PDB chain selector (RULE-BASED)
# ============================================================
class ComplexSelect(Select):
    def __init__(self, mhc_chain_id, binder_chains, ligand_chains, blacklist):
        self.allowed_chains = {mhc_chain_id} | binder_chains | ligand_chains
        self.blacklist = blacklist

    def accept_chain(self, chain):
        return chain.id in self.allowed_chains

    def accept_residue(self, residue):
        # Only exclusion rule
        if residue.get_resname() in self.blacklist:
            return False
        return True

# ============================================================
# === DSSP helix detection (UNCHANGED)
# ============================================================
def dssp_helix_residues(pdb_file, mhc_chain_id, start, end, pdb_id):

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_file)
    model = structure[0]

    helix_res = set()
    try:
        dssp = DSSP(model, pdb_file)
    except Exception as e:
        print(f"[DSSP FAIL] {pdb_id} → {e}")
        return None

    for key in dssp.keys():
        chain_id, res_id = key
        ss = dssp[key][2]

        if ss not in HELIX_CODES:
            continue
        if chain_id != mhc_chain_id:
            continue

        resseq = res_id[1]
        if start <= resseq <= end:
            try:
                helix_res.add(model[chain_id][res_id])
            except KeyError:
                pass

    if not helix_res:
        print(f"[WARN] {pdb_id} → DSSP found NO helices in window {mhc_chain_id}:{start}-{end}")
        return set()

    return helix_res

# ============================================================
# === Load blacklist
# ============================================================
def load_blacklist(csv_path):
    if not os.path.exists(csv_path):
        return set()
    df = pd.read_csv(csv_path)
    return set(df["resname"].astype(str))

# ============================================================
# === Contact detection (UNCHANGED)
# ============================================================
def classify_contact_chains(structure, helix_residues, mhc_chain_id):

    helix_atoms = [a for r in helix_residues for a in r.get_atoms()]
    ns = NeighborSearch(list(structure.get_atoms()))

    contacting = set()

    for atom in helix_atoms:
        for close in ns.search(atom.coord, DISTANCE_CUTOFF):
            res = close.get_parent()
            if not is_aa(res, standard=True):
                continue
            chain = res.get_parent()
            if chain.id != mhc_chain_id:
                contacting.add(chain.id)

    binders, ligands = set(), set()

    for cid in contacting:
        chain = structure[0][cid]
        aa_count = sum(1 for r in chain if is_aa(r, standard=True))
        (binders if aa_count >= MIN_BINDER_AA else ligands).add(cid)

    return binders, ligands

# ============================================================
# === Process single entry
# ============================================================
def process_entry(pdb_file, pdb_id, chain_id, start, end, blacklist):

    structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)

    if chain_id not in structure[0]:
        return [("SKIPPED", pdb_id, chain_id, "ChainNotFound")]

    helix_residues = dssp_helix_residues(pdb_file, chain_id, start, end, pdb_id)

    if helix_residues is None:
        return [("SKIPPED", pdb_id, chain_id, "DSSPFailure")]

    if not helix_residues:
        return [("SKIPPED", pdb_id, chain_id, "NoHelixDSSP")]

    binder_chains, ligand_chains = classify_contact_chains(
        structure, helix_residues, chain_id
    )

    selector = ComplexSelect(chain_id, binder_chains, ligand_chains, blacklist)

    io = PDBIO()
    io.set_structure(structure)
    io.save(
        os.path.join(OUTPUT_PDB_DIR, f"{pdb_id}_filtered_complex.pdb"),
        selector
    )

    results = [("MHC", pdb_id, chain_id, chain_id)]

    for c in binder_chains:
        results.append(("BINDER", pdb_id, chain_id, c))
    for c in ligand_chains:
        results.append(("LIGAND", pdb_id, chain_id, c))

    if not binder_chains and not ligand_chains:
        results.append(("SKIPPED", pdb_id, chain_id, "NoContacts"))

    return results

# ============================================================
# === Main
# ============================================================
def main(csv_path, pdb_dir, pdb_output_dir):

    os.makedirs(pdb_output_dir, exist_ok=True)

    df = pd.read_csv(csv_path)
    blacklist = load_blacklist(BLACKLIST_CSV)

    all_contacts = []

    df_sel = df[df["status"].str.lower() == "primary"]

    for _, row in df_sel.iterrows():
        pdb_id = row["pdb"].split("-")[0]
        pdb_file = os.path.join(pdb_dir, f"{pdb_id}_trimmed_mhc.pdb")

        if not os.path.exists(pdb_file):
            all_contacts.append(("SKIPPED", pdb_id, "A", "MissingPDB"))
            continue

        try:
            all_contacts.extend(
                process_entry(
                    pdb_file, pdb_id, "A",
                    int(row["tstart"]), int(row["tend"]),
                    blacklist
                )
            )
        except Exception as e:
            print(f"[ERROR] {pdb_id} → {e}")
            all_contacts.append(("SKIPPED", pdb_id, "A", "RuntimeError"))
            continue

    # === Unified contacts CSV ===
    contacts = [
        (pdb, chain, tag.lower())
        for tag, pdb, _, chain in all_contacts
        if tag in {"MHC", "BINDER", "LIGAND"}
    ]

    pd.DataFrame(
        contacts,
        columns=["pdb_id", "Chain", "type"]
    ).drop_duplicates().to_csv(CONTACTS_CSV, index=False)

    # === Log CSV ===
    log = {}

    for tag, pdb, _, value in all_contacts:
        log.setdefault(pdb, {"binder": "no", "ligand": "no", "comment": ""})

        if tag == "BINDER":
            log[pdb]["binder"] = "yes"
        elif tag == "LIGAND":
            log[pdb]["ligand"] = "yes"
        elif tag == "SKIPPED":
            log[pdb]["binder"] = "NA"
            log[pdb]["ligand"] = "NA"
            log[pdb]["comment"] = value

    pd.DataFrame(
        [(pdb, v["binder"], v["ligand"], v["comment"]) for pdb, v in log.items()],
        columns=["pdb_id", "binder", "ligand", "comment"]
    ).to_csv(LOG_CSV, index=False)

    print("[DONE] Step5 finished")

if __name__ == "__main__":
    main(INPUT_CSV, PDB_DIR, OUTPUT_PDB_DIR)
