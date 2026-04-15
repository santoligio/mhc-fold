import os
import glob
import pandas as pd
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.PDB import MMCIFParser, PDBIO, Select
from Bio.PDB.DSSP import DSSP
from Bio.PDB.NeighborSearch import NeighborSearch
from tqdm import tqdm


class TrimmedChainSelect(Select):
    """
    Trim a specific chain to a residue range.
    Keep all other chains except duplicates.
    """

    def __init__(self, chain_id, start_res, end_res, duplicate_chains=None):
        self.chain_id = chain_id
        self.start_res = start_res
        self.end_res = end_res
        self.allowed_residues = set()
        self.extra_residues = set()
        self.duplicate_chains = set(duplicate_chains or [])

    def accept_chain(self, chain):
        return chain.id not in self.duplicate_chains

    def accept_residue(self, residue):
        chain_id = residue.get_parent().id
        resseq = residue.id[1]

        # Keep trimmed residues
        if chain_id == self.chain_id:
            if resseq in self.allowed_residues:
                return True

        # Keep extra residues near helices
        if (chain_id, resseq) in self.extra_residues:
            return True

        # Keep everything else (other chains)
        return True

    def collect_allowed_residues(self, structure):
        model = structure[0]
        if self.chain_id not in model:
            return

        chain = model[self.chain_id]
        for residue in chain:
            resseq = residue.id[1]
            if self.start_res <= resseq <= self.end_res:
                self.allowed_residues.add(resseq)

def process_structure(
    cif_path,
    chain_id,
    start_res,
    end_res,
    output_path,
    duplicate_chains=None
):
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("structure", cif_path)

    selector = TrimmedChainSelect(
        chain_id=chain_id,
        start_res=start_res,
        end_res=end_res,
        duplicate_chains=duplicate_chains
    )
    selector.collect_allowed_residues(structure)

    model = structure[0]

    ##########################################################
    extra_residues = set()
    
    for chain_id in model:
        if chain_id == "A":
            for residue in chain_id:
            if residue not in selector.allowed_residues:
                if not is_aa(res, standard=True):
                    extra_residues.add(residue)

    selector.extra_residues = extra_residues

    ##########################################################
    
    if not selector.allowed_residues:
        return f"[FAILED] {os.path.basename(output_path)} | no residues in range"

    expected = end_res - start_res + 1
    n_trimmed = len(selector.allowed_residues)

    if n_trimmed != expected:
        return (
            f"[FAILED] {os.path.basename(output_path)} | "
            f"range mismatch {n_trimmed}/{expected}"
        )

    io = PDBIO()
    io.set_structure(structure)

    with open(output_path, "w") as f:
        f.write(f"HEADER    {os.path.basename(output_path)}\n")
        io.save(f, selector)

    return f"[SAVED] {output_path} | Trimmed residues: {n_trimmed}"


def load_duplicate_chain_map(csv_file):
    df = pd.read_csv(csv_file)
    dup_map = {}

    for _, row in df.iterrows():
        if str(row.get("status", "")).lower() != "duplicate":
            continue

        pdb_id = row["pdb"].split("-")[0].lower()
        chain_id = row["chain"]
        dup_map.setdefault(pdb_id, set()).add(chain_id)

    return dup_map


def extract_structure_id(target):
    parts = target.split("-")

    if parts[0].upper() == "AF":
        return parts[1].upper(), "AFDB"

    return parts[0].lower(), "PDB"


def find_input_cif(input_cif_dir, struct_id, source_hint):
    if source_hint == "AFDB":
        af_matches = glob.glob(
            f"{input_cif_dir}/AF-{struct_id}-*-model_*.cif"
        )
        if af_matches:
            return af_matches[0], "AFDB"

    pdb_matches = glob.glob(
        f"{input_cif_dir}/{struct_id}*_remapped.cif"
    )
    if pdb_matches:
        return pdb_matches[0], "PDB"

    return None, None


def worker(row, input_cif_dir, output_dir, duplicate_map):
    struct_id, source_hint = extract_structure_id(row["pdb"])

    tstart = int(row["tstart"])
    tend   = int(row["tend"])

    cif_path, source = find_input_cif(
        input_cif_dir, struct_id, source_hint
    )
    if cif_path is None:
        return (
            f"[MISSING] No CIF found | source={source_hint} | "
            f"id={struct_id} | target={row['pdb']}"
        )

    output_path = os.path.join(
        output_dir, f"{struct_id}_trimmed_mhc.pdb"
    )

    if source == "AFDB":
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("af", cif_path)
        chain_id = next(structure[0].get_chains()).id

        return process_structure(
            cif_path,
            chain_id,
            tstart,
            tend,
            output_path,
            duplicate_chains=None
        )

    dup_chains = duplicate_map.get(struct_id, set())
    chain_id = row["chain"]

    return process_structure(
        cif_path,
        chain_id,
        tstart,
        tend,
        output_path,
        duplicate_chains=dup_chains
    )


def run_trimming_multithreaded(
    csv_file,
    input_cif_dir,
    output_dir,
    max_workers=8
):
    os.makedirs(output_dir, exist_ok=True)

    df = pd.read_csv(csv_file)
    duplicate_map = load_duplicate_chain_map(csv_file)
    primary_df = df[df["status"].str.lower() == "primary"]

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = [
            executor.submit(
                worker,
                row,
                input_cif_dir,
                output_dir,
                duplicate_map
            )
            for _, row in primary_df.iterrows()
        ]

        for future in tqdm(
            as_completed(futures),
            total=len(futures),
            desc="Trimming structures"
        ):
            print(future.result())


if __name__ == "__main__":

    csv_file = (
        "/mnt/4TB/giovanna/foldseek/version_02/filter/step3/"
        "pdb/pdb_assemblies_remapped.csv"
    )
    input_cif_dir = (
        "/mnt/4TB/giovanna/foldseek/version_02/filter/step3/pdb/remapped_cifs"
    )
    output_dir = (
        "/mnt/4TB/giovanna/foldseek/version_02/filter/step4/"
        "pdb-teste/trimmed_mhc"
    )

    run_trimming_multithreaded(
        csv_file,
        input_cif_dir,
        output_dir,
        max_workers=16
    )
