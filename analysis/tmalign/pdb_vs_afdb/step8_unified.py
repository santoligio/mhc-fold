#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
from itertools import combinations
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

# ============================================================
# USER EDIT HERE
# ============================================================
MODE = "mixed"   # now controlled by afdb_vs_pdb.csv
# ============================================================


# ------------------------------------------------------------
# ID extraction helpers
# ------------------------------------------------------------
def extract_uniprot_id_from_csv(text):
    """AF-Q9MYF8-F1-model_v6 → Q9MYF8"""
    if pd.isna(text):
        return None
    parts = str(text).split("-")
    return parts[1].upper() if len(parts) > 1 else None


def extract_uniprot_id_from_file(text):
    """AF-Q9MYF8_trimmed_mhc.pdb → Q9MYF8"""
    if pd.isna(text):
        return None
    text = os.path.splitext(text)[0]
    text = text.split("_")[0]
    parts = text.split("-")
    return parts[1].upper() if len(parts) > 1 else None


def extract_base_pdb(text):
    """7ABC-assembly1 → 7ABC"""
    if pd.isna(text):
        return None
    return str(text).split("-")[0][:4].upper()


# ------------------------------------------------------------
# TM-align runner
# ------------------------------------------------------------
def run_tmalign(pair):
    file1, file2, input_dir1, input_dir2 = pair
    path1 = os.path.join(input_dir1, file1)
    path2 = os.path.join(input_dir2, file2)

    if not os.path.isfile(path1) or not os.path.isfile(path2):
        return (file1, file2, None)

    try:
        output = subprocess.check_output(
            ["TMalign", path1, path2, "-a", "T"],
            stderr=subprocess.DEVNULL
        ).decode()

        tm_count = 0

        for line in output.splitlines():
            if line.strip().startswith("TM-score="):
                tm_count += 1

                if tm_count == 3:
                    tm_score = float(line.split("=")[1].split()[0])
                    return (file1, file2, tm_score)

        return (file1, file2, None)

    except subprocess.CalledProcessError:
        return (file1, file2, None)

    except subprocess.CalledProcessError:
        return (file1, file2, None)


# ------------------------------------------------------------
def generate_all_pairs(files):
    return combinations(files, 2)


def chunked_iterator(iterable, chunk_size):
    chunk = []
    for item in iterable:
        chunk.append(item)
        if len(chunk) == chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk


# ------------------------------------------------------------
# MAIN FUNCTION
# ------------------------------------------------------------
def parallel_tmalign(
        afdb_dir,
        pdb_dir,
        output_prefix,
        selection_csv,
        num_workers=None,
        chunk_size=10000):

    selection_df = pd.read_csv(selection_csv)

    # -------------------------
    # Split AFDB and PDB entries
    # -------------------------
    afdb_entries = selection_df[selection_df["type"] == "afdb"]
    pdb_entries = selection_df[selection_df["type"] == "pdb"]

    # -------------------------
    # Collect AFDB files
    # -------------------------
    afdb_ids = set(
        afdb_entries["pdb"].map(extract_uniprot_id_from_csv).dropna()
    )

    afdb_files = [
        f for f in os.listdir(afdb_dir)
        if f.endswith(".pdb")
        and extract_uniprot_id_from_file(f) in afdb_ids
    ]

    # -------------------------
    # Collect PDB files
    # -------------------------
    pdb_ids = set(
        pdb_entries["pdb"].apply(extract_base_pdb).dropna()
    )

    pdb_files = [
        f for f in os.listdir(pdb_dir)
        if f.endswith(".pdb")
        and extract_base_pdb(os.path.splitext(f)[0]) in pdb_ids
    ]

    # -------------------------
    # Merge into one structure list
    # Keep track of origin folder
    # -------------------------
    files = []

    for f in afdb_files:
        files.append((f, afdb_dir))

    for f in pdb_files:
        files.append((f, pdb_dir))

    n = len(files)
    total_pairs = n * (n - 1) // 2

    print(f"Mode: MIXED AFDB + PDB")
    print(f"Using {n} total structures")
    print(f"{total_pairs:,} total comparisons")

    if num_workers is None:
        num_workers = max(1, cpu_count() - 1)

    pair_generator = combinations(files, 2)
    global_pbar = tqdm(total=total_pairs, desc="Global progress", position=0)

    chunk_id = 0
    for chunk in chunked_iterator(pair_generator, chunk_size):

        chunk_pairs = [
            (f1[0], f2[0], f1[1], f2[1])
            for f1, f2 in chunk
        ]

        with Pool(processes=num_workers) as pool:
            chunk_results = []

            with tqdm(
                total=len(chunk_pairs),
                desc=f"Chunk {chunk_id}",
                position=1,
                leave=False,
                dynamic_ncols=True
            ) as chunk_pbar:

                for result in pool.imap_unordered(
                        run_tmalign,
                        chunk_pairs,
                        chunksize=100):

                    chunk_results.append(result)
                    global_pbar.update()
                    chunk_pbar.update()

        df_chunk = pd.DataFrame(
            chunk_results,
            columns=["file1", "file2", "tm_score_avg"]
        )

        chunk_filename = f"{output_prefix}_part{chunk_id:04d}.csv"
        df_chunk.to_csv(chunk_filename, index=False)

        print(f"[✓] Saved chunk {chunk_id} to {chunk_filename}")
        chunk_id += 1

    global_pbar.close()
    print(f"\n✅ Finished all {chunk_id} chunks.")


# ------------------------------------------------------------
# Example usage
# ------------------------------------------------------------
if __name__ == "__main__":
    parallel_tmalign(
        afdb_dir="/mnt/4TB/giovanna/foldseek/version_02/filter/step4/afdb/trimmed_mhc",
        pdb_dir="/mnt/4TB/giovanna/foldseek/version_02/filter/step6/pdb/mhc_only",
        output_prefix="/mnt/4TB/giovanna/foldseek/version_02/analysis/tmalign/pdb_vs_afdb/pdb_vs_afdb",
        selection_csv="/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/new_filters/pdb_vs_afdb.csv",
        num_workers=14,
        chunk_size=10000
    )
