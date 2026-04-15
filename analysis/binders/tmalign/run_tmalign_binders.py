#!/usr/bin/env python3

import os
import subprocess
import pandas as pd
from itertools import combinations
from multiprocessing import Pool, cpu_count
from tqdm import tqdm


# ------------------------------------------------------------
# Filename extraction for binder files
# ------------------------------------------------------------
def extract_pdb_and_chain_from_filename(filename):
    """
    Example:
    3qdm_binder_chainD.pdb → (3QDM, D)
    """
    name = os.path.splitext(filename)[0]

    parts = name.split("_")
    if len(parts) < 3:
        return None, None

    pdb_id = parts[0].upper()

    chain_part = parts[-1]  # chainD
    if not chain_part.lower().startswith("chain"):
        return None, None

    chain = chain_part.replace("chain", "").replace("CHAIN", "")
    return pdb_id, chain


# ------------------------------------------------------------
# TM-align runner
# ------------------------------------------------------------
def run_tmalign(pair):
    file1, file2, input_dir = pair
    path1 = os.path.join(input_dir, file1)
    path2 = os.path.join(input_dir, file2)

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
        input_dir,
        output_prefix,
        binder_csv,
        num_workers=None,
        chunk_size=10000):

    binder_df = pd.read_csv(binder_csv)

    # Ensure uppercase
    binder_df["pdb_id"] = binder_df["pdb_id"].str.upper()
    binder_df["new_chain"] = binder_df["new_chain"].astype(str)

    # Create allowed (pdb, chain) pairs
    allowed_pairs = set(
        zip(binder_df["pdb_id"], binder_df["new_chain"])
    )

    # -------------------------
    # File filtering
    # -------------------------
    all_files = [f for f in os.listdir(input_dir) if f.endswith(".pdb")]

    files = []
    for f in all_files:
        pdb_id, chain = extract_pdb_and_chain_from_filename(f)
        if pdb_id is None:
            continue

        if (pdb_id, chain) in allowed_pairs:
            files.append(f)

    files = sorted(files)

    n = len(files)
    total_pairs = n * (n - 1) // 2

    print("Mode: BINDER")
    print(f"Found {len(all_files)} total PDB files")
    print(f"Using {n} binder files after binder CSV filtering")
    print(f"{total_pairs:,} total comparisons")

    if n < 2:
        print("Not enough files to compare.")
        return

    if num_workers is None:
        num_workers = max(1, cpu_count() - 1)

    pair_generator = generate_all_pairs(files)
    global_pbar = tqdm(total=total_pairs, desc="Global progress", position=0)

    chunk_id = 0
    for chunk in chunked_iterator(pair_generator, chunk_size):
        chunk_pairs = [(f1, f2, input_dir) for f1, f2 in chunk]

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
        input_dir="/mnt/4TB/giovanna/foldseek/version_02/filter/step6/pdb/binders",
        output_prefix="/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/tmalign/binders_all-vs-all",
        binder_csv="/mnt/4TB/giovanna/foldseek/version_02/analysis/binders/functional_annotation/binders_chains_analysis_filtered.csv",
        num_workers=14,
        chunk_size=10000
    )
