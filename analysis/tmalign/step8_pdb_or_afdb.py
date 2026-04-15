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
MODE = "afdb"   # options: "afdb" or "pdb"
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
        pdb_csv,
        num_workers=None,
        chunk_size=10000):

    pdb_df = pd.read_csv(pdb_csv)

    # -------------------------
    # MODE-SPECIFIC FILTERING
    # -------------------------
    if MODE == "afdb":
        allowed_ids = set(
            pdb_df["pdb"].map(extract_uniprot_id_from_csv).dropna()
        )

        def file_allowed(fname):
            return extract_uniprot_id_from_file(fname) in allowed_ids

    elif MODE == "pdb":
        allowed_ids = set(
            pdb_df["pdb"].apply(extract_base_pdb).dropna()
        )

        def file_allowed(fname):
            base = extract_base_pdb(os.path.splitext(fname)[0])
            return base in allowed_ids

    else:
        raise ValueError("MODE must be 'afdb' or 'pdb'")

    # -------------------------
    # File filtering
    # -------------------------
    all_files = [f for f in os.listdir(input_dir) if f.endswith(".pdb")]
    files = sorted([f for f in all_files if file_allowed(f)])

    n = len(files)
    total_pairs = n * (n - 1) // 2

    print(f"Mode: {MODE.upper()}")
    print(f"Found {len(all_files)} total PDB files")
    print(f"Using {n} files after CSV filtering")
    print(f"{total_pairs:,} total comparisons")

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
        input_dir="/mnt/4TB/giovanna/foldseek/version_02/filter/step4/afdb/trimmed_mhc",
        output_prefix="/mnt/4TB/giovanna/foldseek/version_02/analysis/tmalign/filtered_afdb_avg/mhc_all-vs-all",
        pdb_csv="/mnt/4TB/giovanna/foldseek/version_02/analysis/functional_annotation/filtered/afdb_models_analysis_filtered.csv",
        num_workers=14,
        chunk_size=10000
    )
