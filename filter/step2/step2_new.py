#!/bin/python3
# ------------------------------
# Conceptualized by Helder Filho
# Modified by Elton Chaves
# Simplified: download + gunzip + remove .gz
# ------------------------------

import os
import gzip
import shutil
import requests
import pandas as pd
from concurrent.futures import ThreadPoolExecutor


def gunzip_and_remove(gz_path):
    """Extract .gz file and delete the original."""
    out_path = gz_path[:-3]  # remove .gz

    with gzip.open(gz_path, 'rb') as f_in:
        with open(out_path, 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)

    os.remove(gz_path)
    print(f"Extracted and removed: {gz_path} -> {out_path}")
    return out_path


def download_assembly(input_csv, pdb_id, output_dir, assembly_number):
    csv_name = os.path.basename(input_csv).lower()

    if 'pdb' in csv_name:
        url = f'https://files.rcsb.org/download/{pdb_id.upper()}-assembly{assembly_number}.cif.gz'
        output_path = os.path.join(
            output_dir, '1_assemblies',
            f'{pdb_id}-assembly{assembly_number}.cif.gz'
        )

        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        if not os.path.exists(output_path):
            r = requests.get(url)

            # -------- FAILSAFE FALLBACK --------
            if r.status_code != 200:
                print(f"Assembly failed, trying entry CIF: {pdb_id}")
                fallback_url = f'https://files.rcsb.org/download/{pdb_id.upper()}.cif.gz'
                output_path = os.path.join(
                    output_dir, '1_assemblies',
                    f'{pdb_id}-assembly1.cif.gz'
                )
                r = requests.get(fallback_url)

            if r.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(r.content)
                print(f"Downloaded: {output_path}")
            else:
                print(f"Failed to download assembly and fallback for: {pdb_id}")
                return None

        else:
            print(f"File exists: {output_path}")

    elif 'afdb' in csv_name:
        url = f'https://alphafold.ebi.ac.uk/files/{pdb_id}.cif'
        output_path = os.path.join(
            output_dir, '1_models',
            f'{pdb_id}.cif'
        )

        os.makedirs(os.path.dirname(output_path), exist_ok=True)

        if not os.path.exists(output_path):
            r = requests.get(url)
            if r.status_code == 200:
                with open(output_path, 'wb') as f:
                    f.write(r.content)
                print(f"Downloaded: {output_path}")
            else:
                print(f"Failed to download: {url}")
                return None
        else:
            print(f"File exists: {output_path}")

    else:
        print(f"Unknown CSV type: {input_csv}")
        return None

    # auto-extract if gz
    if output_path.endswith(".gz"):
        return gunzip_and_remove(output_path)

    return output_path


def process_csv(input_csv, output_dir, threads=8):
    df = pd.read_csv(input_csv)
    os.makedirs(output_dir, exist_ok=True)

    csv_name = os.path.basename(input_csv).lower()

    with ThreadPoolExecutor(max_workers=threads) as executor:

        if 'pdb' in csv_name:
            futures = []
            for _, row in df.iterrows():

                if str(row.get('status', '')) != 'primary':
                    continue

                pdb_assembly = row['pdb']
                pdb_id, assembly_part = pdb_assembly.split('-')
                assembly_number = (
                    assembly_part.replace('assembly', '')
                    .replace('.cif', '')
                )

                futures.append(
                    executor.submit(
                        download_assembly,
                        input_csv,
                        pdb_id,
                        output_dir,
                        assembly_number
                    )
                )

            for f in futures:
                f.result()

        elif 'afdb' in csv_name:
            futures = []
            for _, row in df.iterrows():
                pdb_model = row['pdb']
                futures.append(
                    executor.submit(
                        download_assembly,
                        input_csv,
                        pdb_model,
                        output_dir,
                        None
                    )
                )

            for f in futures:
                f.result()

        else:
            print(f"Unrecognized CSV type: {input_csv}")


if __name__ == "__main__":
    input_csv = {
        'pdb': '/mnt/4TB/giovanna/foldseek/version_02/filter/step1/pdb/pdb_assemblies.csv'
    }

    fs_inp = input_csv['pdb']

    output_dir = "/mnt/4TB/giovanna/foldseek/version_02/filter/step2/pdb-teste"
    process_csv(fs_inp, output_dir, threads=8)
