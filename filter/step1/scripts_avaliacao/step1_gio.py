#!/bin/python3
# ---------------------------
# Written by Elton Chaves
# E-mail: chavesejf@gmail.com
# ---------------------------

import os
import pandas as pd

# ---
fs_db = {
    'pdb': '/mnt/4TB/giovanna/foldseek/version_02/pdb/dbs_pdb_aln',
}

fs_inp = fs_db['pdb']

base_name  = os.path.basename(fs_inp)
odir_name  = os.path.dirname(fs_inp)

df = pd.read_csv(fs_inp, delimiter=r"\s+")
df.columns = [
    'query','target','fident','alnlen','mismatch','gapopen',
    'qstart','qend','tstart','tend','evalue','bits',
    'alntmscore','qtmscore','ttmscore','lddt','lddtfull','prob'
]

# Filtros iniciais (inalterados)
df_filtered = df[
    (df['evalue'] <= 0.01) &
    (df['alnlen'] >= 180 / 2)]

df_filtered = df_filtered.drop_duplicates(subset='target', keep='first')

records = []
pdb_list = set()   # agora como set (sem alterar lógica)

for _, row in df_filtered.iterrows():
    entry = row['target']
    split_entry = 'pdb' in odir_name

    if split_entry:
        try:
            pdb_assembly, chain_id = entry.split('_')
        except ValueError:
            continue

        pdb_id = pdb_assembly.split('-')[0]

        # filtros estruturais PRIMEIRO
        if len(chain_id) != 1:
            continue

        # deduplicação SOMENTE após entrada válida
        if pdb_id in pdb_list:
            continue

        pdb_list.add(pdb_id)

    else:
        pdb_id = entry
        chain_id = 'NoChainInfo'

        if pdb_id in pdb_list:
            continue

        pdb_list.add(pdb_id)

    records.append({
        'pdb': pdb_id,
        'chain': chain_id,
        'tstart': row['tstart'],
        'tend': row['tend']
    })

pdb_files = pd.DataFrame.from_records(records)

pdb_files.to_csv(
    '/mnt/4TB/giovanna/foldseek/version_02/filter/step1/avaliacao_scripts/gio_assemblies.csv',
    index=False
)
