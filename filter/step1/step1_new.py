#!/bin/python3
# ---------------------------
# Written by Elton Chaves
# E-mail: chavesejf@gmail.com
# ---------------------------

import os
import pandas as pd

# ---
fs_db = {
    'pdb':  '/mnt/4TB/giovanna/foldseek/version_02/pdb/dbs_pdb_aln',
    #'afdb': '/mnt/4TB/giovanna/foldseek/version_02/alphafold/3mre_afdb_aln'
}

# ---
fs_inp = fs_db['pdb']

# ---
df         = pd.read_csv(fs_inp, delimiter="\\s+")
df.columns = [
    'query','target','fident','alnlen','mismatch','gapopen',
    'qstart','qend','tstart','tend','evalue','bits',
    'alntmscore','qtmscore','ttmscore','lddt','lddtfull','prob'
]

# Filtros iniciais
df_filtered = df[
    (df['evalue'] <= 0.01) &
    (df['alnlen'] >= 180 / 2)
]

# Garante organização dos targets em ordem de melhor alinhamento e remove repetições
df_filtered = (
    df_filtered
    .sort_values('evalue')
    .drop_duplicates(subset='target', keep='first')
)

# Lista para armazenar resultados
records = []
seen_pdbs = set()
all_targets = set(df_filtered['target'])

# Loop otimizado
for _, row in df_filtered.iterrows():
    entry       = row['target']
    split_entry = '-assembly' in entry

    # Remove PDBs que não são assembly1 
    if split_entry:
        pdb_id = entry.split('_')[0]
        chain_id = entry.split('_')[1]

        if not pdb_id.endswith('-assembly1'):
            pdb_code = pdb_id.split('-')[0]
            base_prefix = f"{pdb_code}-assembly1"
            if any(t.startswith(base_prefix) for t in all_targets):
                continue

    # Se existe o modelo principal, remove os demais
    else:
        pdb_id   = entry
        chain_id = 'NoChainInfo'

        parts = entry.split('-')
        if len(parts) == 5:
            base_entry = f"{parts[0]}-{parts[1]}-{parts[3]}-{parts[4]}"

            if base_entry in all_targets:
                continue
        
    # Guarda as informações de estruturas sem repetição 
    try:
        if pdb_id not in seen_pdbs:
            seen_pdbs.add(pdb_id)
            status = "primary"
  
        else:
            status = "duplicate"

        records.append({
            'pdb': pdb_id,
            'chain': chain_id,
            'tstart': row['tstart'],
            'tend': row['tend'],
            'status': status
        })

    except IndexError:
        print(f"{pdb_id} - Incorrect formatting")
        continue  # ignora entradas mal formatadas

# Converter lista para DataFrame
pdb_files = pd.DataFrame.from_records(records)
                
# Salvar
pdb_files.to_csv(
    '/mnt/4TB/giovanna/foldseek/version_02/filter/step1/pdb-teste/pdb_assemblies.csv',
    index=False)
