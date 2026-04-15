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
    #'afdb': '/media/eltonjfc/files-01_nvme2tb/projetos/helderv/2025-pmhc_fold/_pbv/dbf_afup50min_aln'
    }

# ---
fs_inp = fs_db['pdb']

# ---
base_name  = os.path.basename(fs_inp)
odir_name  = os.path.dirname(fs_inp)
df         = pd.read_csv(fs_inp, delimiter="\s+")
df.columns = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend','evalue',
              'bits','alntmscore','qtmscore','ttmscore','lddt','lddtfull','prob']

# Filtros iniciais
df_filtered = df[(df['qtmscore'] > 0.5) &
                 (df['evalue'] <= 0.01) &
                 (df['alnlen'] >= 180 / 2)]
df_filtered = df_filtered.drop_duplicates(subset='target', keep='first')

# Lista para armazenar resultados
records = []

# Loop otimizado
for _, row in df_filtered.iterrows():
    entry       = row['target']
    split_entry = odir_name.__contains__('pdb')

    if split_entry: 
        pdb_id   = entry.split('_')[0]
        chain_id = entry.split('_')[1]
        
        if not pdb_id.endswith('-assembly1'):
            continue
        
        if len(chain_id) != 1:
            continue
    
    else:
        pdb_id   = entry
        chain_id = 'NoChainInfo'

    try:
        records.append({
            'pdb': pdb_id,
            'chain': chain_id, 
            'tstart': row['tstart'],
            'tend': row['tend']
        })
    except IndexError:
        continue  # ignora entradas mal formatadas

# Converter lista para DataFrame
pdb_files = pd.DataFrame.from_records(records)

# Salvar
pdb_files.to_csv(f'/mnt/4TB/giovanna/foldseek/version_02/analysis/step1/teste/elton_assemblies.csv', index=False)
