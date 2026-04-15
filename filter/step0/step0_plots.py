import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Standard figure size
FIGSIZE = (7, 5)  # width, height in inches

# -----------------------------
# Plot 1: Query TM-Score vs log E-value
# -----------------------------
fs_db = { 'pdb': '/mnt/4TB/giovanna/foldseek/version_02/pdb/dbs_pdb_aln',
         'afdb': '/mnt/4TB/giovanna/foldseek/version_02/alphafold/3mre_afdb_aln' }

df1 = pd.read_csv(fs_db['pdb'], delimiter=r'\s+')
df1.columns = ['query','target','fident','alnlen','mismatch','gapopen','qstart','qend','tstart','tend',
               'evalue','bits','alntmscore','qtmscore','ttmscore','lddt','lddtfull','prob']

df1 = df1[df1['evalue'] > 0]  # avoid log(0)
df1['log_evalue'] = np.log(df1['evalue'])

fig, ax = plt.subplots(figsize=FIGSIZE)
scatter = ax.scatter(df1['qtmscore'], df1['log_evalue'],
                     c=df1['prob'], cmap='viridis',
                     s=1.0, alpha=0.3)

ax.set_xlabel('Query TM-score', fontsize=12)
ax.set_ylabel('log E-value', fontsize=12)

# reference lines
ax.axhline(y=np.log(0.1), color='gray', linestyle='--', alpha=0.5)
ax.axhline(y=np.log(0.01), color='blue', linestyle='--')
ax.axhline(y=np.log(0.001), color='gray', linestyle='--', alpha=0.5)
ax.axvline(x=0.5, color='gray', linestyle='--', alpha=0.5)

cbar = fig.colorbar(scatter, ax=ax)
cbar.set_label("Probability")

plt.tight_layout()
plt.savefig('/mnt/4TB/giovanna/foldseek/version_02/filter/step0/pdb/foldseek_pdb_tm.png', dpi=300)
plt.show()

# -----------------------------
# Plot 2: Alignment Length vs log E-value
# -----------------------------
fig, ax = plt.subplots(figsize=FIGSIZE)
scatter = ax.scatter(df1['alnlen'], df1['log_evalue'],
                     c=df1['prob'], cmap='viridis',
                     s=1.0, alpha=0.3)

ax.set_xlabel('Alignment Length', fontsize=12)
ax.set_ylabel('log E-value', fontsize=12)

# reference lines
ax.axhline(y=np.log(0.01), color='blue', linestyle='--')
ax.axvline(x=90, color='green', linestyle='--')

cbar = fig.colorbar(scatter, ax=ax)
cbar.set_label("Probability")

plt.tight_layout()
plt.savefig('/mnt/4TB/giovanna/foldseek/version_02/filter/step0/pdb/foldseek_pdb_alnlen.png', dpi=300)
plt.show()
