#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.cluster import DBSCAN
import plotly.express as px
from pathlib import Path
import numpy as np

# ============================================================
# Configuration
# ============================================================
MODE = "AFDB"   # change to "PDB" when needed

# ============================================================
# File paths
# ============================================================
BASE_DIR = Path('/mnt/4TB/giovanna/foldseek/version_02/analysis/pca/new_filter1')

PCA_CSV = BASE_DIR / (
    'afdb_pca_result_mhc.csv' if MODE == "AFDB" else 'pdb_pca_result_mhc.csv'
)

ANNOTATIONS_CSV = Path(
    '/mnt/4TB/giovanna/foldseek/version_02/analysis/'
    'functional_annotation/new_filters/'
    f'{"afdb" if MODE == "AFDB" else "pdb"}_mhc_annotations_length_175_185.csv'
)

# NEW: file containing reviewed structures (ONLY UniProt_ID will be used)
REVIEWED_PCA_CSV = Path(
    '/mnt/4TB/giovanna/foldseek/version_02/analysis/pca/new_filter2/afdb_pca_result_mhc.csv'
)

PLOTS_DIR = BASE_DIR / 'plots2'
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

# Static outputs
PLOT_GENE_OUTPUT = PLOTS_DIR / f'{MODE.lower()}_mhc_pca_plot.png'
PLOT_CLASS_OUTPUT = PLOTS_DIR / f'{MODE.lower()}_mhc_pca_plot_class.png'
PLOT_SUPERCLASS_OUTPUT = PLOTS_DIR / f'{MODE.lower()}_mhc_pca_plot_superclass.png'

# Interactive outputs
PLOTLY_GENE_HTML = PLOTS_DIR / f'{MODE.lower()}_mhc_pca_plot_interactive.html'
PLOTLY_CLASS_HTML = PLOTS_DIR / f'{MODE.lower()}_mhc_pca_plot_class_interactive.html'
PLOTLY_SUPERCLASS_HTML = PLOTS_DIR / f'{MODE.lower()}_mhc_pca_plot_superclass_interactive.html'

LABELS_SIZE = 16
TICKS_SIZE = 14

# ============================================================
# Load data
# ============================================================
pca_df = pd.read_csv(PCA_CSV)
annotations_df = pd.read_csv(ANNOTATIONS_CSV)

# ============================================================
# Load reviewed UniProt IDs (ONLY from new_filter2 PCA file)
# ============================================================
reviewed_ids = set()

if MODE == "AFDB" and REVIEWED_PCA_CSV.exists():
    reviewed_pca_df = pd.read_csv(REVIEWED_PCA_CSV)

    reviewed_ids = set(
        reviewed_pca_df['UniProt_ID']
        .astype(str)
        .str.strip()
        .str.upper()
    )

# ============================================================
# Harmonize identifiers
# ============================================================
if MODE == "PDB":
    pca_df['merge_id'] = (
        pca_df['PDB']
        .astype(str)
        .str.strip()
        .str.upper()
    )

    annotations_df['merge_id'] = (
        annotations_df['pdb_id']
        .astype(str)
        .str.strip()
        .str.upper()
    )

else:  # AFDB
    pca_df['merge_id'] = (
        pca_df['UniProt_ID']
        .astype(str)
        .str.strip()
        .str.upper()
    )

    annotations_df['merge_id'] = (
        annotations_df['uniprot_id']
        .astype(str)
        .str.strip()
        .str.upper()
    )

# ============================================================
# Merge PCA with annotations
# ============================================================
merged_df = pd.merge(
    pca_df,
    annotations_df[
        [
            'merge_id',
            'gene_name',
            'mapped_class',
            'mapped_superclass'
        ]
    ],
    on='merge_id',
    how='inner'
)

# ============================================================
# Mark reviewed structures (based ONLY on new_filter2 PCA IDs)
# ============================================================
merged_df['Reviewed'] = merged_df['merge_id'].isin(reviewed_ids)

print("Number of reviewed True:",
      merged_df['Reviewed'].sum())

# ============================================================
# Clustering (structure-based only)
# ============================================================
db = DBSCAN(eps=0.2, min_samples=5)
merged_df['cluster'] = db.fit_predict(
    merged_df[['PC1', 'PC2']]
).astype(str)

# ============================================================
# Generic PCA plotting helpers
# ============================================================
def plot_static_pca(df, color_col, output_path, title):
    fig, ax = plt.subplots(figsize=(8, 6))

    categories = sorted(df[color_col].dropna().unique())
    cmap = mpl.colormaps['turbo']
    colors = dict(zip(categories,
                      cmap(np.linspace(0, 1, len(categories)))))

    # ---- Plot ALL unreviewed first ----
    for cat in categories:
        subset = df[(df[color_col] == cat) & (~df['Reviewed'])]

        ax.scatter(
            subset['PC1'],
            subset['PC2'],
            label=cat,
            color=colors[cat],
            alpha=0.8,
            s=50,
            edgecolor='k',
            linewidth=0.3,
            zorder=1
        )

    # ---- Plot ALL reviewed on top ----
    reviewed_df = df[df['Reviewed']]

    ax.scatter(
        reviewed_df['PC1'],
        reviewed_df['PC2'],
        c=reviewed_df[color_col].map(colors),
        marker='^',
        s=80,
        edgecolor='black',
        linewidth=0.5,
        zorder=5
    )

    ax.set_xlabel('PC1', fontsize=LABELS_SIZE)
    ax.set_ylabel('PC2', fontsize=LABELS_SIZE)
    ax.tick_params(labelsize=TICKS_SIZE)
    ax.grid(linestyle='--')

    ax.legend(
        title=color_col,
        fontsize=7,
        title_fontsize=8,
        loc='center left',
        bbox_to_anchor=(1.02, 0.5)
    )

    ax.set_title(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_interactive_pca(df, color_col, output_path, title):
    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color=color_col,
        symbol='Reviewed',
        hover_data=['merge_id', 'cluster'],
        title=title,
        width=900,
        height=700
    )

    fig.update_layout(
        legend_title_text=color_col,
        template='plotly_white'
    )

    fig.write_html(output_path)

# ============================================================
# PCA plots
# ============================================================
plot_static_pca(
    merged_df,
    color_col='gene_name',
    output_path=PLOT_GENE_OUTPUT,
    title=f'{MODE} MHC PCA — colored by gene_name'
)

plot_interactive_pca(
    merged_df,
    color_col='gene_name',
    output_path=PLOTLY_GENE_HTML,
    title=f'{MODE} MHC PCA — colored by gene_name'
)

plot_static_pca(
    merged_df,
    color_col='mapped_class',
    output_path=PLOT_CLASS_OUTPUT,
    title=f'{MODE} MHC PCA — colored by mapped_class'
)

plot_interactive_pca(
    merged_df,
    color_col='mapped_class',
    output_path=PLOTLY_CLASS_HTML,
    title=f'{MODE} MHC PCA — colored by mapped_class'
)

plot_static_pca(
    merged_df,
    color_col='mapped_superclass',
    output_path=PLOT_SUPERCLASS_OUTPUT,
    title=f'{MODE} MHC PCA — colored by mapped_superclass'
)

plot_interactive_pca(
    merged_df,
    color_col='mapped_superclass',
    output_path=PLOTLY_SUPERCLASS_HTML,
    title=f'{MODE} MHC PCA — colored by mapped_superclass'
)

print(f"[{MODE}] Static gene plot saved to: {PLOT_GENE_OUTPUT}")
print(f"[{MODE}] Interactive gene plot saved to: {PLOTLY_GENE_HTML}")
print(f"[{MODE}] Static class plot saved to: {PLOT_CLASS_OUTPUT}")
print(f"[{MODE}] Interactive class plot saved to: {PLOTLY_CLASS_HTML}")
print(f"[{MODE}] Static superclass plot saved to: {PLOT_SUPERCLASS_OUTPUT}")
print(f"[{MODE}] Interactive superclass plot saved to: {PLOTLY_SUPERCLASS_HTML}")
