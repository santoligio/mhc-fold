#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.cluster import DBSCAN
import plotly.express as px
from pathlib import Path
import numpy as np
import re

# ============================================================
# Configuration
# ============================================================
# Mode can be: "PDB" or "AFDB"
MODE = "PDB"   # change to "PDB" when needed

# ============================================================
# File paths
# ============================================================
BASE_DIR = Path('/mnt/4TB/giovanna/foldseek/version_02/analysis/pca/pdb')

PCA_CSV = BASE_DIR / (
    'afdb_pca_result_mhc.csv' if MODE == "AFDB" else 'pdb_pca_result_mhc.csv'
)

ANNOTATIONS_CSV = Path(
    '/mnt/4TB/giovanna/foldseek/version_02/analysis/'
    'functional_annotation/filtered/'
    f'{"afdb" if MODE == "AFDB" else "pdb"}_mhc_annotations_filtered.csv'
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
# Clustering (structure-based only)
# ============================================================
db = DBSCAN(eps=0.2, min_samples=5)
merged_df['cluster'] = db.fit_predict(
    merged_df[['PC1', 'PC2']]
).astype(str)

# ============================================================
# Read explained variance from log
# ============================================================
LOG_FILE = BASE_DIR / 'step10.log'
explained_var = [1.0, 0.0]  # fallback

if LOG_FILE.exists():
    with open(LOG_FILE) as f:
        content = f.read()
        match = re.search(r'Explained variance:\s*\[(.*?)\]', content)
        if match:
            explained_var = [float(x) for x in match.group(1).split()]

explained_var_percent = [f"{x*100:.1f}%" for x in explained_var]

# ============================================================
# Generic PCA plotting helpers
# ============================================================
def plot_static_pca(df, color_col, output_path):
    fig, ax = plt.subplots(figsize=(8, 6))

    categories = sorted(df[color_col].dropna().unique())
    cmap = plt.get_cmap('tab20')
    colors = cmap(np.linspace(0, 1, len(categories)))

    for cat, color in zip(categories, colors):
        subset = df[df[color_col] == cat]

        ax.scatter(
            subset['PC1'],
            subset['PC2'],
            label=cat,
            color=color,
            alpha=0.8,
            s=50,
            edgecolor='k',
            linewidth=0.3
        )

    ax.set_xlabel(f'PC1 ({explained_var_percent[0]})', fontsize=LABELS_SIZE)
    ax.set_ylabel(f'PC2 ({explained_var_percent[1]})', fontsize=LABELS_SIZE)
    ax.tick_params(labelsize=TICKS_SIZE)
    ax.grid(linestyle='--', alpha=0.3)  # lighter grid lines

    # Legend
    handles, labels = ax.get_legend_handles_labels()
    # Remove any unwanted legend entries (like "reviewed/unreviewed")
    filtered = [(h, l) for h, l in zip(handles, labels) if 'reviewed' not in l.lower()]
    if filtered:
        h, l = zip(*filtered)
        ax.legend(
            h,
            l,
            title=color_col,
            fontsize=7,
            title_fontsize=8,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5)
        )

    # Add n = number of points
    ax.text(
        0.02, 0.95, f'n = {len(df)}',
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top'
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_interactive_pca(df, color_col, output_path):
    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color=color_col,
        hover_data=['merge_id', 'cluster'],
        width=900,
        height=700
    )

    fig.update_layout(
        legend_title_text=color_col,
        template='plotly_white',
        xaxis_title=f'PC1 ({explained_var_percent[0]})',
        yaxis_title=f'PC2 ({explained_var_percent[1]})',
    )

    # Add n as annotation
    fig.add_annotation(
        text=f"n = {len(df)}",
        xref="paper", yref="paper",
        x=0.02, y=0.95,
        showarrow=False,
        font=dict(size=12)
    )

    fig.write_html(output_path)
# ============================================================
# PCA plots — gene_name
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

# ============================================================
# PCA plots — mapped_class
# ============================================================
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

# ============================================================
# PCA plots — mapped_superclass
# ============================================================
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
