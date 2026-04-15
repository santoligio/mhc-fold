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
# Manual color mappings for mapped_class
# ============================================================
MANUAL_CLASS_COLORS = {
    "HLA Class Ia": "#9467bd",
    "HLA Class Ib": "#8c564b",
    "HLA Class Ia SC": "#c5b0d5",
    "MR1": "#e377c2",
    "MR1 SC": "#f7b6d2",
    "ZAG": "#bcbd22",
    "HFE": "#ff9896",
    "FcRn": "#d62728",
    "MIC": "#c49c94",
    "CD1": "#ffbb78",
    "EPCR": "#98df8a",
    "CD1 chimera": "#ff7f0e",
    "ULBP": "#c7c7c7",
    "UL18": "#7f7f7f",
    "CD1 SC": "#aec7e8",
    "OMCP": "#2ca02c",
    "2L": "#1f77b4"
}

# ============================================================
# Manual color mappings for mapped_superclass
# ============================================================
MANUAL_SUPERCLASS_COLORS = {
    "Classical MHC-I": "#3182bd",
    "Non-classical MHC-I": "#756bb1",
    "Classical MHC-I SC": "#9ecae1",
    "MHC-I like": "#fd8d3c",
    "MHC-I like SC": "#fdae6b",
    "MHC-I like chimera": "#fdd0a2",
    "MHC mimetics": "#74c476"
}

# ============================================================
# Helper function to get deterministic colors
# ============================================================
def get_manual_colors(categories, manual_mapping):
    """
    Returns a dict mapping each category to a hex color.
    Categories not in manual_mapping will cycle through tab20 as fallback.
    """
    tab20 = mpl.cm.get_cmap('tab20')
    fallback_colors = [mpl.colors.to_hex(tab20(i)) for i in range(20)]
    color_map = {}
    for i, cat in enumerate(sorted(categories)):
        if cat in manual_mapping:
            color_map[cat] = manual_mapping[cat]
        else:
            color_map[cat] = fallback_colors[i % 20]
    return color_map

# ============================================================
# File paths
# ============================================================
BASE_DIR = Path('/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/pca/pdb_vs_afdb')

# Combined PCA (contains BOTH PDB and AFDB)
PCA_CSV = BASE_DIR / 'pdb_vs_afdb_pca_result_mhc.csv'

# Separate annotation files
PDB_ANNOTATIONS = Path(
    '/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/'
    'functional_annotation/filtered/'
    'pdb_mhc_annotations_filtered.csv'
)

AFDB_ANNOTATIONS = Path(
    '/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/'
    'functional_annotation/new_filters/'
    'afdb_mhc_annotations_reviewed.csv'
)

PLOTS_DIR = BASE_DIR / 'plots_combined'
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

LABELS_SIZE = 16
TICKS_SIZE = 14

# ============================================================
# Load PCA
# ============================================================
pca_df = pd.read_csv(PCA_CSV)

pca_df['merge_id'] = (
    pca_df['UniProt_ID']
    .astype(str)
    .str.strip()
    .str.upper()
)

# Determine structure type (robust PDB detection: exactly 4 alphanumeric chars)
pca_df['StructureType'] = np.where(
    pca_df['merge_id'].str.fullmatch(r'[A-Z0-9]{4}'),
    'PDB',
    'AFDB'
)

# ============================================================
# Load and harmonize annotations
# ============================================================
pdb_ann = pd.read_csv(PDB_ANNOTATIONS)
afdb_ann = pd.read_csv(AFDB_ANNOTATIONS)

pdb_ann['merge_id'] = (
    pdb_ann['pdb_id']
    .astype(str)
    .str.strip()
    .str.upper()
)

afdb_ann['merge_id'] = (
    afdb_ann['uniprot_id']
    .astype(str)
    .str.strip()
    .str.upper()
)

# Keep same columns
cols = ['merge_id', 'gene_name', 'mapped_class', 'mapped_superclass']

pdb_ann = pdb_ann[cols]
afdb_ann = afdb_ann[cols]

# Combine annotations
annotations_df = pd.concat([pdb_ann, afdb_ann], ignore_index=True)

# ============================================================
# Merge PCA + annotations
# ============================================================
merged_df = pd.merge(
    pca_df,
    annotations_df,
    on='merge_id',
    how='inner'
)

# ============================================================
# Clustering
# ============================================================
db = DBSCAN(eps=0.2, min_samples=5)
merged_df['cluster'] = db.fit_predict(
    merged_df[['PC1', 'PC2']]
).astype(str)

# ============================================================
# Read explained variance from log if exists
# ============================================================
LOG_FILE = BASE_DIR / 'step9.log'
explained_var = [1.0, 0.0]  # fallback

if LOG_FILE.exists():
    with open(LOG_FILE) as f:
        content = f.read()
        match = re.search(r'Explained variance:\s*\[(.*?)\]', content)
        if match:
            explained_var = [float(x) for x in match.group(1).split()]

explained_var_percent = [f"{x*100:.1f}%" for x in explained_var]

# ============================================================
# Static plot helper
# ============================================================
def plot_static_pca(df, color_col, filename):
    fig, ax = plt.subplots(figsize=(8, 6))

    categories = sorted(df[color_col].dropna().unique())

    # Determine color mapping
    if color_col == 'mapped_class':
        color_map = get_manual_colors(categories, MANUAL_CLASS_COLORS)
        legend_title = "Class"
    elif color_col == 'mapped_superclass':
        color_map = get_manual_colors(categories, MANUAL_SUPERCLASS_COLORS)
        legend_title = "Superclass"

        legend_order = [
            "Classical MHC-I",
            "Classical MHC-I SC",
            "Non-classical MHC-I",
            "MHC-I like",
            "MHC-I like SC",
            "MHC-I like chimera",
            "MHC mimetics"
        ]
        # keep only categories that exist
        categories = [c for c in legend_order if c in categories]
    else:
        cmap = mpl.cm.get_cmap('turbo')
        color_map = dict(zip(
            categories,
            [mpl.colors.to_hex(cmap(i / max(len(categories)-1, 1)))
             for i in range(len(categories))]
        ))
        legend_title = "Protein identity"

    # Plot points
    for cat in categories:
        first = True  # ensures only ONE legend entry per category

        for structure_type, marker in [('PDB', 'o'), ('AFDB', '^')]:
            subset = df[
                (df[color_col] == cat) &
                (df['StructureType'] == structure_type)
            ]

            if not subset.empty:
                ax.scatter(
                    subset['PC1'],
                    subset['PC2'],
                    color=color_map[cat],
                    marker=marker,
                    s=50,
                    alpha=0.8,
                    edgecolor='k',
                    linewidth=0.3,
                    label=cat if first else None
                )
                first = False

    ax.set_xlabel(f'PC1 ({explained_var_percent[0]})', fontsize=LABELS_SIZE)
    ax.set_ylabel(f'PC2 ({explained_var_percent[1]})', fontsize=LABELS_SIZE)
    ax.tick_params(labelsize=TICKS_SIZE)
    ax.grid(linestyle='--', alpha=0.3)

    # Unique legend
    handles = [
        plt.Line2D(
            [0], [0],
            marker='o',
            color='w',
            label=cat,
            markerfacecolor=color_map[cat],
            markeredgecolor='k',
            markersize=8,
            markeredgewidth=0.3
        )
        for cat in categories
    ]

    leg = ax.legend(
        handles=handles,
        title=legend_title,
        fontsize=7,
        title_fontsize=8,
        loc='center left',
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    leg.get_frame().set_edgecolor('white')
    
    # Add n
    ax.text(
        0.02, 0.95, f'n = {len(df)}',
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top'
    )

    plt.tight_layout()
    plt.savefig(PLOTS_DIR / filename, dpi=300, bbox_inches='tight')
    plt.close()
    
# ============================================================
# Interactive plot helper
# ============================================================
def plot_interactive_pca(df, color_col, filename):
    categories = sorted(df[color_col].dropna().unique())

    # Determine color mapping
    if color_col == 'mapped_class':
        color_map = get_manual_colors(categories, MANUAL_CLASS_COLORS)
        legend_title = "Class"
    elif color_col == 'mapped_superclass':
        color_map = get_manual_colors(categories, MANUAL_SUPERCLASS_COLORS)
        legend_title = "Superclass"
    else:  # gene_name or others
        cmap = mpl.cm.get_cmap('turbo')
        color_map = dict(zip(categories, [mpl.colors.to_hex(cmap(i / max(len(categories)-1, 1))) for i in range(len(categories))]))
        legend_title = "Protein identity"

    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color=color_col,
        symbol='StructureType',
        symbol_map={
            'PDB': 'circle',
            'AFDB': 'triangle-up'
        },
        color_discrete_map=color_map,
        hover_data=['merge_id', 'cluster'],
        width=900,
        height=700
    )

    fig.update_traces(marker=dict(size=10, line=dict(width=0.5, color='black')))
    
    fig.update_layout(
        template='plotly_white',
        legend_title_text=legend_title,
        legend=dict(
            bgcolor="white",
            bordercolor="white",
            borderwidth=2,
            traceorder='normal'
        ),
        xaxis_title=f'PC1 ({explained_var_percent[0]})',
        yaxis_title=f'PC2 ({explained_var_percent[1]})'
    )

    # Add n = number of points
    fig.add_annotation(
        text=f"n = {len(df)}",
        xref="paper", yref="paper",
        x=0.02, y=0.95,
        showarrow=False,
        font=dict(size=12)
    )

    fig.write_html(PLOTS_DIR / filename)

# ============================================================
# Generate plots
# ============================================================
plot_static_pca(merged_df, 'gene_name', 'combined_mhc_pca.png')
plot_static_pca(merged_df, 'mapped_class', 'combined_mhc_pca_class.png')
plot_static_pca(merged_df, 'mapped_superclass', 'combined_mhc_pca_superclass.png')

plot_interactive_pca(merged_df, 'gene_name', 'combined_mhc_pca__interactive.html')
plot_interactive_pca(merged_df, 'mapped_class', 'combined_mhc_pca_class_interactive.html')
plot_interactive_pca(merged_df, 'mapped_superclass', 'combined_mhc_pca_superclass_interactive.html')

print("All plots generated successfully.")
