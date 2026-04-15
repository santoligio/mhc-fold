#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn.cluster import DBSCAN
import plotly.express as px
from pathlib import Path
import numpy as np
import re  # added for explained variance reading

# ============================================================
# Tab20 hex codes dictionary (for reference)
# ============================================================
tab20 = mpl.cm.get_cmap('tab20')
TAB20_HEX = {f'color{i}': mpl.colors.to_hex(tab20(i)) for i in range(20)}

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
    "CD chimera": "#ff7f0e",
    "ULBP": "#c7c7c7",
    "UL18": "#7f7f7f",
    "CD SC": "#aec7e8",
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
    color_map = {}
    fallback_colors = list(TAB20_HEX.values())
    for i, cat in enumerate(sorted(categories)):
        if cat in manual_mapping:
            color_map[cat] = manual_mapping[cat]
        else:
            color_map[cat] = fallback_colors[i % 20]
    return color_map

# ============================================================
# Paths
# ============================================================
BASE_DIR = Path('/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/pca/reviewed+length')

PCA_CSV = BASE_DIR / 'combined_pca_result_mhc.csv'
ANNOTATIONS_CSV = Path(
    '/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/'
    'functional_annotation/filtered/'
    'afdb_mhc_annotations_filtered.csv'
)
REVIEWED_MODELS_CSV = Path(
    '/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/'
    'functional_annotation/new_filters/'
    'afdb_models_analysis_reviewed.csv'
)

PLOTS_DIR = BASE_DIR / 'plots_union'
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

PLOT_GENE_OUTPUT = PLOTS_DIR / 'mhc_pca_plot.png'
PLOT_CLASS_OUTPUT = PLOTS_DIR / 'mhc_pca_plot_class.png'
PLOT_SUPERCLASS_OUTPUT = PLOTS_DIR / 'mhc_pca_plot_superclass.png'

PLOTLY_GENE_HTML = PLOTS_DIR / 'mhc_pca_plot_interactive.html'
PLOTLY_CLASS_HTML = PLOTS_DIR / 'mhc_pca_plot_class_interactive.html'
PLOTLY_SUPERCLASS_HTML = PLOTS_DIR / 'mhc_pca_plot_superclass_interactive.html'

LABELS_SIZE = 16
TICKS_SIZE = 14

# ============================================================
# Load data
# ============================================================
pca_df = pd.read_csv(PCA_CSV)
annotations_df = pd.read_csv(ANNOTATIONS_CSV)
reviewed_df = pd.read_csv(REVIEWED_MODELS_CSV)

# ============================================================
# Extract UniProt from reviewed file
# ============================================================
reviewed_uniprots = set(
    reviewed_df['pdb']
    .astype(str)
    .str.split('-')
    .str[1]
    .str.strip()
    .str.upper()
)

# ============================================================
# Harmonize identifiers
# ============================================================
pca_df['merge_id'] = pca_df['UniProt_ID'].astype(str).str.strip().str.upper()
annotations_df['merge_id'] = annotations_df['uniprot_id'].astype(str).str.strip().str.upper()

# ============================================================
# Merge PCA with annotation
# ============================================================
merged_df = pd.merge(
    pca_df,
    annotations_df[['merge_id','gene_name','mapped_class','mapped_superclass']],
    on='merge_id',
    how='inner'
)

# ============================================================
# Mark reviewed
# ============================================================
merged_df['Reviewed'] = merged_df['merge_id'].isin(reviewed_uniprots)
print("Number of reviewed True:", merged_df['Reviewed'].sum())

# ============================================================
# Clustering
# ============================================================
db = DBSCAN(eps=0.2, min_samples=5)
merged_df['cluster'] = db.fit_predict(merged_df[['PC1', 'PC2']]).astype(str)

# ============================================================
# Read explained variance from log
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
# Plotting helpers
# ============================================================
def plot_static_pca(df, color_col, output_path):
    fig, ax = plt.subplots(figsize=(8, 6))
    categories = sorted(df[color_col].dropna().unique())

    # Choose manual color mapping if available
    if color_col == 'mapped_class':
        colors = get_manual_colors(categories, MANUAL_CLASS_COLORS)
    elif color_col == 'mapped_superclass':
        colors = get_manual_colors(categories, MANUAL_SUPERCLASS_COLORS)
                
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
        cmap = mpl.colormaps['turbo']
        colors = dict(zip(categories, cmap(np.linspace(0, 1, len(categories)))))

    # --- Plot points ---
    for cat in categories:
        subset = df[df[color_col] == cat]

        # Unreviewed circles
        unreviewed_subset = subset[~subset['Reviewed']]
        if not unreviewed_subset.empty:
            ax.scatter(
                unreviewed_subset['PC1'],
                unreviewed_subset['PC2'],
                c=[colors[cat]] * len(unreviewed_subset),
                marker='o',
                s=50,
                edgecolor='k',
                linewidth=0.3,
                alpha=0.8,
                zorder=1
            )

        # Reviewed triangles
        reviewed_subset = subset[subset['Reviewed']]
        if not reviewed_subset.empty:
            ax.scatter(
                reviewed_subset['PC1'],
                reviewed_subset['PC2'],
                c=[colors[cat]] * len(reviewed_subset),
                marker='^',
                s=80,
                edgecolor='k',
                linewidth=0.3,
                alpha=0.8,
                zorder=5
            )

    # --- Axis labels with explained variance ---
    ax.set_xlabel(f'PC1 ({explained_var_percent[0]})', fontsize=LABELS_SIZE)
    ax.set_ylabel(f'PC2 ({explained_var_percent[1]})', fontsize=LABELS_SIZE)

    # --- Grid lines ---
    ax.grid(True, linestyle='--', alpha=0.3)

    # --- Legend for categories ---
    legend_title = (
        "Protein identity" if color_col == "gene_name"
        else "Class" if color_col == "mapped_class"
        else "Superclass"
    )

    handles = [
        plt.Line2D([0],[0], marker='o', color='w', label=cat,
                   markerfacecolor=colors[cat],
                   markeredgecolor='k',
                   markersize=8,
                   markeredgewidth=0.3)
        for cat in categories
    ]

    leg = ax.legend(handles=handles, title=legend_title,
                    fontsize=7, title_fontsize=8,
                    loc='center left', bbox_to_anchor=(1.02, 0.5),
                    frameon=True)
    leg.get_frame().set_edgecolor('white')
    leg.get_frame().set_linewidth(1.5)
    leg.get_frame().set_facecolor('white')
    leg.get_frame().set_alpha(1.0)

    # --- Add n = number of points (top left) ---
    ax.text(
        0.02, 0.95, f"n = {len(df)}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment='top',
        horizontalalignment='left'
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def plot_interactive_pca(df, color_col, output_path):
    categories = sorted(df[color_col].dropna().unique())

    if color_col == 'mapped_class':
        color_map = get_manual_colors(categories, MANUAL_CLASS_COLORS)
        legend_title = "Class"
    elif color_col == 'mapped_superclass':
        color_map = get_manual_colors(categories, MANUAL_SUPERCLASS_COLORS)
        legend_title = "Superclass"
    else:
        turbo = mpl.cm.get_cmap('turbo')
        all_colors = [mpl.colors.to_hex(turbo(i / max(len(categories)-1, 1))) for i in range(len(categories))]
        color_map = dict(zip(categories, all_colors))
        legend_title = "Protein identity"

    fig = px.scatter(
        df,
        x='PC1',
        y='PC2',
        color=color_col,
        color_discrete_map=color_map,
        hover_data=['merge_id', 'cluster'],
        width=900,
        height=700
    )

    fig.update_layout(
        legend_title_text=legend_title,
        legend=dict(
            bgcolor="white",
            bordercolor="white",
            borderwidth=2
        ),
        template='plotly_white',
        xaxis_title=f'PC1 ({explained_var_percent[0]})',
        yaxis_title=f'PC2 ({explained_var_percent[1]})'
    )

    fig.write_html(output_path)

# ============================================================
# Generate plots
# ============================================================
plot_static_pca(merged_df, 'gene_name', PLOT_GENE_OUTPUT)
plot_interactive_pca(merged_df, 'gene_name', PLOTLY_GENE_HTML)

plot_static_pca(merged_df, 'mapped_class', PLOT_CLASS_OUTPUT)
plot_interactive_pca(merged_df, 'mapped_class', PLOTLY_CLASS_HTML)

plot_static_pca(merged_df, 'mapped_superclass', PLOT_SUPERCLASS_OUTPUT)
plot_interactive_pca(merged_df, 'mapped_superclass', PLOTLY_SUPERCLASS_HTML)

print("All plots generated successfully.")
