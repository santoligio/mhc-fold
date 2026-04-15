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
# Tab20 hex codes dictionary (for reference)
# ============================================================
tab20 = mpl.cm.get_cmap('tab20')
TAB20_HEX = {f'color{i}': mpl.colors.to_hex(tab20(i)) for i in range(20)}

# ============================================================
# Manual color mappings for mapped_class
# ============================================================

MANUAL_CLASS_COLORS = {
    "HLA Class Ia": "#9467bd",       # tab20 color0
    "HLA Class Ib": "#8c564b",       # tab20 color1
    "HLA Class Ia SC": "#c5b0d5",    # tab20 color2
    "MR1": "#e377c2",                # tab20 color3
    "MR1 SC": "#f7b6d2",             # tab20 color4
    "ZAG": "#bcbd22",                # tab20 color5
    "HFE": "#ff9896",                # tab20 color6
    "FcRn": "#d62728",               # tab20 color7
    "MIC": "#c49c94",                # tab20 color8
    "CD1": "#ffbb78",                # tab20 color9
    "EPCR": "#98df8a",               # tab20 color10
    "CD chimera": "#ff7f0e",         # tab20 color11
    "ULBP": "#c7c7c7",               # tab20 color12
    "UL18": "#7f7f7f",               # tab20 color13
    "CD SC": "#aec7e8",              # tab20 color14
    "OMCP": "#2ca02c",               # tab20 color15
    "2L": "#1f77b4"                  # tab20 color16
}

# ============================================================
# Manual color mappings for mapped_superclass
# ============================================================

MANUAL_SUPERCLASS_COLORS = {
    "Classical MHC-I": "#3182bd",          # tab20 color0
    "Non-classical MHC-I": "#756bb1",      # tab20 color1
    "Classical MHC-I SC": "#9ecae1",       # tab20 color2
    "MHC-I like": "#fd8d3c",               # tab20 color3
    "MHC-I like SC": "#fdae6b",            # tab20 color4
    "MHC-I like chimera": "#fdd0a2",       # tab20 color5
    "MHC mimetics": "#74c476"              # tab20 color6
}

# ============================================================
# Helper function to get deterministic colors
# ============================================================
def get_manual_colors(categories, manual_mapping):
    """
    Returns a dict mapping each category to a hex color.
    Categories not in manual_mapping will cycle through TAB20_HEX as fallback.
    """
    color_map = {}
    fallback_colors = list(TAB20_HEX.values())
    for i, cat in enumerate(sorted(categories)):
        if cat in manual_mapping:
            color_map[cat] = manual_mapping[cat]
        else:
            color_map[cat] = fallback_colors[i % 20]
    return color_map

# ============================================================
# Configuration
# ============================================================
# Mode can be: "PDB" or "AFDB"
MODE = "PDB"   # change to "PDB" when needed

# ============================================================
# File paths
# ============================================================
BASE_DIR = Path('/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/pca/pdb')

PCA_CSV = BASE_DIR / (
    'afdb_pca_result_mhc.csv' if MODE == "AFDB" else 'pdb_pca_result_mhc.csv'
)

ANNOTATIONS_CSV = Path(
    '/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/'
    'functional_annotation/filtered/'
    f'{"afdb" if MODE == "AFDB" else "pdb"}_mhc_annotations_filtered.csv'
)

PLOTS_DIR = BASE_DIR / 'plots'
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
# Generic PCA plotting helpers
# ============================================================

def plot_static_pca(df, color_col, output_path):
    fig, ax = plt.subplots(figsize=(8, 6))

    categories = sorted(df[color_col].dropna().unique())

    # --- Determine which color mapping to use ---
    if color_col == 'gene_name':
        turbo = mpl.cm.get_cmap('turbo')
        # Generate n colors evenly spaced over the colormap
        all_colors = [mpl.colors.to_hex(turbo(i / max(len(categories)-1, 1))) for i in range(len(categories))]
        color_map = dict(zip(categories, all_colors))

        # Legend order: purely alphabetical
        legend_order = sorted(categories)
        legend_title = "Protein identity"

    elif color_col == 'mapped_class':
        color_map = get_manual_colors(categories, MANUAL_CLASS_COLORS)

        # Legend order: purely alphabetical
        legend_order = sorted(categories)
        legend_title = "Class"

    elif color_col == 'mapped_superclass':
        color_map = get_manual_colors(categories, MANUAL_SUPERCLASS_COLORS)
        # Specific order for superclass
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
        legend_order = [c for c in legend_order if c in categories]
        legend_title = "Superclass"

    else:
        cmap = plt.get_cmap('tab20')
        fallback_colors = cmap(np.linspace(0, 1, len(categories)))
        color_map = dict(zip(categories, fallback_colors))
        legend_order = categories
        legend_title = color_col

    # --- Plot each category ---
    for cat in categories:
        subset = df[df[color_col] == cat]
        ax.scatter(
            subset['PC1'],
            subset['PC2'],
            label=cat,
            color=color_map[cat],
            alpha=0.8,
            s=50,
            edgecolor='k',
            linewidth=0.3
        )

    ax.set_xlabel(f'PC1 ({explained_var_percent[0]})', fontsize=LABELS_SIZE)
    ax.set_ylabel(f'PC2 ({explained_var_percent[1]})', fontsize=LABELS_SIZE)
    ax.tick_params(labelsize=TICKS_SIZE)
    ax.grid(linestyle='--', alpha=0.3)

    # --- Legend with white box ---
    handles, labels = ax.get_legend_handles_labels()

    # Optional filtering
    filtered = [(h, l) for h, l in zip(handles, labels) if 'reviewed' not in l.lower()]
    
    if filtered:
        h, l = zip(*filtered)  # keep order, don't convert to dict
        # reorder according to legend_order, keep only existing
        legend_order_filtered = [c for c in legend_order if c in l]
        h = [h[l.index(c)] for c in legend_order_filtered]
        l = legend_order_filtered
        leg = ax.legend(
            h, l,
            title=legend_title,
            fontsize=7,
            title_fontsize=8,
            loc='center left',
            bbox_to_anchor=(1.02, 0.5),
            frameon=True
        )
        
        # White outline
        leg.get_frame().set_edgecolor('white')
        leg.get_frame().set_linewidth(1.5)
        leg.get_frame().set_facecolor('white')
        leg.get_frame().set_alpha(1.0)    
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
    categories = sorted(df[color_col].dropna().unique())

    # --- Determine which color mapping to use ---
    if color_col == 'gene_name':
        # combine tab20, tab20b, tab20c
        tab20 = mpl.cm.get_cmap('tab20')
        tab20b = mpl.cm.get_cmap('tab20b')
        tab20c = mpl.cm.get_cmap('tab20c')
        all_colors = (
            [mpl.colors.to_hex(tab20(i)) for i in range(20)] +
            [mpl.colors.to_hex(tab20b(i)) for i in range(20)] +
            [mpl.colors.to_hex(tab20c(i)) for i in range(20)]
        )
        color_map = dict(zip(categories, all_colors[:len(categories)]))
        legend_order = sorted(categories)  # full alphabetical
        legend_title = "Protein identity"

    elif color_col == 'mapped_class':
        color_map = get_manual_colors(categories, MANUAL_CLASS_COLORS)
        legend_order = sorted(categories)  # full alphabetical
        legend_title = "Class"

    elif color_col == 'mapped_superclass':
        color_map = get_manual_colors(categories, MANUAL_SUPERCLASS_COLORS)
        # Manual order
        manual_order = [
            "Classical MHC-I",
            "Classical MHC-I SC",
            "Non-classical MHC-I",
            "MHC-I like",
            "MHC-I like SC"
        ]
        legend_order = [c for c in manual_order if c in categories]
        legend_title = "Superclass"

    else:
        color_map = {c: None for c in categories}
        legend_order = categories
        legend_title = color_col

    # --- Interactive plot ---
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

    # --- Update layout with white legend box ---
    fig.update_layout(
        legend_title_text=legend_title,
        legend=dict(
            bgcolor="white",
            bordercolor="white",
            borderwidth=2,
            traceorder='normal'  # keeps order in legend as given by legend_order
        ),
        template='plotly_white',
        xaxis_title=f'PC1 ({explained_var_percent[0]})',
        yaxis_title=f'PC2 ({explained_var_percent[1]})',
    )

    # --- Add n as annotation ---
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
)

plot_interactive_pca(
    merged_df,
    color_col='gene_name',
    output_path=PLOTLY_GENE_HTML,
)

# ============================================================
# PCA plots — mapped_class
# ============================================================
plot_static_pca(
    merged_df,
    color_col='mapped_class',
    output_path=PLOT_CLASS_OUTPUT,
)

plot_interactive_pca(
    merged_df,
    color_col='mapped_class',
    output_path=PLOTLY_CLASS_HTML,
)

# ============================================================
# PCA plots — mapped_superclass
# ============================================================
plot_static_pca(
    merged_df,
    color_col='mapped_superclass',
    output_path=PLOT_SUPERCLASS_OUTPUT,
)

plot_interactive_pca(
    merged_df,
    color_col='mapped_superclass',
    output_path=PLOTLY_SUPERCLASS_HTML,
)

print(f"[{MODE}] Static gene plot saved to: {PLOT_GENE_OUTPUT}")
print(f"[{MODE}] Interactive gene plot saved to: {PLOTLY_GENE_HTML}")
print(f"[{MODE}] Static class plot saved to: {PLOT_CLASS_OUTPUT}")
print(f"[{MODE}] Interactive class plot saved to: {PLOTLY_CLASS_HTML}")
print(f"[{MODE}] Static superclass plot saved to: {PLOT_SUPERCLASS_OUTPUT}")
print(f"[{MODE}] Interactive superclass plot saved to: {PLOTLY_SUPERCLASS_HTML}")
