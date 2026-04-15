#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import plotly.express as px
from pathlib import Path
import re

# ============================================================
# Paths
# ============================================================

BASE_DIR = Path('/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/binders/pca')

PCA_CSV = BASE_DIR / 'pca_result_binders.csv'

ANNOTATIONS_CSV = Path(
    '/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/binders/'
    'functional_annotation/'
    'binders_annotations_filtered.csv'
)

PLOTS_DIR = BASE_DIR / 'plots'
PLOTS_DIR.mkdir(parents=True, exist_ok=True)

PLOT_CLASS_OUTPUT = PLOTS_DIR / 'binders_pca_class.png'
PLOT_SUPERCLASS_OUTPUT = PLOTS_DIR / 'binders_pca_superclass.png'

PLOTLY_CLASS_HTML = PLOTS_DIR / 'binders_pca_class_interactive.html'
PLOTLY_SUPERCLASS_HTML = PLOTS_DIR / 'binders_pca_superclass_interactive.html'

LABELS_SIZE = 16
TICKS_SIZE = 14

# ============================================================
# Load data
# ============================================================

pca_df = pd.read_csv(PCA_CSV)
annotations_df = pd.read_csv(ANNOTATIONS_CSV)

# ============================================================
# Build merge key
# PCA file: Binder = pdb_chain (e.g. 1AO7_D)
# ============================================================

annotations_df["Binder"] = (
    annotations_df["pdb_id"].astype(str).str.strip().str.upper()
    + "_"
    + annotations_df["new_chain"].astype(str).str.strip()
)

# Ensure PCA binder format is consistent
pca_df["Binder"] = pca_df["Binder"].astype(str).str.strip().str.upper()

# ============================================================
# Merge
# ============================================================

merged_df = pd.merge(
    pca_df,
    annotations_df[["Binder", "mapped_class", "mapped_superclass"]],
    on="Binder",
    how="inner"
)

# ============================================================
# Remove unknown entries
# ============================================================

merged_df = merged_df[
    (merged_df["mapped_class"] != "?") &
    (merged_df["mapped_superclass"] != "?")
].copy()

# ============================================================
# Read explained variance
# ============================================================

LOG_FILE = BASE_DIR / 'run_binders_pca.log'
explained_var = [1.0, 0.0]

if LOG_FILE.exists():
    with open(LOG_FILE) as f:
        content = f.read()
        match = re.search(r'Explained variance:\s*\[(.*?)\]', content)
        if match:
            explained_var = [float(x) for x in match.group(1).split()]

explained_var_percent = [f"{x*100:.1f}%" for x in explained_var]

# ============================================================
# Static PCA plot
# ============================================================

def plot_static_pca(df, color_col, output_path):

    fig, ax = plt.subplots(figsize=(8, 6))
    categories = sorted(df[color_col].unique())

    # CLASS → tab20
    if color_col == "mapped_class":
        legend_title = "Class"
        
        cmap1 = mpl.cm.get_cmap("tab20b").colors
        cmap2 = mpl.cm.get_cmap("tab20").colors
        combined_colors = list(cmap1) + list(cmap2)
        
    else:
        legend_title = "Superclass"

        # Distinct qualitative palette (not tab-based)
        
        cmap1 = mpl.cm.get_cmap("tab20c").colors
        combined_colors = list(cmap1)

    colors = [
        mpl.colors.to_hex(combined_colors[i % len(combined_colors)])
        for i in range(len(categories))
    ]

    color_map = dict(zip(categories, colors))    
    
    for cat in categories:
        subset = df[df[color_col] == cat]
        ax.scatter(
            subset["PC1"],
            subset["PC2"],
            label=cat,
            color=color_map[cat],
            alpha=0.85,
            s=55,
            edgecolor="k",
            linewidth=0.3
        )

    ax.set_xlabel(f"PC1 ({explained_var_percent[0]})", fontsize=LABELS_SIZE)
    ax.set_ylabel(f"PC2 ({explained_var_percent[1]})", fontsize=LABELS_SIZE)
    ax.tick_params(labelsize=TICKS_SIZE)
    ax.grid(linestyle="--", alpha=0.3)

    leg = ax.legend(
        title=legend_title,
        fontsize=8,
        title_fontsize=9,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=True
    )

    leg.get_frame().set_edgecolor("white")
    leg.get_frame().set_facecolor("white")
    leg.get_frame().set_alpha(1.0)

    ax.text(
        0.02, 0.95,
        f"n = {len(df)}",
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top"
    )

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close()

# ============================================================
# Interactive PCA plot
# ============================================================

def plot_interactive_pca(df, color_col, output_path):

    categories = sorted(df[color_col].unique())

    if color_col == "mapped_class":
        legend_title = "Class"
        
        cmap1 = mpl.cm.get_cmap("tab20b").colors
        cmap2 = mpl.cm.get_cmap("tab10").colors
        combined_colors = list(cmap1) + list(cmap2)
        
    else:
        legend_title = "Superclass"

        combined_colors = mpl.cm.get_cmap("Dark2").colors

    colors = [
        mpl.colors.to_hex(combined_colors[i % len(combined_colors)])
        for i in range(len(categories))
    ]

    color_map = dict(zip(categories, colors))
    
    fig = px.scatter(
        df,
        x="PC1",
        y="PC2",
        color=color_col,
        color_discrete_map=color_map,
        hover_data=["Binder"],
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
        template="plotly_white",
        xaxis_title=f"PC1 ({explained_var_percent[0]})",
        yaxis_title=f"PC2 ({explained_var_percent[1]})"
    )

    fig.add_annotation(
        text=f"n = {len(df)}",
        xref="paper", yref="paper",
        x=0.02, y=0.95,
        showarrow=False,
        font=dict(size=12)
    )

    fig.write_html(output_path)

# ============================================================
# Generate plots
# ============================================================

plot_static_pca(merged_df, "mapped_class", PLOT_CLASS_OUTPUT)
plot_interactive_pca(merged_df, "mapped_class", PLOTLY_CLASS_HTML)

plot_static_pca(merged_df, "mapped_superclass", PLOT_SUPERCLASS_OUTPUT)
plot_interactive_pca(merged_df, "mapped_superclass", PLOTLY_SUPERCLASS_HTML)

print("Static class plot saved to:", PLOT_CLASS_OUTPUT)
print("Interactive class plot saved to:", PLOTLY_CLASS_HTML)
print("Static superclass plot saved to:", PLOT_SUPERCLASS_OUTPUT)
print("Interactive superclass plot saved to:", PLOTLY_SUPERCLASS_HTML)
