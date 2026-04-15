#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import matplotlib.colors as mcolors

# ------------------------------------------------------------
# Input files
# ------------------------------------------------------------
CSV_PDB_MAIN   = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/pdb_info_helices.csv"
CSV_AFDB_MAIN  = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/afdb_info_helices.csv"

CSV_CONTACTS   = "/mnt/c/Users/gio/Documents/foldseek/version_02/filter/step5/pdb/mhc_contacts.csv"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
df_pdb_main  = pd.read_csv(CSV_PDB_MAIN)
df_afdb_main = pd.read_csv(CSV_AFDB_MAIN)
df_contacts  = pd.read_csv(CSV_CONTACTS)

# ------------------------------------------------------------
# ---- PDB processing ----
# ------------------------------------------------------------
df_pdb_main["pdb_id"] = df_pdb_main["file"].str[:4].str.upper()
df_contacts["pdb_id"] = df_contacts["pdb_id"].astype(str).str.upper()

# Identify bound PDB structures (must contain ligand entry)
bound_pdb_ids = set(
    df_contacts[df_contacts["type"].str.lower() == "ligand"]["pdb_id"]
)

df_pdb_main["source"] = "PDB"
df_pdb_main["binding_status"] = df_pdb_main["pdb_id"].apply(
    lambda x: "bound" if x in bound_pdb_ids else "unbound"
)

# ------------------------------------------------------------
# ---- AFDB processing ----
# ------------------------------------------------------------
df_afdb_main["source"] = "AFDB"
df_afdb_main["binding_status"] = "unbound"

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------
df = pd.concat([df_pdb_main, df_afdb_main], ignore_index=True)

df = df.dropna(subset=["helix_pct_total", "helix_com_dist"])

# ------------------------------------------------------------
# Color mapping (Bound vs Unbound)
# ------------------------------------------------------------
color_map = {
    "bound": "#d62728",     # red
    "unbound": "#1f77b4"    # blue
}

# ------------------------------------------------------------
# Plot (Matplotlib)
# ------------------------------------------------------------
from matplotlib.lines import Line2D

plt.figure(figsize=(8, 6))

for status in ["bound", "unbound"]:
    sub = df[df["binding_status"] == status]

    # PDB (circles)
    sub_pdb = sub[sub["source"] == "PDB"]
    plt.scatter(
        sub_pdb["helix_com_dist"],
        sub_pdb["helix_pct_total"],
        marker="o",
        s=40,
        alpha=0.8,
        color=color_map[status],
        edgecolors="black",
        linewidths=0.4
    )

    # AFDB (triangles)
    sub_afdb = sub[sub["source"] == "AFDB"]
    plt.scatter(
        sub_afdb["helix_com_dist"],
        sub_afdb["helix_pct_total"],
        marker="^",
        s=55,
        alpha=0.8,
        color=color_map[status],
        edgecolors="black",
        linewidths=0.4
    )

plt.xlabel("Helix COM Distance (Å)")
plt.ylabel("Helical Residues (%)")

# ------------------------------------------------------------
# Legend (Bound / Unbound only)
# ------------------------------------------------------------
legend_elements = [
    Line2D([0], [0],
           marker='o',
           linestyle='None',
           label='Bound',
           markerfacecolor=color_map["bound"],
           markeredgecolor='black',
           markeredgewidth=0.4,
           markersize=8),

    Line2D([0], [0],
           marker='o',
           linestyle='None',
           label='Unbound',
           markerfacecolor=color_map["unbound"],
           markeredgecolor='black',
           markeredgewidth=0.4,
           markersize=8),
]

plt.legend(
    handles=legend_elements,
    frameon=True
)

plt.tight_layout()

plt.savefig(
    "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/bound_unbound_combined.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# Interactive Plot (Plotly)
# ------------------------------------------------------------
df["symbol"] = df["source"].map({
    "PDB": "circle",
    "AFDB": "triangle-up"
})

fig = go.Figure()

for status in ["bound", "unbound"]:
    sub = df[df["binding_status"] == status]

    fig.add_trace(
        go.Scatter(
            x=sub["helix_com_dist"],
            y=sub["helix_pct_total"],
            mode="markers",
            marker=dict(
                color=color_map[status],
                size=9,
                symbol=sub["symbol"],
                line=dict(width=1, color="black")
            ),
            name=status.capitalize(),
            legendgroup=status,
            showlegend=True,
            customdata=sub[["source", "file"]].values,
            hovertemplate=(
                "<b>Status:</b> %{text}<br>"
                "<b>Source:</b> %{customdata[0]}<br>"
                "<b>File:</b> %{customdata[1]}<br>"
                "<b>Helix COM Distance:</b> %{x:.2f} Å<br>"
                "<b>Helical Residues:</b> %{y:.2f}%<extra></extra>"
            ),
            text=[status] * len(sub)
        )
    )

fig.update_layout(
    xaxis_title="Helix COM Distance (Å)",
    yaxis_title="Helical Residues (%)",
    template="simple_white",
    width=900,
    height=600,
    legend=dict(
        title="Binding Status",
        x=1.02,
        y=1,
        xanchor="left",
        yanchor="top",
        bordercolor="black",
        borderwidth=1
    ),
    margin=dict(r=200)
)

fig.write_html(
    "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/bound_unbound_combined.html"
)

print("Bound vs Unbound plots saved.")
