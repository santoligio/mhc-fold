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

CSV_PDB_ANNOT  = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/filtered/pdb_mhc_annotations_filtered.csv"
CSV_AFDB_ANNOT = "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/functional_annotation/new_filters/afdb_mhc_annotations_reviewed.csv"

# ------------------------------------------------------------
# Load data
# ------------------------------------------------------------
df_pdb_main   = pd.read_csv(CSV_PDB_MAIN)
df_afdb_main  = pd.read_csv(CSV_AFDB_MAIN)
df_pdb_annot  = pd.read_csv(CSV_PDB_ANNOT)
df_afdb_annot = pd.read_csv(CSV_AFDB_ANNOT)

# ------------------------------------------------------------
# ---- PDB MERGE ----
# ------------------------------------------------------------
df_pdb_main["pdb_id"] = df_pdb_main["file"].str[:4].str.upper()
df_pdb_annot["pdb_id"] = df_pdb_annot["pdb_id"].astype(str).str.upper()

df_pdb = pd.merge(
    df_pdb_main,
    df_pdb_annot[["pdb_id", "mapped_class"]],
    on="pdb_id",
    how="left"
)

df_pdb["source"] = "PDB"

# ------------------------------------------------------------
# ---- AFDB MERGE ----
# ------------------------------------------------------------
def extract_uniprot(filename):
    prefix = str(filename).split("_")[0]   # AF-P04439
    return prefix.split("-")[1].upper()

df_afdb_main["uniprot_id"] = df_afdb_main["file"].apply(extract_uniprot)
df_afdb_annot["uniprot_id"] = df_afdb_annot["uniprot_id"].astype(str).str.upper()

df_afdb = pd.merge(
    df_afdb_main,
    df_afdb_annot[["uniprot_id", "mapped_class"]],
    on="uniprot_id",
    how="left"
)

df_afdb["source"] = "AFDB"

# ------------------------------------------------------------
# Combine
# ------------------------------------------------------------
df = pd.concat([df_pdb, df_afdb], ignore_index=True)

df = df.dropna(subset=["helix_pct_total", "helix_com_dist", "mapped_class"])

# ------------------------------------------------------------
# Color mapping per class
# ------------------------------------------------------------
classes = sorted(df["mapped_class"].unique())
cmap = plt.get_cmap("tab20")

class_to_color = {
    cls: cmap(i % 20)
    for i, cls in enumerate(classes)
}

# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------
plt.figure(figsize=(8, 6))

for cls in classes:
    sub = df[df["mapped_class"] == cls]

    # PDB (circles)
    sub_pdb = sub[sub["source"] == "PDB"]
    plt.scatter(
        sub_pdb["helix_com_dist"],
        sub_pdb["helix_pct_total"],
        marker="o",
        s=40,
        alpha=0.8,
        color=class_to_color[cls],
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
        color=class_to_color[cls],
        edgecolors="black",
        linewidths=0.4
    )

# ------------------------------------------------------------
# Formatting
# ------------------------------------------------------------
plt.xlabel("Helix COM Distance (Å)")
plt.ylabel("Helical Residues (%)")

# ------------------------------------------------------------
# Clean Legend (Class Colors Only)
# ------------------------------------------------------------

plt.tight_layout()

plt.savefig(
    "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/pdb_afdb_combined.png",
    dpi=300,
    bbox_inches="tight"
)

plt.close()

# ------------------------------------------------------------
# Interactive Plot (Plotly)
# ------------------------------------------------------------

# Convert matplotlib colors to hex for plotly
class_to_color_hex = {
    cls: mcolors.to_hex(class_to_color[cls])
    for cls in classes
}

# Symbol mapping
df["symbol"] = df["source"].map({
    "PDB": "circle",
    "AFDB": "triangle-up"
})

fig = go.Figure()

for cls in classes:
    sub = df[df["mapped_class"] == cls]

    fig.add_trace(
        go.Scatter(
            x=sub["helix_com_dist"],
            y=sub["helix_pct_total"],
            mode="markers",
            marker=dict(
                color=class_to_color_hex[cls],
                size=9,
                symbol=sub["symbol"],
                line=dict(width=1, color="black")
            ),
            name=cls,
            legendgroup=cls,
            showlegend=True,
            customdata=sub[["source", "file"]].values,
            hovertemplate=(
                "<b>Class:</b> %{text}<br>"
                "<b>Source:</b> %{customdata[0]}<br>"
                "<b>File:</b> %{customdata[1]}<br>"
                "<b>Helix COM Distance:</b> %{x:.2f} Å<br>"
                "<b>Helical Residues:</b> %{y:.2f}%<extra></extra>"
            ),
            text=[cls] * len(sub)
        )
    )

fig.update_layout(
    xaxis_title="Helix COM Distance (Å)",
    yaxis_title="Helical Residues (%)",
    template="simple_white",
    width=900,
    height=600,
    legend=dict(
        title="MHC Class",
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
    "/mnt/c/Users/gio/Documents/foldseek/version_02/analysis/helices/pdb_afdb_combined.html"
)

print("Interactive plot saved.")
