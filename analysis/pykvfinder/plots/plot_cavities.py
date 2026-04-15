#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from pathlib import Path


# ------------------------------------------------------------
# Manual color mapping
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# ID extraction
# ------------------------------------------------------------

def extract_af_uniprot(structure_name):
    try:
        return structure_name.split("_")[0].split("-")[1]
    except Exception:
        return None


def extract_pdb_id(structure_name):
    return structure_name[:4].upper()


# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()

    ap.add_argument("--pdb-csv", required=True,
                    help="Cavities CSV for PDB structures")
    ap.add_argument("--afdb-csv", required=True,
                    help="Cavities CSV for AFDB structures")

    ap.add_argument("--pdb-annotations", required=True,
                    help="PDB annotations CSV (must contain pdb_id, mapped_class)")
    ap.add_argument("--afdb-annotations", required=True,
                    help="AFDB annotations CSV (must contain uniprot_id, mapped_class)")

    ap.add_argument("--output", default="cavity_2D_plot.png")

    ap.add_argument("--interactive-output",
                    default="cavity_2D_plot.html",
                    help="Output HTML for interactive Plotly version")

    args = ap.parse_args()

    # ------------------------------------------------------------
    # Load cavity CSVs
    # ------------------------------------------------------------

    df_pdb = pd.read_csv(args.pdb_csv)
    df_af = pd.read_csv(args.afdb_csv)

    required_cols = {"structure", "cavity_id", "volume", "avg_depth"}
    if not required_cols.issubset(df_pdb.columns):
        raise SystemExit("PDB CSV missing required columns")

    if not required_cols.issubset(df_af.columns):
        raise SystemExit("AFDB CSV missing required columns")

    # ------------------------------------------------------------
    # Aggregate cavities (sum per structure)
    # ------------------------------------------------------------

    df_pdb = (
        df_pdb.groupby("structure", as_index=False)
        .agg({
            "volume": "sum",
            "avg_depth": "mean"
        })
    )

    df_af = (
        df_af.groupby("structure", as_index=False)
        .agg({
            "volume": "sum",
            "avg_depth": "mean"
        })
    )

    # ------------------------------------------------------------
    # Extract IDs
    # ------------------------------------------------------------

    df_pdb["pdb_id"] = df_pdb["structure"].apply(extract_pdb_id)
    df_af["uniprot_id"] = df_af["structure"].apply(extract_af_uniprot)

    # ------------------------------------------------------------
    # Load annotations
    # ------------------------------------------------------------

    df_pdb_ann = pd.read_csv(args.pdb_annotations)
    df_af_ann = pd.read_csv(args.afdb_annotations)

    # Normalize PDB annotation IDs (first 4 letters only)
    df_pdb_ann["pdb_id"] = df_pdb_ann["pdb_id"].str[:4].str.upper()

    # ------------------------------------------------------------
    # Merge mapped_class
    # ------------------------------------------------------------

    df_pdb = df_pdb.merge(
        df_pdb_ann[["pdb_id", "mapped_class"]],
        on="pdb_id",
        how="left"
    )

    df_af = df_af.merge(
        df_af_ann[["uniprot_id", "mapped_class"]],
        on="uniprot_id",
        how="left"
    )

    # Remove entries without class
    df_pdb = df_pdb.dropna(subset=["mapped_class"])
    df_af = df_af.dropna(subset=["mapped_class"])

    # Add source flag
    df_pdb["source"] = "PDB"
    df_af["source"] = "AFDB"

    # Combine
    df_all = pd.concat([df_pdb, df_af], ignore_index=True)

    # ------------------------------------------------------------
    # Color mapping (manual)
    # ------------------------------------------------------------

    classes = sorted(df_all["mapped_class"].unique())

    class_to_color = {
        cls: MANUAL_CLASS_COLORS.get(cls, "#333333")
        for cls in classes
    }

    # ------------------------------------------------------------
    # Static Plot (Matplotlib)
    # ------------------------------------------------------------

    plt.figure(figsize=(8, 6))

    for cls in classes:
        df_cls = df_all[df_all["mapped_class"] == cls]

        # PDB circles
        df_pdb_cls = df_cls[df_cls["source"] == "PDB"]
        plt.scatter(
            df_pdb_cls["volume"],
            df_pdb_cls["avg_depth"],
            color=class_to_color[cls],
            marker="o",
            edgecolor="black",
            linewidth=0.4,
            label=cls,
            alpha=0.8
        )

        # AFDB triangles
        df_af_cls = df_cls[df_cls["source"] == "AFDB"]
        plt.scatter(
            df_af_cls["volume"],
            df_af_cls["avg_depth"],
            color=class_to_color[cls],
            marker="^",
            edgecolor="black",
            linewidth=0.5,
            alpha=0.8
        )

    plt.xlabel("Volume (summed cavities)")
    plt.ylabel("Average depth (summed cavities)")

    # Remove duplicate legend entries
    handles, labels = plt.gca().get_legend_handles_labels()
    unique = dict(zip(labels, handles))

    legend = plt.legend(
        unique.values(),
        unique.keys(),
        title="Class",
        loc="upper left",
        bbox_to_anchor=(1.02, 1),
        borderaxespad=0,
        frameon=True
    )

    legend.get_frame().set_facecolor("white")   # background
    legend.get_frame().set_edgecolor("white")   # 👈 white outline
    legend.get_frame().set_linewidth(1.5)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.show()

    print(f"Plot saved to {args.output}")

    # ------------------------------------------------------------
    # Interactive Plot (Plotly)
    # ------------------------------------------------------------

    df_all["symbol"] = df_all["source"].map({
        "PDB": "circle",
        "AFDB": "triangle-up"
    })

    fig = go.Figure()

    for cls in classes:
        df_cls = df_all[df_all["mapped_class"] == cls]

        fig.add_trace(
            go.Scatter(
                x=df_cls["volume"],
                y=df_cls["avg_depth"],
                mode="markers",
                marker=dict(
                    color=class_to_color[cls],
                    size=9,
                    symbol=df_cls["symbol"],
                    line=dict(width=1, color="black")
                ),
                name=cls,
                legendgroup=cls,
                showlegend=True,
                hovertemplate=(
                    "<b>Class:</b> %{text}<br>"
                    "<b>Source:</b> %{customdata[0]}<br>"
                    "<b>Structure:</b> %{customdata[1]}<br>"
                    "<b>Volume:</b> %{x}<br>"
                    "<b>Avg Depth:</b> %{y}<extra></extra>"
                ),
                text=[cls] * len(df_cls),
                customdata=df_cls[["source", "structure"]].values
            )
        )

    fig.update_layout(
        xaxis_title="Volume (summed cavities)",
        yaxis_title="Average depth (summed cavities)",
        template="simple_white",
        width=900,
        height=600,
        legend=dict(
            title="Class",
            x=1.02,
            y=1,
            xanchor="left",
            yanchor="top",
            bordercolor="black",
            borderwidth=1
        ),
        margin=dict(r=200)
    )

    fig.write_html(args.interactive_output)
    print(f"Interactive plot saved to {args.interactive_output}")


if __name__ == "__main__":
    main()
