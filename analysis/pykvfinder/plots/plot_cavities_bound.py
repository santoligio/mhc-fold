#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from matplotlib.lines import Line2D


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

    ap.add_argument("--pdb-csv", required=True)
    ap.add_argument("--afdb-csv", required=True)

    ap.add_argument("--pdb-annotations", required=True)
    ap.add_argument("--afdb-annotations", required=True)

    ap.add_argument("--contacts", required=True,
                    help="mhc_contacts.csv file")

    ap.add_argument("--output", default="cavity_bound_plot.png")
    ap.add_argument("--interactive-output",
                    default="cavity_bound_plot.html")

    args = ap.parse_args()

    # ------------------------------------------------------------
    # Load cavity CSVs
    # ------------------------------------------------------------

    df_pdb = pd.read_csv(args.pdb_csv)
    df_af = pd.read_csv(args.afdb_csv)
    df_contacts = pd.read_csv(args.contacts)

    required_cols = {"structure", "cavity_id", "volume", "avg_depth"}
    if not required_cols.issubset(df_pdb.columns):
        raise SystemExit("PDB CSV missing required columns")

    if not required_cols.issubset(df_af.columns):
        raise SystemExit("AFDB CSV missing required columns")

    # ------------------------------------------------------------
    # Aggregate cavities
    # ------------------------------------------------------------

    df_pdb = df_pdb.groupby("structure", as_index=False).agg({
        "volume": "sum",
        "avg_depth": "mean"
    })

    df_af = df_af.groupby("structure", as_index=False).agg({
        "volume": "sum",
        "avg_depth": "mean"
    })

    # ------------------------------------------------------------
    # Extract IDs
    # ------------------------------------------------------------

    df_pdb["pdb_id"] = df_pdb["structure"].apply(extract_pdb_id)
    df_af["uniprot_id"] = df_af["structure"].apply(extract_af_uniprot)

    df_contacts["pdb_id"] = df_contacts["pdb_id"].str.upper()

    # ------------------------------------------------------------
    # Identify bound PDB structures (ligand entry required)
    # ------------------------------------------------------------

    bound_pdb_ids = set(
        df_contacts[df_contacts["type"].str.lower() == "ligand"]["pdb_id"]
    )

    # ------------------------------------------------------------
    # Assign binding status
    # ------------------------------------------------------------

    df_pdb["binding_status"] = df_pdb["pdb_id"].apply(
        lambda x: "bound" if x in bound_pdb_ids else "unbound"
    )

    df_af["binding_status"] = "unbound"

    df_pdb["source"] = "PDB"
    df_af["source"] = "AFDB"

    # Combine
    df_all = pd.concat([df_pdb, df_af], ignore_index=True)

    # ------------------------------------------------------------
    # Color mapping
    # ------------------------------------------------------------

    color_map = {
        "bound": "#d62728",
        "unbound": "#1f77b4"
    }

    # ------------------------------------------------------------
    # Static Plot (Matplotlib)
    # ------------------------------------------------------------

    plt.figure(figsize=(8, 6))

    for status in ["bound", "unbound"]:
        df_status = df_all[df_all["binding_status"] == status]

        df_pdb_s = df_status[df_status["source"] == "PDB"]
        plt.scatter(
            df_pdb_s["volume"],
            df_pdb_s["avg_depth"],
            color=color_map[status],
            marker="o",
            edgecolor="black",
            linewidth=0.4,
            alpha=0.8
        )

        df_af_s = df_status[df_status["source"] == "AFDB"]
        plt.scatter(
            df_af_s["volume"],
            df_af_s["avg_depth"],
            color=color_map[status],
            marker="^",
            edgecolor="black",
            linewidth=0.4,
            alpha=0.8
        )

    plt.xlabel("Volume (summed cavities)")
    plt.ylabel("Average depth (summed cavities)")

    # Legend (Bound / Unbound only)
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

    plt.legend(handles=legend_elements, frameon=True)

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.close()

    print(f"Plot saved to {args.output}")

    # ------------------------------------------------------------
    # Interactive Plot (Plotly)
    # ------------------------------------------------------------

    df_all["symbol"] = df_all["source"].map({
        "PDB": "circle",
        "AFDB": "triangle-up"
    })

    fig = go.Figure()

    for status in ["bound", "unbound"]:
        df_status = df_all[df_all["binding_status"] == status]

        fig.add_trace(
            go.Scatter(
                x=df_status["volume"],
                y=df_status["avg_depth"],
                mode="markers",
                marker=dict(
                    color=color_map[status],
                    size=9,
                    symbol=df_status["symbol"],
                    line=dict(width=1, color="black")
                ),
                name=status.capitalize(),
                legendgroup=status,
                showlegend=True,
                hovertemplate=(
                    "<b>Status:</b> %{text}<br>"
                    "<b>Source:</b> %{customdata[0]}<br>"
                    "<b>Structure:</b> %{customdata[1]}<br>"
                    "<b>Volume:</b> %{x}<br>"
                    "<b>Avg Depth:</b> %{y}<extra></extra>"
                ),
                text=[status] * len(df_status),
                customdata=df_status[["source", "structure"]].values
            )
        )

    fig.update_layout(
        xaxis_title="Volume (summed cavities)",
        yaxis_title="Average depth (summed cavities)",
        template="simple_white",
        width=900,
        height=600,
        legend=dict(
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
