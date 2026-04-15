#!/usr/bin/env python3

import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cavities", required=True,
                    help="filter_cavities_out CSV file")
    ap.add_argument("--annotations", required=True,
                    help="pdb_mhc_annotations CSV file")
    args = ap.parse_args()

    # ------------------------------------------------------------
    # Load files
    # ------------------------------------------------------------
    df_cav = pd.read_csv(args.cavities)
    df_ann = pd.read_csv(args.annotations)

    # ------------------------------------------------------------
    # Basic validation
    # ------------------------------------------------------------
    if "structure" not in df_cav.columns or "volume" not in df_cav.columns:
        raise SystemExit("Cavities CSV must contain 'structure' and 'volume' columns")

    if "pdb_id" not in df_ann.columns or "mapped_class" not in df_ann.columns:
        raise SystemExit("Annotations CSV must contain 'pdb_id' and 'mapped_class' columns")

    # ------------------------------------------------------------
    # Extract first 4 letters
    # ------------------------------------------------------------
    df_cav["pdb_id"] = df_cav["structure"].astype(str).str[:4].str.upper()
    df_ann["pdb_id"] = df_ann["pdb_id"].astype(str).str[:4].str.upper()

    # ------------------------------------------------------------
    # Merge mapped_class
    # ------------------------------------------------------------
    df = pd.merge(
        df_cav,
        df_ann[["pdb_id", "mapped_class"]],
        on="pdb_id",
        how="left"
    )

    df = df.dropna(subset=["mapped_class", "volume"])

    # ------------------------------------------------------------
    # Compute median per class
    # ------------------------------------------------------------
    for cls, group in df.groupby("mapped_class"):

        print("--------------------------------------------------")
        print(f"Class: {cls}")
        
        # Compute median
        median_volume = group["volume"].median()
        
        # Compute absolute distance to median
        group = group.copy()
        group["distance_to_median"] = (group["volume"] - median_volume).abs()
    
        # Sort by closeness to median
        closest3 = group.sort_values("distance_to_median").head(3)
        
        print(f"Median volume: {median_volume:.3f}")
        
        for rank, (_, row) in enumerate(closest3.iterrows(), start=1):
            print(
                f"Rank {rank}: {row['structure']}  |  "
                f"Volume: {row['volume']:.3f}  |  "
                f"Distance to median: {row['distance_to_median']:.3f}"
            )

            print("--------------------------------------------------")


if __name__ == "__main__":
    main()
