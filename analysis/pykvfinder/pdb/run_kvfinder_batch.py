#!/usr/bin/env python3
"""
Batch-run pyKVFinder on aligned PDBs in a directory (no box required).

Supports filtering by:
  - AFDB (UniProt ID)
  - PDB (4-letter PDB ID)

Use --filter-type to choose.
"""

import os
os.environ["KMP_WARNINGS"] = "0"

import argparse
from pathlib import Path
import traceback
import pyKVFinder
import pandas as pd


# ---------------------------------------------------------
# ID extraction helpers
# ---------------------------------------------------------

def extract_uniprot_from_filename(filename: str):
    """
    Example AFDB filename:
    AF-Q9Y2R2-F1_model_v4_aligned.pdb

    Extracts: Q9Y2R2
    """
    try:
        return filename.split('_')[0].split('-')[1]
    except Exception:
        return None


def extract_pdbid_from_filename(filename: str):
    """
    Example PDB filename:
    1abc_aligned.pdb

    Extracts: 1ABC
    """
    try:
        return filename[:4].upper()
    except Exception:
        return None


# ---------------------------------------------------------
# KVFinder execution
# ---------------------------------------------------------

def run_one(
    pdb_path: Path,
    out_dir: Path,
    probe_in: float,
    probe_out: float,
    volume_cutoff: float,
    removal_distance: float,
    include_depth: bool,
    include_hydropathy: bool,
    hydrophobicity_scale: str,
    overwrite: bool,
):
    stem = pdb_path.stem
    this_out = out_dir / stem
    this_out.mkdir(parents=True, exist_ok=True)

    results_toml = this_out / "results.toml"
    cavity_pdb = this_out / "cavity.pdb"

    if not overwrite and results_toml.exists() and cavity_pdb.exists():
        return "skip", str(results_toml), str(cavity_pdb)

    results = pyKVFinder.run_workflow(
        str(pdb_path),
        include_depth=include_depth,
        include_hydropathy=include_hydropathy,
        hydrophobicity_scale=hydrophobicity_scale,
        probe_in=probe_in,
        probe_out=probe_out,
        volume_cutoff=volume_cutoff,
        removal_distance=removal_distance,
    )

    results.export_all(
        fn=str(results_toml),
        output=str(cavity_pdb),
        include_frequencies_pdf=False,
    )

    return "ok", str(results_toml), str(cavity_pdb)


# ---------------------------------------------------------
# Main
# ---------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in-dir", required=True)
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--pattern", default="*_aligned.pdb")

    ap.add_argument("--filter-csv", required=True,
                    help="CSV file containing column 'uniprot_id' or 'pdb_id'")
    ap.add_argument("--filter-type", required=True,
                    choices=["afdb", "pdb"],
                    help="Filtering mode")

    ap.add_argument("--probe-in", type=float, default=1.4)
    ap.add_argument("--probe-out", type=float, default=10.0)
    ap.add_argument("--volume-cutoff", type=float, default=1.0)
    ap.add_argument("--removal-distance", type=float, default=1.8)

    ap.add_argument("--include-depth", action="store_true", default=True)
    ap.add_argument("--no-depth", action="store_false", dest="include_depth")

    ap.add_argument("--include-hydropathy", action="store_true", default=True)
    ap.add_argument("--no-hydropathy", action="store_false", dest="include_hydropathy")

    ap.add_argument("--hydrophobicity-scale", default="EisenbergWeiss")

    ap.add_argument("--overwrite", action="store_true")
    ap.add_argument("--fail-fast", action="store_true")
    args = ap.parse_args()

    in_dir = Path(args.in_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    filter_csv = Path(args.filter_csv).resolve()

    if not in_dir.is_dir():
        raise SystemExit(f"ERROR: --in-dir does not exist: {in_dir}")

    if not filter_csv.exists():
        raise SystemExit(f"ERROR: --filter-csv not found: {filter_csv}")

    out_dir.mkdir(parents=True, exist_ok=True)

    # -----------------------------------------------------
    # Load filter list
    # -----------------------------------------------------
    df_filter = pd.read_csv(filter_csv)

    if args.filter_type == "afdb":
        if "uniprot_id" not in df_filter.columns:
            raise SystemExit("ERROR: 'uniprot_id' column not found in CSV")
        allowed_ids = set(df_filter["uniprot_id"].astype(str))

    else:  # pdb
        if "pdb" not in df_filter.columns:
            raise SystemExit("ERROR: 'pdb' column not found in CSV")

        # Normalize CSV values → first 4 chars uppercase
        allowed_ids = set(
            df_filter["pdb"]
            .astype(str)
            .str[:4]
            .str.upper()
        )

    pdbs = sorted(in_dir.glob(args.pattern))
    if not pdbs:
        raise SystemExit(f"No files matched pattern '{args.pattern}' in {in_dir}")

    ok = skip = err = filtered_out = 0
    errors = []

    for pdb_path in pdbs:
        try:
            if args.filter_type == "afdb":
                file_id = extract_uniprot_from_filename(pdb_path.name)
            else:
                file_id = extract_pdbid_from_filename(pdb_path.name)

            if not file_id or file_id not in allowed_ids:
                filtered_out += 1
                print(f"[filtered] {pdb_path.name}")
                continue

            status, results_toml, cavity_pdb = run_one(
                pdb_path=pdb_path,
                out_dir=out_dir,
                probe_in=args.probe_in,
                probe_out=args.probe_out,
                volume_cutoff=args.volume_cutoff,
                removal_distance=args.removal_distance,
                include_depth=args.include_depth,
                include_hydropathy=args.include_hydropathy,
                hydrophobicity_scale=args.hydrophobicity_scale,
                overwrite=args.overwrite,
            )

            if status == "ok":
                ok += 1
                print(f"[ok]   {pdb_path.name}")
            else:
                skip += 1
                print(f"[skip] {pdb_path.name}")

        except Exception as e:
            err += 1
            msg = f"[err]  {pdb_path.name}: {e}"
            print(msg)
            tb = traceback.format_exc()
            errors.append((str(pdb_path), msg, tb))
            if args.fail_fast:
                print(tb)
                raise

    print("\nSummary:")
    print(f"  ok:        {ok}")
    print(f"  skip:      {skip}")
    print(f"  filtered:  {filtered_out}")
    print(f"  err:       {err}")

    if errors:
        log_path = out_dir / "kvfinder_errors.log"
        with open(log_path, "w") as f:
            for pdb, msg, tb in errors:
                f.write(pdb + "\n")
                f.write(msg + "\n")
                f.write(tb + "\n")
                f.write("-" * 80 + "\n")
        print(f"\nError log written to: {log_path}")

if __name__ == "__main__":
    main()
