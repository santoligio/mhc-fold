#!/usr/bin/env python3
"""
Filter pyKVFinder cavities to keep only those near an aligned reference peptide (CA-only PDB),
then write:
  - per-structure filtered cavity pdb
  - a global CSV summary (all structures, all kept cavities)

Assumptions (matching your outputs):
- KVFinder output root contains subfolders, each with:
    results.toml
    cavity.pdb
- cavity.pdb contains cavity points with cavity ID encoded as the residue name (e.g., KAA, KAB, ...)
  (This matches your snippet and also the keys in results.toml.)

Distance rule (your idea):
- Keep cavity ID if ANY cavity point is within < cutoff Å from ANY peptide CA coordinate.

Example:
  python filter_kvfinder_groove_cavities.py \
    --kvfinder-root /path/to/kvfinder_out \
    --peptide-ca-pdb /path/to/ref_peptide_CA_aligned.pdb \
    --cutoff 3.0
"""

import argparse
import csv
from pathlib import Path
from typing import Dict, List, Tuple, Optional

import numpy as np

# TOML reading: Python 3.11+ has tomllib; otherwise use tomli
try:
    import tomllib  # type: ignore
except Exception:  # pragma: no cover
    tomllib = None
    import tomli  # type: ignore


def load_toml(path: Path) -> dict:
    data = path.read_bytes()
    if tomllib is not None:
        return tomllib.loads(data.decode("utf-8"))
    return tomli.loads(data.decode("utf-8"))


def read_peptide_ca_coords(peptide_pdb: Path) -> np.ndarray:
    """
    Read CA coordinates from a PDB (assumed CA-only, but we still filter atom name == 'CA').
    Returns array shape (M, 3).
    """
    coords = []
    with peptide_pdb.open("r") as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            coords.append((x, y, z))
    if not coords:
        raise ValueError(f"No CA atoms found in peptide PDB: {peptide_pdb}")
    return np.array(coords, dtype=float)


def parse_cavity_pdb_points_by_id(cavity_pdb: Path) -> Dict[str, np.ndarray]:
    """
    Parse KVFinder cavity points and group by cavity ID.

    In your cavity.pdb snippet, the cavity ID is the "residue name" field, e.g.:
      ATOM ... HA  KAA  259 ...
                 ^^^
    So we group points by resname (columns 17-20 in PDB format).

    Returns: {cavity_id: coords(N,3)}
    """
    pts: Dict[str, List[Tuple[float, float, float]]] = {}
    with cavity_pdb.open("r") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            # PDB columns (fixed-width):
            # resname: 17-20 (0-based 17:20)
            cav_id = line[17:20].strip()
            if not cav_id:
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            pts.setdefault(cav_id, []).append((x, y, z))

    return {k: np.array(v, dtype=float) for k, v in pts.items()}


def min_dist_points_to_points(A: np.ndarray, B: np.ndarray) -> float:
    """
    Compute min Euclidean distance between two point clouds A (N,3) and B (M,3).
    Uses scipy cKDTree if available; falls back to numpy brute-force.
    """
    try:
        from scipy.spatial import cKDTree  # type: ignore
        tree = cKDTree(B)
        dists, _ = tree.query(A, k=1)
        return float(np.min(dists))
    except Exception:
        # brute force: (N,M,3) can be heavy; but cavity points are usually manageable
        diff = A[:, None, :] - B[None, :, :]
        d2 = np.sum(diff * diff, axis=2)
        return float(np.sqrt(np.min(d2)))


def filter_cavities_near_peptide(
    cavity_points: Dict[str, np.ndarray],
    peptide_ca: np.ndarray,
    cutoff: float,
) -> Tuple[List[str], Dict[str, float]]:
    """
    Returns:
      kept_ids: list of cavity IDs kept
      min_dists: dict cavity_id -> min distance to peptide CA
    """
    kept = []
    min_dists: Dict[str, float] = {}
    for cav_id, pts in cavity_points.items():
        dmin = min_dist_points_to_points(pts, peptide_ca)
        min_dists[cav_id] = dmin
        if dmin < cutoff:
            kept.append(cav_id)
    kept.sort()
    return kept, min_dists


def write_filtered_cavity_pdb(in_pdb: Path, out_pdb: Path, keep_ids: set):
    """
    Write a cavity PDB containing only ATOM lines whose resname matches kept IDs.
    Keep non-ATOM lines (headers) if present.
    """
    with in_pdb.open("r") as fin, out_pdb.open("w") as fout:
        for line in fin:
            if line.startswith("ATOM"):
                cav_id = line[17:20].strip()
                if cav_id in keep_ids:
                    fout.write(line)
            else:
                fout.write(line)


def get_metrics_from_results(results: dict, cav_id: str) -> dict:
    """
    Pull cavity metrics from results.toml structure like:
      RESULTS:
        VOLUME: {KAA: ...}
        AREA: ...
        MAX_DEPTH: ...
        AVG_DEPTH: ...
        AVG_HYDROPATHY: ...
    """
    R = results.get("RESULTS", {})
    def g(section: str) -> Optional[float]:
        return R.get(section, {}).get(cav_id, None)

    return {
        "volume": g("VOLUME"),
        "area": g("AREA"),
        "max_depth": g("MAX_DEPTH"),
        "avg_depth": g("AVG_DEPTH"),
        "avg_hydropathy": g("AVG_HYDROPATHY"),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--kvfinder-root", required=True, help="Root directory produced by your kvfinder_batch.py (contains per-PDB subfolders).")
    ap.add_argument("--peptide-ca-pdb", required=True, help="Aligned reference peptide PDB containing CA atoms.")
    ap.add_argument("--cutoff", type=float, default=3.0, help="Distance cutoff in Å (default: 3.0).")
    ap.add_argument("--results-name", default="results.toml", help="File name inside each subfolder (default: results.toml).")
    ap.add_argument("--cavity-name", default="cavity.pdb", help="File name inside each subfolder (default: cavity.pdb).")
    ap.add_argument("--out-csv", default="groove_cavities_summary.csv", help="CSV filename written inside kvfinder-root.")
    ap.add_argument("--write-filtered-pdb", action="store_true", help="Write filtered cavity PDB per structure (recommended).")
    ap.add_argument("--filtered-suffix", default="cavity_groove.pdb", help="Name for filtered cavity PDB inside each folder.")
    args = ap.parse_args()

    root = Path(args.kvfinder_root).resolve()
    pep = Path(args.peptide_ca_pdb).resolve()
    if not root.is_dir():
        raise SystemExit(f"ERROR: --kvfinder-root not found: {root}")
    if not pep.is_file():
        raise SystemExit(f"ERROR: --peptide-ca-pdb not found: {pep}")

    peptide_ca = read_peptide_ca_coords(pep)

    rows_out = []
    n_struct = 0
    n_with_kept = 0
    

    for sub in sorted([p for p in root.iterdir() if p.is_dir()]):
        results_path = sub / args.results_name
        cavity_path = sub / args.cavity_name

        if not results_path.is_file() or not cavity_path.is_file():
            continue

        n_struct += 1
        results = load_toml(results_path)

        cavity_points = parse_cavity_pdb_points_by_id(cavity_path)
        kept_ids, min_dists = filter_cavities_near_peptide(cavity_points, peptide_ca, args.cutoff)

        if kept_ids:
            n_with_kept += 1

        if args.write_filtered_pdb:
            out_pdb = sub / args.filtered_suffix
            write_filtered_cavity_pdb(cavity_path, out_pdb, set(kept_ids))

        # Add rows to global CSV
        struct_id = sub.name

        if kept_ids:
            for cav_id in kept_ids:
                m = get_metrics_from_results(results, cav_id)
                rows_out.append({
                    "structure": struct_id,
                    "cavity_id": cav_id,
                    "min_dist_to_peptide_ca": min_dists.get(cav_id, None),
                    **m,
                    "results_toml": str(results_path),
                    "cavity_pdb": str(cavity_path),
                })
        else:
            # No cavity matched the peptide proximity rule -> still report structure
            rows_out.append({
                "structure": struct_id,
                "cavity_id": "NA",
                "min_dist_to_peptide_ca": "NA",
                "volume": "NA",
                "area": "NA",
                "max_depth": "NA",
                "avg_depth": "NA",
                "avg_hydropathy": "NA",
                "results_toml": str(results_path),
                "cavity_pdb": str(cavity_path),
            })
            
    # Write global CSV
    out_csv = root / args.out_csv
    fieldnames = [
        "structure", "cavity_id",
        "min_dist_to_peptide_ca",
        "volume", "area", "max_depth", "avg_depth", "avg_hydropathy",
        "results_toml", "cavity_pdb",
    ]
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    print(f"[done] Structures scanned: {n_struct}")
    print(f"[done] Structures with >=1 kept cavity: {n_with_kept}")
    print(f"[done] Kept cavity rows written: {len(rows_out)}")
    print(f"[done] CSV: {out_csv}")


if __name__ == "__main__":
    main()

