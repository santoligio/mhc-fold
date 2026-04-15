#!/usr/bin/env python3
"""
Auto-detect alpha1/alpha2-like groove helices in MHC-like proteins using Biopython + DSSP (mkdssp),
WITHOUT requiring a chain ID.

Single merge parameter:
- max_break = maximum gap (in residues) allowed between sequential helix segments to be merged into
  the same helix COMPONENT.
  Example: helix segments (118-128) and (133-154) have gap = 133-128-1 = 4.
  If --max_break >= 4 (e.g., 20), they are merged into one component (118-154).

If after merging we still have > max_components helix components, we SKIP pair/geometry computations
and report status=TOO_MANY_COMPONENTS in the CSV.

Outputs include:
- helix1/helix2 ranges + lengths (from the two selected helix components)
- sep_res (gap between components, in residue-index space)
- helix_com_dist (Å), helix_mean_nn_ca_dist (Å)
- whole-chain helix % and groove helix contribution

Requirements:
  conda install -c conda-forge biopython dssp numpy
"""

from __future__ import annotations

import os
import re
import csv
import shutil
import argparse
import tempfile
import pandas as pd
from pathlib import Path
from typing import Dict, Tuple, Optional, List, Any

import numpy as np
from Bio.PDB import PDBParser, MMCIFIO, Select
from Bio.PDB.DSSP import DSSP
from Bio.PDB.Polypeptide import is_aa


DSSP_HELIX = {"H", "G", "I"}
DUMMY_CRYST1 = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1\n"

VALID_ONE = {"C", "N", "O", "H", "S", "P", "F", "I", "K"}
VALID_TWO = {
    "BR","CL","NA","MG","ZN","FE","MN","CU","NI","CO","SI","AL","SR","CD","CS","HG","PT",
    "PB","SB","AG","AU","BA","BI","TI","ZR","OS","IR","MO","SE","KR","XE","RB","V","CR",
    "W","U","LI","BE","AS","GA","GE","SN","TE"
}


# ----------------------------
# Parsing + sanitization
# ----------------------------

def load_structure(path: Path):
    parser = PDBParser(QUIET=True, PERMISSIVE=True)
    return parser.get_structure(path.stem, str(path))


def _infer_element_from_name(atom_name: str) -> Optional[str]:
    if not atom_name:
        return None
    name = atom_name.strip().upper()
    lead = name[0]
    if lead in {"C", "N", "O", "S", "H", "P"}:
        return lead
    cand2 = name[:2]
    cand1 = name[:1]
    if cand2 in VALID_TWO:
        return cand2
    if cand1 in VALID_ONE:
        return cand1
    return None


class _DropUnknownElementSelect(Select):
    def __init__(self, drop_unknown: bool):
        super().__init__()
        self.drop_unknown = drop_unknown

    def accept_atom(self, atom):
        el = (atom.element or "").strip().upper()
        ok = (len(el) == 1 and el in VALID_ONE) or (len(el) == 2 and el in VALID_TWO)
        if not ok:
            inferred = _infer_element_from_name(atom.get_name())
            if inferred:
                atom.element = inferred
                return True
            return not self.drop_unknown
        return True


def _write_sanitized_temp_cif(original_path: Path, model_index: int, drop_unknown: bool) -> Path:
    structure = load_structure(original_path)
    models = list(structure.get_models())

    if model_index < 0 or model_index >= len(models):
        raise RuntimeError(f"Model index {model_index} out of range for {original_path.name}")

    model = models[model_index]

    io = MMCIFIO()
    io.set_structure(model)

    # Use NamedTemporaryFile and close fd immediately
    
    tmp_file = tempfile.NamedTemporaryFile(prefix="sanitized_", suffix=".cif", delete=False)
    tmp_path = Path(tmp_file.name)
    tmp_file.close()  # closes the open fd to avoid "too many open files"                                                                                    

    io.save(str(tmp_path), select=_DropUnknownElementSelect(drop_unknown))
    return tmp_path


# ----------------------------
# Geometry helpers
# ----------------------------

def get_ca_coords_for_reskeys(chain, reskeys_set: set) -> np.ndarray:
    coords = []
    for res in chain.get_residues():
        if res.id[0] != " ":
            continue
        key = (int(res.id[1]), (res.id[2] or "").strip())
        if key not in reskeys_set:
            continue
        if "CA" in res:
            coords.append(res["CA"].get_coord())
    if not coords:
        return np.zeros((0, 3), dtype=float)
    return np.asarray(coords, dtype=float)


def com(coords: np.ndarray) -> np.ndarray:
    if coords.size == 0:
        return np.array([np.nan, np.nan, np.nan], dtype=float)
    return coords.mean(axis=0)


def mean_nearest_neighbor_dist(A: np.ndarray, B: np.ndarray) -> float:
    if A.size == 0 or B.size == 0:
        return float("nan")
    diff = A[:, None, :] - B[None, :, :]
    dist = np.sqrt((diff * diff).sum(axis=2))
    return float(dist.min(axis=1).mean())


# ----------------------------
# Helix segment logic
# ----------------------------

def contiguous_segments(indices: List[int]) -> List[Tuple[int, int]]:
    if not indices:
        return []
    indices = sorted(indices)
    segs = []
    start = prev = indices[0]
    for i in indices[1:]:
        if i == prev + 1:
            prev = i
        else:
            segs.append((start, prev))
            start = prev = i
    segs.append((start, prev))
    return segs


def merge_segments_with_gaps(segs: List[Tuple[int, int]], max_gap: int) -> List[Tuple[int, int]]:
    """
    Merge sequential segments if gap <= max_gap, where:
      gap = next_start - prev_end - 1
    """
    if not segs:
        return []
    segs = sorted(segs)
    merged = [segs[0]]
    for s, e in segs[1:]:
        ps, pe = merged[-1]
        gap = s - pe - 1
        if gap <= max_gap:
            merged[-1] = (ps, e)
        else:
            merged.append((s, e))
    return merged


def pick_two_components_max_total_len(segs: List[Tuple[int, int]]) -> Optional[Tuple[Tuple[int, int], Tuple[int, int]]]:
    """
    Pick the pair of components that maximizes total length.
    Tie-break: larger separation, then earlier start.
    """
    if len(segs) < 2:
        return None

    best = None
    best_score = None

    for i in range(len(segs)):
        for j in range(i + 1, len(segs)):
            s1, e1 = segs[i]
            s2, e2 = segs[j]
            if s2 < s1:
                (s1, e1, s2, e2) = (s2, e2, s1, e1)

            len1 = e1 - s1 + 1
            len2 = e2 - s2 + 1
            total_len = len1 + len2
            sep_res = s2 - e1

            score = (total_len, sep_res, -s1)
            if best is None or score > best_score:
                best = ((s1, e1), (s2, e2))
                best_score = score

    return best


def reskey_label(reskey: Tuple[int, str]) -> str:
    resseq, icode = reskey
    return f"{resseq}{icode}".strip()


def components_to_string(segs: List[Tuple[int, int]]) -> str:
    # 0-based index segments for debugging
    return ";".join([f"{s}-{e}" for s, e in segs])


# ----------------------------
# Chain evaluation
# ----------------------------

def evaluate_chain(model, dssp, chain_id: str,
                   min_chain_len: int,
                   min_len: int,
                   min_merged_len: int,
                   max_break: int,
                   max_components: int) -> Optional[Dict[str, Any]]:
    chain_obj = None
    for ch in model:
        if ch.id == chain_id:
            chain_obj = ch
            break
    if chain_obj is None:
        return None

    residues = [r for r in chain_obj.get_residues() if is_aa(r, standard=True)]
    if len(residues) < min_chain_len:
        return None

    ordered = [(int(r.id[1]), (r.id[2] or "").strip()) for r in residues]
    ordered_set = set(ordered)

    ss_by_res: Dict[Tuple[int, str], str] = {}
    for (ch, resid), d in dssp.property_dict.items():
        hetfield, resseq, icode = resid
        if ch != chain_id or hetfield != " ":
            continue
        key = (int(resseq), (icode or "").strip())
        if key in ordered_set:
            ss_by_res[key] = d[2]

    n_res_total = len(ordered)
    helix_count_total = sum(1 for key in ordered if ss_by_res.get(key, "C") in DSSP_HELIX)
    helix_pct_total = round(100.0 * helix_count_total / n_res_total, 2) if n_res_total else 0.0

    helix_pos = [i for i, key in enumerate(ordered) if ss_by_res.get(key, "C") in DSSP_HELIX]
    raw_segs = contiguous_segments(helix_pos)

    # Filter short RAW segments
    raw_segs = [s for s in raw_segs if (s[1] - s[0] + 1) >= min_len]
    
    # Merge segments
    merged = merge_segments_with_gaps(raw_segs, max_gap=max_break)
    
    #Filter short MERGED components
    segs = [s for s in merged if (s[1] - s[0] + 1) >= min_merged_len]

    if len(segs) < 2:
        return None

    
    # Single merge step: merge helix segments into components if gap <= max_break
    if len(segs) < 2:
        return None

    components_str = components_to_string(segs)

    # If too many components, skip computations and report
    if len(segs) > max_components:
        return {
            "chain": chain_id,
            "status_override": "TOO_MANY_COMPONENTS",
            "n_components": len(segs),
            "components_idx": components_str,
            "n_res_total": n_res_total,
            "helix_count_total": helix_count_total,
            "helix_pct_total": helix_pct_total,
        }

    picked = pick_two_components_max_total_len(segs)
    if picked is None:
        return None
    (h1s, h1e), (h2s, h2e) = picked

    helix1_len = (h1e - h1s + 1)
    helix2_len = (h2e - h2s + 1)
    groove_helix_count = helix1_len + helix2_len
    groove_helix_pct = round(100.0 * groove_helix_count / n_res_total, 2) if n_res_total else 0.0

    sep_res = h2s - h1e

    helix1_keys = set(ordered[h1s:h1e + 1])
    helix2_keys = set(ordered[h2s:h2e + 1])

    h1_ca = get_ca_coords_for_reskeys(chain_obj, helix1_keys)
    h2_ca = get_ca_coords_for_reskeys(chain_obj, helix2_keys)

    h1_com = com(h1_ca)
    h2_com = com(h2_ca)

    helix_com_dist = float(np.linalg.norm(h1_com - h2_com)) \
        if np.isfinite(h1_com).all() and np.isfinite(h2_com).all() else float("nan")
    helix_mean_nn_ca_dist = mean_nearest_neighbor_dist(h1_ca, h2_ca)

    helix_com_dist_out = round(helix_com_dist, 3) if helix_com_dist == helix_com_dist else ""
    helix_mean_nn_out = round(helix_mean_nn_ca_dist, 3) if helix_mean_nn_ca_dist == helix_mean_nn_ca_dist else ""

    return {
        "chain": chain_id,
        "status_override": "",
        "n_components": len(segs),
        "components_idx": components_str,

        "helix1_start": reskey_label(ordered[h1s]),
        "helix1_end": reskey_label(ordered[h1e]),
        "helix2_start": reskey_label(ordered[h2s]),
        "helix2_end": reskey_label(ordered[h2e]),
        "helix1_len": helix1_len,
        "helix2_len": helix2_len,

        "sep_res": sep_res,
        "sep_penalty": 0,
        "n_helix_segments_kept": len(segs),

        "helix_com_dist": helix_com_dist_out,
        "helix_mean_nn_ca_dist": helix_mean_nn_out,

        "n_res_total": n_res_total,
        "helix_count_total": helix_count_total,
        "helix_pct_total": helix_pct_total,
        "groove_helix_count": groove_helix_count,
        "groove_helix_pct": groove_helix_pct,
    }


def pick_best_chain(model, dssp, args) -> Optional[Dict[str, Any]]:
    best = None
    for chain in model.get_chains():
        cand = evaluate_chain(
            model=model, dssp=dssp, chain_id=chain.id,
            min_chain_len=args.min_chain_len,
            min_len=args.min_len,
            min_merged_len=args.min_merged_len,
            max_break=args.max_break,
            max_components=args.max_components,
        )

        if cand is None:
            continue

        cand_status = cand.get("status_override", "")
        if best is None:
            best = cand
            continue

        best_status = best.get("status_override", "")

        if best_status == "TOO_MANY_COMPONENTS" and cand_status != "TOO_MANY_COMPONENTS":
            best = cand
            continue
        if cand_status == "TOO_MANY_COMPONENTS" and best_status != "TOO_MANY_COMPONENTS":
            continue

        if cand_status != "TOO_MANY_COMPONENTS" and best_status != "TOO_MANY_COMPONENTS":
            if (cand.get("groove_helix_count", -1) > best.get("groove_helix_count", -1)) or (
                cand.get("groove_helix_count", -1) == best.get("groove_helix_count", -1)
                and cand.get("helix_pct_total", -1) > best.get("helix_pct_total", -1)
            ):
                best = cand
        else:
            if (cand.get("n_components", 999) < best.get("n_components", 999)) or (
                cand.get("n_components", 999) == best.get("n_components", 999)
                and cand.get("helix_pct_total", -1) > best.get("helix_pct_total", -1)
            ):
                best = cand

    return best


def iter_pdb_files(in_dir: Path, glob_pat: str) -> List[Path]:
    return sorted(in_dir.glob(glob_pat))


def process_file(path: Path, args) -> Dict[str, Any]:
    structure = load_structure(path)
    models = list(structure.get_models())
    if args.model < 0 or args.model >= len(models):
        return {"file": path.name, "status": "FAILED", "error": f"Model {args.model} not found"}

    model = models[args.model]

    tmp_cif = _write_sanitized_temp_cif(path, model_index=args.model, drop_unknown=args.drop_unknown_elements)
    try:
        dssp = DSSP(model, str(tmp_cif), dssp=args.dssp_exec)
    except Exception as e:
        return {"file": path.name, "status": "DSSP_FAILED", "error": str(e)}
    finally:
        try:
            os.unlink(tmp_cif)
        except OSError:
            pass

    best = pick_best_chain(model, dssp, args)
    if best is None:
        return {"file": path.name, "status": "NO_PAIR_FOUND", "error": ""}

    status_override = best.get("status_override", "")
    if status_override:
        return {
            "file": path.name,
            "status": status_override,
            "chain": best.get("chain", ""),
            "n_components": best.get("n_components", ""),
            "components_idx": best.get("components_idx", ""),
            "n_res_total": best.get("n_res_total", ""),
            "helix_count_total": best.get("helix_count_total", ""),
            "helix_pct_total": best.get("helix_pct_total", ""),
            "error": "",
        }

    return {"file": path.name, "status": "OK", **best, "error": ""}


def main():
    ap = argparse.ArgumentParser(description="Auto-detect MHC groove helices (alpha1/alpha2) using DSSP with sanitization.")
    ap.add_argument("--in_dir", required=True, help="Directory with PDB files.")
    ap.add_argument("--out_csv", required=True, help="Output CSV.")
    ap.add_argument("--index_csv", required=True, help="CSV file containing PDB IDs to include.")

    ap.add_argument("--glob", default="*.pdb", help="Glob pattern (default: *.pdb).")
    ap.add_argument("--model", type=int, default=0, help="Model index (default: 0).")
    ap.add_argument("--dssp_exec", default="mkdssp", help="mkdssp executable (default: mkdssp).")
    ap.add_argument("--drop_unknown_elements", action="store_true",
                    help="Drop atoms whose element cannot be inferred (prevents mkdssp failures).")

    ap.add_argument("--min_chain_len", type=int, default=60, help="Minimum standard-AA chain length to consider.")
    ap.add_argument("--min_len", type=int, default=12, help="Minimum helix segment length (before merging).")
    ap.add_argument("--min_merged_len", type=int, default=15, help="Minimum helix component length AFTER merging (default: 15).")
    
    # Now: max_break is the single parameter controlling merging of sequential helix segments into components
    ap.add_argument("--max_break", type=int, default=20,
                    help="Maximum gap (in residues) between sequential helix segments to merge into the same component.")

    ap.add_argument("--max_components", type=int, default=3,
                    help="If > this number of helix components remain after merging, skip calculations and report status.")

    # kept for compatibility (unused)
    ap.add_argument("--expected_sep", type=int, default=50, help="(Unused now) kept for compatibility.")
    ap.add_argument("--sep_min", type=int, default=20, help="(Unused now) kept for compatibility.")
    ap.add_argument("--sep_max", type=int, default=80, help="(Unused now) kept for compatibility.")

    args = ap.parse_args()

    if shutil.which(args.dssp_exec) is None:
        raise SystemExit(f"ERROR: '{args.dssp_exec}' not found in PATH. Try: conda install -c conda-forge dssp")

    indir = Path(args.in_dir)

    # Read assemblies CSV 
    index_df = pd.read_csv(args.index_csv)
    
    # Normalize IDs for matching
    def normalize_id(name: str):
        name = str(name).strip()

        # AFDB → extract AF-UNIPROT
        m = re.search(r"(AF-[A-Z0-9]+)", name.upper())
        if m:
            return m.group(1)

        # PDB → first 4 characters
        return name[:4].lower()
        
    index_df["norm_id"] = index_df["pdb"].apply(normalize_id)
    allowed_ids = set(index_df["norm_id"])
        
    # Get all PDB files
    all_files = iter_pdb_files(indir, args.glob)
        
    files = [f for f in all_files if normalize_id(f.stem) in allowed_ids]    
    rows = [process_file(f, args) for f in files]

    fieldnames = [
        "file", "status", "chain",
        "n_components", "components_idx",

        "helix1_start", "helix1_end", "helix1_len",
        "helix2_start", "helix2_end", "helix2_len",
        "sep_res", "sep_penalty", "n_helix_segments_kept",

        "helix_com_dist", "helix_mean_nn_ca_dist",

        "n_res_total", "helix_count_total", "helix_pct_total",
        "groove_helix_count", "groove_helix_pct",

        "error"
    ]

    out = Path(args.out_csv)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    ok = sum(1 for r in rows if r.get("status") == "OK")
    print(f"Wrote {len(rows)} rows to {out} ({ok} OK)")


if __name__ == "__main__":
    main()

