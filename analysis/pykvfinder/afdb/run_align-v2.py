#!/usr/bin/env python3
"""
Batch-align all PDBs in a directory onto a reference using TM-align (matrix file)
and apply the rigid transform with Biopython, saving *_aligned.pdb.

Key points:
- We run: TMalign mobile.pdb ref.pdb -m <matrix_file>
- The matrix file may contain either:
    (A) 3x4 rows:  t u11 u12 u13
                  t u21 u22 u23
                  t u31 u32 u33
  or (B) 5 columns with an index first (like your example):
                  m t u(m,1) u(m,2) u(m,3)
- TM-align defines: X = t + U * x   (column-vector convention)
- Biopython Atom.transform expects right-multiplying rotation:
    X = x * R + t
  so we must pass R = U.T and t unchanged.

Usage:
  python align_tmalign_biopython_full.py --ref ref.pdb --in-dir pdbs/ --out-dir aligned/ --debug
"""

import argparse
import subprocess
from pathlib import Path

import numpy as np
from Bio.PDB import PDBParser, PDBIO


def run_tmalign_write_matrix(tmalign_exe: str, mobile_pdb: str, ref_pdb: str, matrix_path: str) -> str:
    """Run TM-align and write the matrix to matrix_path using -m. Returns stdout."""
    cmd = [tmalign_exe, mobile_pdb, ref_pdb, "-m", matrix_path]
    p = subprocess.run(cmd, capture_output=True, text=True)
    if p.returncode != 0:
        raise RuntimeError(
            f"TM-align failed\nCMD: {' '.join(cmd)}\nSTDERR:\n{p.stderr}\nSTDOUT:\n{p.stdout}"
        )
    return p.stdout


def rotation_quality(U: np.ndarray):
    """det(U) should be ~+1; U U^T should be ~I."""
    det = float(np.linalg.det(U))
    ortho_err = float(np.linalg.norm(U @ U.T - np.eye(3)))
    return det, ortho_err

def read_tmalign_matrix_file(matrix_path: str):
    """
    Parse TM-align matrix file written by '-m'.
    Robust: only reads the 3 lines AFTER the rotation-matrix header.
    """

    with open(matrix_path, "r") as f:
        lines = f.readlines()

    start = None
    for i, line in enumerate(lines):
        if "rotation matrix" in line.lower():
            start = i + 1
            break

    if start is None:
        raise ValueError("Rotation matrix header not found in TM-align output.")

    matrix_rows = []

    for line in lines[start:]:
        parts = line.strip().split()

        # We expect either:
        # 1 t u11 u12 u13
        # or
        # t u11 u12 u13

        try:
            nums = [float(x) for x in parts]
        except ValueError:
            continue

        if len(nums) == 5:
            # indexed format
            matrix_rows.append(nums[1:5])
        elif len(nums) == 4:
            # plain format
            matrix_rows.append(nums)
        else:
            continue

        if len(matrix_rows) == 3:
            break

    if len(matrix_rows) != 3:
        raise ValueError("Failed to extract 3 rotation rows from matrix file.")

    M = np.array(matrix_rows, dtype=float)  # shape (3,4)
    tvec = M[:, 0]
    U = M[:, 1:4]

    return U, tvec


def transform_structure_biopython(in_pdb: str, out_pdb: str, U: np.ndarray, t: np.ndarray):
    """
    TM-align:   X = t + U * x    (column-vector)
    Biopython:  X = x * R + t    (right-multiplying rotation)
    So: R = U.T
    """
    R = U.T

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("mob", in_pdb)

    for atom in structure.get_atoms():
        atom.transform(R, t)

    io = PDBIO()
    io.set_structure(structure)
    io.save(out_pdb)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ref", required=True, help="Reference PDB (Chain_2).")
    ap.add_argument("--in-dir", required=True, help="Directory with input PDBs to align (Chain_1).")
    ap.add_argument("--out-dir", required=True, help="Directory for *_aligned.pdb outputs.")
    ap.add_argument("--tmalign", default="TMalign", help="TM-align executable (default: TMalign).")
    ap.add_argument("--pattern", default="*.pdb", help="Glob pattern for input PDBs (default: *.pdb).")
    ap.add_argument("--skip-existing", action="store_true", help="Skip if *_aligned.pdb exists.")
    ap.add_argument("--debug", action="store_true", help="Print matrix diagnostics.")
    args = ap.parse_args()

    ref = Path(args.ref).resolve()
    in_dir = Path(args.in_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    # Always save matrices here
    mat_dir = out_dir / "matrices"
    mat_dir.mkdir(parents=True, exist_ok=True)

    pdbs = sorted(in_dir.glob(args.pattern))
    if not pdbs:
        raise SystemExit(f"No files matched {args.pattern} in {in_dir}")

    for pdb in pdbs:
        pdb = pdb.resolve()
        if pdb == ref:
            continue

        out_pdb = out_dir / f"{pdb.stem}_aligned.pdb"
        if args.skip_existing and out_pdb.exists():
            print(f"[skip] {pdb.name}")
            continue

        matrix_path = mat_dir / f"{pdb.stem}_tmalign.matrix.txt"

        try:
            _stdout = run_tmalign_write_matrix(args.tmalign, str(pdb), str(ref), str(matrix_path))

            U, t = read_tmalign_matrix_file(str(matrix_path))
            det, ortho_err = rotation_quality(U)

            if args.debug:
                print(f"\n=== {pdb.name} ===")
                print(f"Matrix saved to: {matrix_path}")
                print("U (as parsed):\n", U)
                print("t (as parsed):\n", t)
                print(f"det(U) = {det:.6f}")
                print(f"||U U^T - I|| = {ortho_err:.3e}")
                if ortho_err > 1e-2:
                    print("WARNING: Rotation is not orthonormal -> parsing is likely wrong for this TM-align build.")

            transform_structure_biopython(str(pdb), str(out_pdb), U, t)
            print(f"[ok]   {pdb.name} -> {out_pdb.name}")

        except Exception as e:
            print(f"[err]  {pdb.name}: {e}")
            print(f"       Matrix (if written) is at: {matrix_path}")

if __name__ == "__main__":
    main()

