"""One read-only Newton defect correction of a captured CFDX face residual.

The original V18 CSR Jacobian is used as an INEXACT inverse. The candidate
pressure is not committed: rerun the native V79 probe to recompute the real
Rhie-Chow residual on the original face connectivity before deciding anything.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.sparse.linalg import splu

from scripts.bfs_pressure_precision_oracle_v79 import extended_residual


def read_split(path, n, nnz):
    with Path(path).open() as f:
        if f.readline().strip() != "CFDX_PRESSURE_SPLIT_V79":
            raise ValueError("invalid V79 split magic")
        if f.readline().split() != [str(n), str(nnz)]:
            raise ValueError("split not for this pressure matrix")
        data = np.loadtxt(f, ndmin=2)
    if data.shape != (n, 2) or not np.all(np.isfinite(data)):
        raise ValueError("invalid split data")
    return data[:, 0], data[:, 1]


def read_defect(path, n, nnz):
    with Path(path).open() as f:
        if f.readline().strip() != "CFDX_PRESSURE_FACE_DEFECT_V79":
            raise ValueError("invalid face-defect magic")
        if f.readline().split() != [str(n), str(nnz)]:
            raise ValueError("defect not for this matrix")
        defect = np.loadtxt(f, ndmin=1)
    if defect.shape != (n,) or not np.all(np.isfinite(defect)):
        raise ValueError("invalid or incomplete physical face defect")
    return defect


def refine(matrix, rhs, split, defect_path, output_dir):
    matrix = Path(matrix)
    rhs = Path(rhs)
    output_dir = Path(output_dir)
    a = mmread(matrix).tocsr()
    a.sum_duplicates()
    a.sort_indices()
    b = np.loadtxt(rhs, ndmin=1)
    n = len(b)
    high, low = read_split(split, n, a.nnz)
    defect = read_defect(defect_path, n, a.nnz)
    lu = splu(a.tocsc(), permc_spec="COLAMD")
    increment = lu.solve(defect)
    if not np.all(np.isfinite(increment)):
        raise ValueError("invalid Newton correction")
    extended = (
        high.astype(np.longdouble)
        + low.astype(np.longdouble)
        + increment.astype(np.longdouble)
    )
    new_high = extended.astype(np.float64)
    new_low = (extended - new_high.astype(np.longdouble)).astype(np.float64)
    _, csr_before = extended_residual(a, b, high, low)
    _, csr_after = extended_residual(a, b, new_high, new_low)
    face_before = float(np.linalg.norm(defect) / np.linalg.norm(b))
    output_dir.mkdir(parents=True, exist_ok=True)
    with (output_dir / "oracle_face_corrected_v79.split").open("w") as f:
        f.write("CFDX_PRESSURE_SPLIT_V79\n")
        f.write(f"{n} {a.nnz}\n")
        for h, l in zip(new_high, new_low):
            f.write(f"{h:.17g} {l:.17g}\n")
    rep = {
        "schema": "CFDX-FACE-DEFECT-CORRECTION-V79/1",
        "csr_before": csr_before,
        "csr_after": csr_after,
        "physical_face_before": face_before,
        "physical_face_after": "REQUIRES_NATIVE_REEVALUATION",
        "inexact_jacobian": True,
        "original_gate": 1e-11,
        "native_linear_certificates": False,
        "production_commit_authorized": False,
    }
    (output_dir / "face_defect_refinement_v79.json").write_text(
        json.dumps(rep, indent=2) + "\n"
    )
    return rep


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("matrix", "rhs", "split", "face-defect", "output-dir"):
        p.add_argument("--" + name, required=True, type=Path)
    args = p.parse_args()
    print(
        json.dumps(
            refine(
                args.matrix, args.rhs, args.split, args.face_defect, args.output_dir
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
