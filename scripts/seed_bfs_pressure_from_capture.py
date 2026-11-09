"""Create a *non-qualified initial guess* from a captured CFDX pressure solve.

Uses the same CSR/RHS as the real post-boundary, post-gauge SIMPLE operator.
This is an initialization diagnostic, not a replacement for the physical
pressure equation or its original true-residual/conservation gates.

The output carries pressure ``fields/scalar/p`` and, if absent, the explicit
initial U vector so that --initial-velocity-field --initial-pressure-field
can import both without changing the mesh, case BCs, or numerical options.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import struct
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _double_from_hex(token: str) -> float:
    return struct.unpack("=d", struct.pack("=Q", int(token, 16)))[0]


def read_capture(path: Path) -> tuple[csr_matrix, np.ndarray, dict]:
    with path.open(encoding="ascii") as stream:
        if stream.readline().strip() != "CFDX_LINEAR_REPLAY_V1":
            raise ValueError("not a CFDX pressure CSR capture")
        metadata = list(map(int, stream.readline().split()))
        rows, columns, nnz, _budget = map(int, stream.readline().split())
        tolerance = _double_from_hex(stream.readline().strip())
        if rows <= 0 or rows != columns or nnz < rows or nnz > rows * 64:
            raise ValueError("invalid pressure capture dimensions")
        row_offsets = np.fromiter(
            (int(stream.readline().strip(), 16) for _ in range(rows + 1)),
            dtype=np.int64,
            count=rows + 1,
        )
        column_ids = np.empty(nnz, dtype=np.int32)
        matrix_values = np.empty(nnz, dtype=np.float64)
        for j in range(nnz):
            column, value = stream.readline().split()
            column_ids[j] = int(column, 16)
            matrix_values[j] = _double_from_hex(value)
        rhs = np.empty(rows, dtype=np.float64)
        initial = np.empty(rows, dtype=np.float64)
        for j in range(rows):
            bvalue, xvalue = stream.readline().split()
            rhs[j] = _double_from_hex(bvalue)
            initial[j] = _double_from_hex(xvalue)
        if (
            row_offsets[0] != 0
            or row_offsets[-1] != nnz
            or np.any(row_offsets[1:] < row_offsets[:-1])
            or np.any(column_ids < 0)
            or np.any(column_ids >= rows)
            or not np.isfinite(matrix_values).all()
            or not np.isfinite(rhs).all()
            or not np.isfinite(initial).all()
        ):
            raise ValueError("invalid or non-finite pressure CSR capture")
        matrix = csr_matrix(
            (matrix_values, column_ids, row_offsets), shape=(rows, rows)
        )
        # Do not silently treat a non-pressure or incompatible capture as a
        # certificate. Metadata is retained for later provenance checks.
        return matrix, rhs, {"metadata": metadata, "target_tolerance": tolerance}


def build_seed(case: Path, capture: Path, output: Path, *, blend: float = 1.0) -> dict:
    if not (np.isfinite(blend) and 0.0 < blend <= 1.0):
        raise ValueError("blend must be within (0,1]")
    case, capture, output = case.resolve(), capture.resolve(), output.resolve()
    if not case.is_file() or not capture.is_file() or output == case or output.exists():
        raise ValueError("input capture/case missing, or output would overwrite a file")
    A, b, meta = read_capture(capture)
    with h5py.File(case) as src:
        cells = int(src.attrs["n_cells"])
        if cells != A.shape[0]:
            raise ValueError("captured pressure matrix does not match mesh cell count")
        if "case_setup_json" not in src.attrs:
            raise ValueError("pressure initialization requires a complete CFDX case")
        setup = json.loads(src.attrs["case_setup_json"])
        raw_u = setup["initial_condition"].get("velocity_vector")
        if (
            not isinstance(raw_u, list)
            or len(raw_u) != 3
            or not all(np.isfinite(component) for component in raw_u)
        ):
            raise ValueError("invalid case-declared initial velocity")
        if (
            "fields" in src
            and "scalar" in src["fields"]
            and "p" in src["fields/scalar"]
        ):
            raise ValueError("input already has a pressure field; refuse overwrite")
    lu = splu(A.tocsc())
    pressure = lu.solve(b)
    pressure *= blend
    if not np.isfinite(pressure).all():
        raise ValueError("non-finite pressure seed")
    residual = A @ pressure - blend * b
    relative = float(np.linalg.norm(residual) / max(np.linalg.norm(blend * b), 1e-300))
    if not np.isfinite(relative) or relative > 1e-5:
        raise ArithmeticError("pressure seed is too inaccurate even for initialization")
    output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(case, output)
    try:
        with h5py.File(output, "r+") as dst:
            fields = dst.require_group("fields")
            scalar = fields.require_group("scalar")
            scalar.create_dataset("p", data=pressure, compression="gzip")
            vectors = fields.require_group("vector")
            if "U" not in vectors:
                velocity = np.broadcast_to(np.array(raw_u)[None, :], (cells, 3))
                vectors.create_dataset("U", data=velocity, compression="gzip")
            dst.attrs["pressure_initialization_status"] = "DIAGNOSTIC_NOT_QUALIFIED"
            dst.attrs["pressure_initialization_capture_sha256"] = sha256(capture)
    except BaseException:
        output.unlink(missing_ok=True)
        raise
    return {
        "status": "DIAGNOSTIC_INITIAL_ONLY_NOT_QUALIFIED",
        "case_sha256": sha256(case),
        "capture_sha256": sha256(capture),
        "output_sha256": sha256(output),
        "pressure_rows": cells,
        "pressure_norm_l2": float(np.linalg.norm(pressure)),
        "pressure_min": float(pressure.min()),
        "pressure_max": float(pressure.max()),
        "seed_linear_true_relative": relative,
        "capture_target_tolerance": meta["target_tolerance"],
        "blend": blend,
        "no_changed_physical_bc_or_tolerance": True,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", type=Path, required=True)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--blend", type=float, default=1.0)
    args = parser.parse_args()
    report = build_seed(args.case, args.capture, args.output, blend=args.blend)
    args.output.with_suffix(".pressure-seed.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
