"""Independent pressure-roundoff audit on an *unchanged* CFDX CSR capture.

The long-double accumulator is an experiment on an initial pressure-correction
system; it is NOT a replacement for the production FV pressure equation,
pressure-face flux update, or any of the strict CFDX convergence gates.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.sparse.linalg import splu

if __package__:
    from .seed_bfs_pressure_from_capture import read_capture, sha256
else:
    from seed_bfs_pressure_from_capture import read_capture, sha256


def extended_residual(A, b, x):
    """Compute b - A*x with explicit extended-precision CSR products and sum."""
    if A.shape[0] != A.shape[1] or b.size != A.shape[0] or x.size != b.size:
        raise ValueError("incompatible CSR/rhs/solution shapes")
    if not (np.isfinite(A.data).all() and np.isfinite(b).all() and np.isfinite(x).all()):
        raise ValueError("non-finite pressure data")
    wide = np.longdouble
    if np.finfo(wide).eps >= np.finfo(np.float64).eps:
        raise RuntimeError("this platform has no usable extended-precision accumulator")
    # CSR row_offsets are strictly increasing for the pressure matrices
    # emitted by CFDX. For an optional arbitrary CSR test with empty rows,
    # process only nonempty rows so reduceat cannot cross a row boundary.
    products = A.data.astype(wide) * x.astype(wide)[A.indices]
    out = np.zeros(b.size, dtype=wide)
    nonempty = np.flatnonzero(np.diff(A.indptr))
    if nonempty.size:
        out[nonempty] = np.add.reduceat(products, A.indptr[nonempty])
    return b.astype(wide) - out


def audit(capture: Path, passes: int = 2) -> dict:
    if passes < 0 or passes > 8:
        raise ValueError("refinement passes must be in [0,8]")
    A, b, meta = read_capture(capture)
    lu = splu(A.tocsc())
    x = lu.solve(b).astype(np.longdouble)
    ref = np.linalg.norm(b.astype(np.longdouble))
    if not ref or not np.isfinite(ref):
        raise ValueError("zero/nonfinite pressure RHS norm is not supported in this audit")
    results = []
    for step in range(passes + 1):
        wide_r = extended_residual(A, b, x)
        restored_double = np.asarray(x, dtype=np.float64)
        results.append({
            "pass": step,
            "extended_precision_true_relative": float(np.linalg.norm(wide_r) / ref),
            "stored_double_true_relative": float(
                np.linalg.norm(b - A @ restored_double) / np.linalg.norm(b)
            ),
            "extended_pressure_norm_l2": float(np.linalg.norm(x)),
        })
        if step != passes:
            # Factorization remains the same: only residual accumulation and
            # solution update have increased precision in this reference audit.
            dx = lu.solve(np.asarray(wide_r, dtype=np.float64))
            if not np.isfinite(dx).all():
                raise ValueError("direct correction contains non-finite entries")
            x += dx.astype(np.longdouble)
    return {
        "status": "DIAGNOSTIC_ONLY_NOT_PHYSICAL_CERTIFICATE",
        "capture_sha256": sha256(capture),
        "rows": int(A.shape[0]),
        "nonzeros": int(A.nnz),
        "rhs_norm_l2": float(np.linalg.norm(b)),
        "target_tolerance": meta["target_tolerance"],
        "longdouble_mantissa_bits": int(np.finfo(np.longdouble).nmant),
        "steps": results,
        "production_double_pressure_correction_unchanged": True,
        "physical_convergence_proven": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--passes", type=int, default=2)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit(args.capture, args.passes)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(json.dumps(report, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
