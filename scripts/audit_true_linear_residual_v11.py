"""Independent, original-system linear-residual certificate (never marks CFD as qualified).

Use a CFDX_PRESSURE_DUMP_PREFIX .mtx/.rhs pair. Recompute the direct-solver
residual in two accumulation precisions and distinguish residual norms from
backward error. An accurate backward error does NOT establish mass conservation.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.sparse.linalg import splu


def audit(prefix: Path, target: float = 1.0e-11) -> dict:
    if not np.isfinite(target) or target <= 0:
        raise ValueError("target must be positive and finite")
    matrix = mmread(str(prefix) + ".mtx").tocsr()
    rhs = np.loadtxt(str(prefix) + ".rhs", ndmin=1)
    n = matrix.shape[0]
    if n == 0 or matrix.shape != (n, n) or rhs.shape != (n,):
        raise ValueError("matrix or RHS has invalid shape")
    if not np.all(np.isfinite(matrix.data)) or not np.all(np.isfinite(rhs)):
        raise ValueError("matrix or RHS is not finite")
    rhs_norm2 = float(np.linalg.norm(rhs))
    if not (rhs_norm2 > 0):
        raise ValueError("zero RHS: relative benchmark undefined")
    x = splu(matrix.tocsc()).solve(rhs)
    r64 = (matrix @ x) - rhs
    # Recompute A*x-b with extended-precision row products and summation.
    # The factors and x remain FP64; this is a residual-audit, not an FP128 solve.
    rlong = np.zeros(n, dtype=np.longdouble)
    componentwise_denominator = np.zeros(n, dtype=np.longdouble)
    for i in range(n):
        begin, end = matrix.indptr[i : i + 2]
        a = matrix.data[begin:end].astype(np.longdouble)
        xi = x[matrix.indices[begin:end]].astype(np.longdouble)
        rlong[i] = np.sum(a * xi, dtype=np.longdouble) - np.longdouble(rhs[i])
        componentwise_denominator[i] = np.abs(np.longdouble(rhs[i])) + np.sum(
            np.abs(a * xi), dtype=np.longdouble
        )
    relative64 = float(np.linalg.norm(r64) / rhs_norm2)
    relative_ext = float(
        np.sqrt(np.sum(rlong * rlong, dtype=np.longdouble)) / np.longdouble(rhs_norm2)
    )
    matrix_inf = float(np.max(np.asarray(abs(matrix).sum(axis=1)).ravel()))
    x_inf = float(np.linalg.norm(x, np.inf))
    b_inf = float(np.linalg.norm(rhs, np.inf))
    residual_inf = float(np.max(np.abs(rlong)))
    normwise_backward = residual_inf / max(
        matrix_inf * x_inf + b_inf, np.finfo(float).tiny
    )
    positive_den = componentwise_denominator > 0
    componentwise_error = (
        float(
            np.max(
                np.abs(rlong[positive_den]) / componentwise_denominator[positive_den]
            )
        )
        if np.any(positive_den)
        else 0.0
    )
    if np.any(~positive_den & (rlong != 0)):
        componentwise_error = float("inf")
    return {
        "matrix_size": n,
        "nonzero_entries": int(matrix.nnz),
        "float64_l2_relative_true_residual": relative64,
        "extended_accumulation_l2_relative_true_residual": relative_ext,
        "true_residual_inf": residual_inf,
        "normwise_backward_error": normwise_backward,
        "componentwise_backward_error": componentwise_error,
        "target_relative_residual": target,
        "linear_status": "CERTIFIED" if relative_ext <= target else "NOT_CERTIFIED",
        "physical_status": "NOT_EVALUATED",
        "warning": "Direct solve is not exact; extended accumulation does not change the precision of x or A",
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("prefix", type=Path)
    p.add_argument("--tolerance", type=float, default=1.0e-11)
    p.add_argument("--json", type=Path)
    a = p.parse_args()
    result = audit(a.prefix, a.tolerance)
    content = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if a.json:
        a.json.parent.mkdir(parents=True, exist_ok=True)
        a.json.write_text(content)
    print(content, end="")
    return 0 if result["linear_status"] == "CERTIFIED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
