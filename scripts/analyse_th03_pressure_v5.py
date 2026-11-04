"""Independent SciPy certification of CFDX's ORIGINAL pressure-correction matrix.

Input: CFDX_PRESSURE_DUMP_PREFIX=/path/prefix from the production solver.
The report never substitutes a direct residual for a failed CFDX solve.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.sparse.linalg import LinearOperator, onenormest, splu


def analyse(prefix: Path, target: float = 1e-11) -> dict:
    if not np.isfinite(target) or target <= 0:
        raise ValueError('target must be finite and positive')
    matrix = mmread(str(prefix) + '.mtx').tocsc()
    rhs = np.loadtxt(str(prefix) + '.rhs', ndmin=1)
    if matrix.shape[0] != matrix.shape[1] or matrix.shape[0] != len(rhs):
        raise ValueError('invalid matrix/right-hand-side dimensions')
    if not np.all(np.isfinite(matrix.data)) or not np.all(np.isfinite(rhs)):
        raise ValueError('nonfinite matrix or rhs')
    diag = matrix.diagonal()
    symmetry = matrix - matrix.T
    symmetry_max = float(np.max(np.abs(symmetry.data))) if symmetry.nnz else 0.
    lu = splu(matrix)
    correction = lu.solve(rhs)
    residual = matrix @ correction - rhs
    rhs_norm = np.linalg.norm(rhs)
    if rhs_norm == 0:
        raise ValueError('zero RHS is not a pressure benchmark')
    physical_rel = float(np.linalg.norm(residual) / rhs_norm)
    inverse = LinearOperator(matrix.shape, matvec=lu.solve,
        rmatvec=lambda x: lu.solve(x, trans='T'))
    condition_1 = float(onenormest(matrix) * onenormest(inverse))
    matrix_norm = float(np.max(np.asarray(abs(matrix).sum(axis=1)).ravel()))
    correction_norm = float(np.linalg.norm(correction, np.inf))
    backward_error = float(np.linalg.norm(residual, np.inf) /
        max(matrix_norm * correction_norm + np.linalg.norm(rhs, np.inf), 1e-300))
    return {
        'rows': int(matrix.shape[0]), 'nonzeros': int(matrix.nnz),
        'symmetric_max_absolute_error': symmetry_max,
        'all_positive_diagonal': bool(np.all(diag > 0)),
        'diagonal_dynamic_range': float(diag.max() / diag.min())
            if np.all(diag > 0) else None,
        'condition_1_estimate': condition_1,
        'direct_solver_true_relative_residual': physical_rel,
        'direct_solver_backward_error': backward_error,
        'requested_true_relative_tolerance': target,
        'status': 'LINEAR_CERTIFIED' if physical_rel <= target else 'LINEAR_NOT_CERTIFIED',
        'note': 'Physical benchmark qualification still requires nonlinear conservation, y+, wall stress, and three-grid convergence',
    }


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('prefix', type=Path)
    p.add_argument('--tolerance', type=float, default=1e-11)
    p.add_argument('--output', type=Path)
    args = p.parse_args()
    report = analyse(args.prefix, args.tolerance)
    payload = json.dumps(report, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(payload)
    print(payload, end='')
    return 0 if report['status'] == 'LINEAR_CERTIFIED' else 2


if __name__ == '__main__':
    raise SystemExit(main())
