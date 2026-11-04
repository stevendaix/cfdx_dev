#!/usr/bin/env python3
"""Independent sparse LU + iterative refinement precision audit, NOT a pressure solver patch.

The certificate is the unscaled true Euclidean residual ||A*x-b||2/||b||2.
High-precision summation is also reported, but DOES NOT replace the FP64
production residual gate, nor any mass-conservation or physical criterion.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.sparse.linalg import splu


def inspect(mtx: Path, rhs_path: Path, tolerance: float = 1e-11,
            iterations: int = 8) -> dict:
    if not (math.isfinite(tolerance) and tolerance > 0) or iterations < 0:
        raise ValueError('invalid pressure audit controls')
    A = mmread(str(mtx)).tocsr()
    b = np.atleast_1d(np.loadtxt(rhs_path, dtype=float))
    if A.shape[0] != A.shape[1] or A.shape[0] != b.size or not np.all(np.isfinite(b)):
        raise ValueError('pressure matrix and rhs inconsistent')
    if not np.all(np.isfinite(A.data)) or not np.linalg.norm(b):
        raise ValueError('pressure audit requires a finite nonzero right-hand side')
    lu = splu(A.tocsc())
    x = lu.solve(b)
    scale = np.linalg.norm(b)
    def extended_residual(solution: np.ndarray) -> np.ndarray:
        xx = solution.astype(np.longdouble)
        bb = b.astype(np.longdouble)
        vals = A.data.astype(np.longdouble)
        for row in range(A.shape[0]):
            a, z = A.indptr[row:row+2]
            bb[row] -= np.sum(vals[a:z] * xx[A.indices[a:z]], dtype=np.longdouble)
        return bb
    observations = []
    for step in range(iterations + 1):
        native = b - A @ x
        more_precise = extended_residual(x)
        ratio_fp64 = float(np.linalg.norm(native) / scale)
        ratio_long = float(np.sqrt(np.sum(more_precise**2, dtype=np.longdouble)) / scale)
        observations.append({
            'refinement': step,
            'relative_true_residual_fp64': ratio_fp64,
            'relative_true_residual_extended': ratio_long,
            'passed_original_fp64_gate': ratio_fp64 <= tolerance,
        })
        if step < iterations:
            correction = lu.solve(more_precise.astype(float))
            if not np.all(np.isfinite(correction)):
                raise ValueError('pressure refinement became nonfinite')
            x += correction
    min_fp64 = min(row['relative_true_residual_fp64'] for row in observations)
    return {
        'status': 'DIAGNOSTIC_ONLY_NOT_QUALIFIED',
        'n': int(A.shape[0]),
        'nonzeros': int(A.nnz),
        'tolerance': tolerance,
        'initial_true_residual_fp64': observations[0]['relative_true_residual_fp64'],
        'best_true_residual_fp64': min_fp64,
        'initial_gate_passed': observations[0]['passed_original_fp64_gate'],
        'any_refinement_passed_original_gate': any(r['passed_original_fp64_gate'] for r in observations),
        'observations': observations,
        'physical_qualification': 'NOT_QUALIFIED',
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--matrix', required=True, type=Path)
    parser.add_argument('--rhs', required=True, type=Path)
    parser.add_argument('--tolerance', type=float, default=1e-11)
    parser.add_argument('--iterations', type=int, default=8)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    report = inspect(args.matrix, args.rhs, args.tolerance, args.iterations)
    payload = json.dumps(report, indent=2) + '\n'
    if args.output:
        args.output.write_text(payload, encoding='utf-8')
    print(payload, end='')


if __name__ == '__main__':
    main()