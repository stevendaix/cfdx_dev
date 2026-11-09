#!/usr/bin/env python3
"""Fail-closed two-sided nullspace and compatible RHS audit of a *captured* CFDX Fp.

Uses LAPACK SVD on a bounded-size CFP01 capture, never modifies a production
matrix/RHS and never promotes a PCD or physical convergence certificate.
No implicit pseudo-inverse is installed into the C++ Schur preconditioner.
For large sparse systems use a separately certified scalable deflation basis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import struct
import sys
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix


class CertificateError(ValueError):
    pass


def read_capture(path: Path, max_rows: int = 1024):
    raw = path.read_bytes()
    if len(raw) < 24 or raw[:8] != b"CFDXFP01":
        raise CertificateError("invalid CFP01 magic or truncated header")
    n, nnz = struct.unpack_from("<QQ", raw, 8)
    if not 1 <= n <= max_rows or not n <= nnz <= 32 * n:
        raise CertificateError("invalid/unsupported CSR dimensions; dense SVD limit exceeded")
    expected = 24 + 4 * (n + 1 + nnz) + 8 * (nnz + 2 * n)
    if len(raw) != expected:
        raise CertificateError("CFP01 payload size mismatch or trailing bytes")
    pos = 24
    def get(dtype, count):
        nonlocal pos
        arr = np.frombuffer(raw, dtype=dtype, count=count, offset=pos).copy()
        pos += np.dtype(dtype).itemsize * count
        return arr
    rows = get("<u4", n + 1)
    cols = get("<u4", nnz)
    vals = get("<f8", nnz)
    rhs = get("<f8", n)
    initial = get("<f8", n)
    if rows[0] != 0 or rows[-1] != nnz or np.any(rows[1:] < rows[:-1]) or np.any(cols >= n):
        raise CertificateError("invalid CSR offsets or columns")
    if not (np.isfinite(vals).all() and np.isfinite(rhs).all() and np.isfinite(initial).all()):
        raise CertificateError("nonfinite operator or vectors")
    for row in range(n):
        segment = cols[rows[row]:rows[row+1]]
        if len(segment) and np.any(segment[1:] <= segment[:-1]):
            raise CertificateError("CSR columns must be strictly sorted and unique")
    A = csr_matrix((vals, cols, rows), shape=(n, n))
    if not np.isfinite(np.linalg.norm(rhs)):
        raise CertificateError("nonfinite RHS norm")
    return A, rhs, initial, hashlib.sha256(raw).hexdigest()


def analyze(path: Path, *, max_rows: int = 1024,
            rank_rtol: float = 1e-12, min_gap: float = 1e3) -> dict:
    if not (0 < rank_rtol < 1e-6 and math.isfinite(rank_rtol)):
        raise CertificateError("invalid numerical rank threshold")
    A, b, initial, sha256 = read_capture(path, max_rows)
    dense = A.toarray()
    u, s, vt = np.linalg.svd(dense, full_matrices=False)
    spectral = float(s[0]); fro = float(np.linalg.norm(dense))
    if spectral == 0 or fro == 0:
        raise CertificateError("zero Fp matrix is not a certified pressure operator")
    threshold = rank_rtol * spectral
    resolved = s > threshold
    rank = int(np.count_nonzero(resolved))
    nullity = len(b) - rank
    if nullity:
        gap = float(s[rank - 1] / max(s[rank], np.finfo(float).tiny)) if rank else 0.
    else:
        gap = math.inf
    right = vt[~resolved].T
    left = u[:, ~resolved]
    right_defect = float(np.linalg.norm(dense @ right) / fro)
    left_defect = float(np.linalg.norm(dense.T @ left) / fro)
    # SVD explicitly supplies an orthonormal basis on BOTH sides. Its left
    # projection is the lower bound on ||Ax-b|| for an exact singular system.
    bnorm = float(np.linalg.norm(b))
    incompatible = left @ (left.T @ b)
    compatible_rhs = b - incompatible
    min_residual = float(np.linalg.norm(incompatible))
    # Construct the minimum-norm *compatible* pseudoinverse only for EVIDENCE.
    coeff = (u[:, resolved].T @ compatible_rhs) / s[resolved]
    x = vt[resolved].T @ coeff
    projected_true = float(np.linalg.norm(dense @ x - compatible_rhs))
    original_true = float(np.linalg.norm(dense @ x - b))
    right_orthogonality = float(np.linalg.norm(right.T @ x))
    numerical_basis_certified = (
        rank > 0 and gap >= min_gap and
        right_defect <= 1e-12 and left_defect <= 1e-12 and
        projected_true <= 1e-8 * max(1., bnorm) and
        np.isfinite(x).all()
    )
    return {
        "schema": "cfdx-pcd-fp-bisided-nullspace/v13",
        "capture_sha256": sha256,
        "matrix": {"rows": len(b), "nnz": A.nnz, "frobenius_norm": fro,
                   "spectral_norm": spectral, "numerical_rank": rank,
                   "numerical_nullity_left": nullity, "numerical_nullity_right": nullity,
                   "smallest_nonzero_singular": float(s[rank - 1]) if rank else 0.,
                   "largest_null_singular": float(s[rank]) if nullity else 0.,
                   "spectral_gap": gap if math.isfinite(gap) else None},
        "rank_relative_tolerance": rank_rtol,
        "left_null_transpose_residual_over_frobenius": left_defect,
        "right_null_residual_over_frobenius": right_defect,
        "rhs_l2": bnorm,
        "incompatible_rhs_component_l2": min_residual,
        "incompatible_rhs_relative": min_residual / bnorm if bnorm else 0.,
        "projected_rhs_l2": float(np.linalg.norm(compatible_rhs)),
        "projected_solution_l2": float(np.linalg.norm(x)),
        "projected_true_residual_l2": projected_true,
        "projected_true_residual_relative_rhs": projected_true / bnorm if bnorm else 0.,
        "original_true_residual_l2": original_true,
        "original_true_residual_relative_rhs": original_true / bnorm if bnorm else 0.,
        "right_gauge_orthogonality": right_orthogonality,
        "numerical_basis_certified": bool(numerical_basis_certified),
        "original_rhs_compatible_at_1e-8": bool(
            numerical_basis_certified and min_residual <= 1e-8 * bnorm),
        "physical_convergence": False,
        "notice": "No physical FV equation or PCD operator modified; SVD is a bounded-size diagnostic",
    }


def export_basis(path: Path, output: Path, *, max_rows: int = 1024) -> None:
    """Portable little-endian CFPB13; C++ rechecks Fp against its exact CSR."""
    A, _, _, _ = read_capture(path, max_rows)
    u, singular, vt = np.linalg.svd(A.toarray(), full_matrices=False)
    rank = int(np.count_nonzero(singular > 1e-12 * singular[0]))
    nullity = len(singular) - rank
    if not (1 <= nullity <= 8):
        raise CertificateError("basis export requires a certified nonempty small nullspace")
    rows = A.indptr.astype("<u4").tobytes()
    cols = A.indices.astype("<u4").tobytes()
    values = A.data.astype("<f8").tobytes()
    h = 14695981039346656037
    for byte in (struct.pack("<QQ", A.shape[0], A.nnz) + rows + cols + values):
        h = ((h ^ byte) * 1099511628211) & ((1 << 64) - 1)
    payload = b"CFDXNB13" + struct.pack("<QQQ", A.shape[0], nullity, h)
    for basis in (u[:, rank:], vt[rank:, :].T):
        for column in range(nullity):
            payload += np.asarray(basis[:, column], dtype="<f8").tobytes()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_bytes(payload)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-rows", type=int, default=1024)
    parser.add_argument("--basis-output", type=Path,
                        help="export certified two-sided basis tied to current CSR")
    parser.add_argument("--require-compatible", action="store_true",
                        help="fail when original Fp rhs cannot satisfy 1e-8")
    args = parser.parse_args(argv)
    try:
        report = analyze(args.capture, max_rows=args.max_rows)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        print(f"PCD_FP_BISIDED_V13 nullity={report['matrix']['numerical_nullity_left']} "
              f"rhs_incompatible={report['incompatible_rhs_relative']:.12g} "
              f"projected_residual={report['projected_true_residual_relative_rhs']:.12g} "
              f"certified={int(report['numerical_basis_certified'])} "
              "physical_qualified=0")
        if not report["numerical_basis_certified"]:
            return 2
        if args.basis_output:
            export_basis(args.capture, args.basis_output, max_rows=args.max_rows)
        if args.require_compatible and not report["original_rhs_compatible_at_1e-8"]:
            return 3
        return 0
    except (OSError, ValueError, np.linalg.LinAlgError) as exc:
        print(f"PCD_FP_BISIDED_V13_ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
