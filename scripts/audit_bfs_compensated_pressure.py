"""Compensated pressure-operator replay on an *actual captured* CSR.

This is NOT a pressure solver for CFDX production, nor a replacement for
Rhie-Chow face-flux validation. It shows whether FP64 rounding of the stored
pressure correction explains a pressure residual floor, while preserving the
original matrix, RHS, tolerance, and positive/negative qualification gates.

The algebraic face balance is a symmetric M-matrix graph decomposition of
A*x. It is not an independent geometric FV face reconstruction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import splu

from seed_bfs_pressure_from_capture import read_capture


def _sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def compensated_residual(
    a: csr_matrix, b: np.ndarray, hi: np.ndarray, lo: np.ndarray
) -> tuple[np.ndarray, np.ndarray, bool]:
    """Compute b-A*(hi+lo) without forming a rounded hi+lo vector.

    Laplacian row differences and row sums are evaluated separately. A is
    symmetrical for the BFS pressure capture, but residual computation itself
    also supports non-symmetric matrices. No external scaling changes A or b.
    """
    n = b.size
    if (a.shape != (n, n) or hi.shape != (n,) or lo.shape != (n,)
        or not all(np.isfinite(array).all() for array in (a.data, b, hi, lo))):
        raise ValueError("non-finite or mismatched pressure inputs")
    dtype = np.longdouble
    ad = a.data.astype(dtype)
    hh = hi.astype(dtype)
    ll = lo.astype(dtype)
    bb = b.astype(dtype)
    residual = np.empty(n, dtype=dtype)
    row_sums = np.zeros(n, dtype=dtype)
    for i in range(n):
        start, end = a.indptr[i : i + 2]
        indices = a.indices[start:end]
        coeff = ad[start:end]
        total = np.sum(coeff, dtype=dtype)
        row_sums[i] = total
        # Explicitly group diagonal cancellation before the small correction.
        value = np.sum(coeff * (hh[indices] - hh[i]), dtype=dtype)
        value += np.sum(coeff * (ll[indices] - ll[i]), dtype=dtype)
        value += total * hh[i] + total * ll[i]
        residual[i] = bb[i] - value
    symmetric = (a - a.T).nnz == 0 and np.all(a.diagonal() > 0)
    if symmetric:
        off = a.copy()
        off.setdiag(0)
        off.eliminate_zeros()
        symmetric = np.all(off.data <= 0)
    return residual, row_sums, bool(symmetric)



def algebraic_face_balance(
    a: csr_matrix, b: np.ndarray, hi: np.ndarray, lo: np.ndarray,
    row_sums: np.ndarray,
) -> dict:
    """Independent CSR edge-flux decomposition, not geometric Rhie-Chow.

    For symmetric pressure Laplacians, the off-diagonal entries define
    equal-and-opposite internal face corrections. The row sums represent
    boundary pressure leakage; physical FV faces still need verification.
    """
    n = len(b)
    divergence = np.zeros(n, dtype=np.longdouble)
    h = hi.astype(np.longdouble)
    l = lo.astype(np.longdouble)
    a_values = a.data.astype(np.longdouble)
    for i in range(n):
        for k in range(a.indptr[i], a.indptr[i + 1]):
            j = int(a.indices[k])
            if j <= i:
                continue
            face_flux = -a_values[k] * ((h[i] - h[j]) + (l[i] - l[j]))
            divergence[i] += face_flux
            divergence[j] -= face_flux
    internal_sum = np.sum(divergence, dtype=np.longdouble)
    for i in range(n):
        divergence[i] += row_sums[i] * h[i] + row_sums[i] * l[i]
    bb = b.astype(np.longdouble)
    norm_b = np.sqrt(np.sum(bb * bb, dtype=np.longdouble))
    defect = bb - divergence
    return {
        "algebraic_face_relative_residual": float(
            np.sqrt(np.sum(defect * defect, dtype=np.longdouble)) / norm_b
        ),
        "internal_face_flux_pair_global_sum": float(internal_sum),
        "boundary_leakage_global_sum": float(np.sum(
            divergence, dtype=np.longdouble) - internal_sum
        ),
        "rhs_global_sum": float(np.sum(bb, dtype=np.longdouble)),
    }


def two_sum_update(hi: np.ndarray, lo: np.ndarray, delta: np.ndarray) -> None:
    """In-place vectorized TwoSum, maintaining an FP64 high/low expansion."""
    if not np.all(np.isfinite(delta)):
        raise ValueError("non-finite LU correction")
    s = hi + delta
    bb = s - hi
    err = (hi - (s - bb)) + (delta - bb)
    t = lo + err
    updated = s + t
    leftover = t - (updated - s)
    if not np.isfinite(updated).all() or not np.isfinite(leftover).all():
        raise ValueError("compensated pressure update overflow")
    hi[:], lo[:] = updated, leftover


def inspect(a: csr_matrix, b: np.ndarray, *, tolerance: float = 1e-11,
            refinements: int = 2, split_output: Path | None = None,
            capture_checksum: str | None = None) -> dict:
    if not (math.isfinite(tolerance) and tolerance > 0.0) or not 0 <= refinements <= 8:
        raise ValueError("invalid pressure replay controls")
    if np.finfo(np.longdouble).eps >= np.finfo(np.float64).eps:
        raise RuntimeError("extended precision is unavailable on this platform")
    a = a.tocsr(copy=True)
    a.sum_duplicates()
    b = np.asarray(b, dtype=np.float64)
    if a.shape[0] != a.shape[1] or a.shape[0] != b.size or not b.size:
        raise ValueError("invalid pressure CSR dimensions")
    if not np.all(np.isfinite(a.data)) or not np.all(np.isfinite(b)):
        raise ValueError("non-finite pressure matrix or RHS")
    b_scale = np.linalg.norm(b)
    if not np.isfinite(b_scale) or b_scale <= 0:
        raise ValueError("pressure replay requires a nonzero RHS")
    factor = splu(a.tocsc())
    hi = factor.solve(b)
    lo = np.zeros_like(hi)
    observations = []
    for iteration in range(refinements + 1):
        residual, row_sums, graph_eligible = compensated_residual(a, b, hi, lo)
        corrected = float(np.sqrt(np.sum(residual * residual, dtype=np.longdouble)) /
                          np.longdouble(b_scale))
        rounded = float(np.linalg.norm(b - a @ hi) / b_scale)
        observations.append({
            "iteration": iteration,
            "stored_fp64_relative_residual": rounded,
            "compensated_original_relative_residual": corrected,
            "compensated_algebraic_gate_passed": corrected <= tolerance,
            "stored_fp64_gate_passed": rounded <= tolerance,
        })
        if iteration == refinements:
            break
        delta = factor.solve(np.asarray(residual, dtype=np.float64))
        two_sum_update(hi, lo, delta)
    original_residual, row_sums, graph_eligible = compensated_residual(a, b, hi, lo)
    if not np.all(np.isfinite(original_residual)):
        raise ValueError("non-finite compensated pressure residual")
    # Both the internal owner/neighbour graph and boundary leakage come from
    # the original CSR. They match Ax algebraically, but this does not claim
    # that the production face fluxes use the retained low-order components.
    balance = float(np.sqrt(np.sum(original_residual**2, dtype=np.longdouble)) /
                    np.longdouble(b_scale))
    face_audit = (algebraic_face_balance(a, b, hi, lo, row_sums)
                  if graph_eligible else None)
    if split_output is not None:
        if not graph_eligible or balance > tolerance or not capture_checksum:
            raise ValueError("cannot export uncertified or unlinked pressure split")
        split_output = Path(split_output)
        if split_output.exists():
            raise FileExistsError("refusing to overwrite an existing split")
        split_output.parent.mkdir(parents=True, exist_ok=True)
        with split_output.open("x", encoding="ascii") as output:
            output.write("CFDX_PRESSURE_SPLIT_V1\n")
            output.write(f"{b.size} {capture_checksum}\n")
            for high, low in zip(hi, lo):
                output.write(f"{float(high).hex()} {float(low).hex()}\n")
    return {
        "schema": "cfdx-bfs-pressure-compensated-audit/v1",
        "status": "ALGEBRAIC_REPLAY_ONLY_NOT_FV_QUALIFIED",
        "n": int(a.shape[0]),
        "nnz": int(a.nnz),
        "original_tolerance_unchanged": tolerance,
        "longdouble_machine_epsilon": float(np.finfo(np.longdouble).eps),
        "pressure_operator_symmetric_m_matrix": graph_eligible,
        "boundary_row_sum_nonzero_count": int(np.count_nonzero(abs(row_sums) >
            1e-13 * np.asarray(abs(a).sum(axis=1)).ravel())),
        "algebraic_compensated_balance_relative": balance,
        "algebraic_face_audit": face_audit,
        "compensated_split_written": str(split_output) if split_output else None,
        "observations": observations,
        "final_split_low_norm": float(np.linalg.norm(lo)),
        "final_split_high_norm": float(np.linalg.norm(hi)),
        "original_fp64_pressure_gate_passed": observations[-1]["stored_fp64_gate_passed"],
        "physical_face_flux_certified": False,
        "sst_yplus_cf_reattachment_certified": False,
        "accept_as_production_pressure_correction": False,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--capture", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--refinements", type=int, default=2)
    parser.add_argument("--split-output", type=Path,
                        help="optional export of high/low FP64 correction for C++ audit")
    args = parser.parse_args()
    matrix, rhs, meta = read_capture(args.capture)
    checksum = args.capture.read_text(encoding="ascii").splitlines()[-1].strip()
    if len(checksum) != 16 or not all(ch in "0123456789abcdef" for ch in checksum.lower()):
        raise ValueError("capture checksum footer is missing")
    report = inspect(matrix, rhs, tolerance=meta["target_tolerance"],
                     refinements=args.refinements, split_output=args.split_output,
                     capture_checksum=checksum)
    report["capture_sha256"] = _sha(args.capture)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": report["status"], "observations": report["observations"],
                      "physical_face_flux_certified": False}, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
