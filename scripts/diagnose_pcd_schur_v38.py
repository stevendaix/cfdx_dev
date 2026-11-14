"""V38: audit the PCD pressure Schur models on ONE real captured saddle A.

The exact finite-volume coupled matrix A and its RHS must come from the same
production state as Mp/Kp/Fp. The CFP01 pressure operator files contain ZERO
placeholder RHS/initial vectors. An independent capture of an actual Fp RHS
is NOT a pressure operator file. No physical convergence is certified here.

This bounded-size forensic audit compares *ideal algebraic inverse actions*,
not finite-budget GMRES/FGMRES actions. Any discrepancy is model error even
if the inner solvers could be made exact. Computational cost is O(n^3).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import warnings
from pathlib import Path

import numpy as np
from scipy.io import mmread
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve
from certify_pcd_fp_bi_nullspace_v13 import read_capture

SCHEMA = "cfdx-pcd-exact-schur-audit/v38"


def digest(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def matrix_market(path: Path) -> np.ndarray:
    raw = mmread(path)
    if hasattr(raw, "toarray"):
        raw = raw.toarray()
    a = np.asarray(raw, dtype=np.float64)
    if a.ndim != 2 or not np.isfinite(a).all():
        raise ValueError(f"invalid/nonfinite MatrixMarket object: {path}")
    return a


def pressure_operator(path: Path, nmax: int) -> np.ndarray:
    a, rhs, x0, _ = read_capture(path, max_rows=nmax)
    if np.any(rhs != 0.0) or np.any(x0 != 0.0):
        raise ValueError("V38 pressure-operator capture requires ZERO placeholder vectors")
    return a.toarray()


def spectral_info(a: np.ndarray, rtol: float = 1e-12) -> dict:
    values = np.linalg.svd(a, compute_uv=False)
    if values.size == 0 or not np.all(np.isfinite(values)) or values[0] <= 0.0:
        raise ValueError("invalid/zero pressure spectral operator")
    rank = int(np.count_nonzero(values > rtol * values[0]))
    last_nonzero = float(values[rank - 1]) if rank else 0.0
    return {
        "rows": int(a.shape[0]),
        "numerical_rank": rank,
        "numerical_nullity": int(a.shape[0] - rank),
        "sigma_max": float(values[0]),
        "sigma_smallest_retained": last_nonzero,
        "sigma_min": float(values[-1]),
        "condition_retained": float(values[0] / last_nonzero),
    }


def normalized_gap(r: np.ndarray, action: np.ndarray, schur: np.ndarray) -> dict:
    norm_r = np.linalg.norm(r)
    if norm_r < 1e-300:
        raise ValueError("zero tested Schur right-hand side")
    im = schur @ action
    denom = float(im @ im)
    alpha = float((r @ im) / denom) if denom > 1e-300 else 0.0
    output = {
        "schur_true_residual_relative": float(np.linalg.norm(im - r) / norm_r),
        "best_scalar_alpha": alpha,
        "schur_true_residual_after_best_scalar": float(
            np.linalg.norm(alpha * im - r) / norm_r),
        "action_amplification_l2": float(np.linalg.norm(action) / norm_r),
    }
    if not all(np.isfinite(v) for v in output.values()):
        raise ValueError("nonfinite PCD pressure action certificate")
    return output


def analyze(coupled_path: Path, coupled_rhs_path: Path,
            mp_path: Path, kp_path: Path, fp_path: Path,
            max_cells: int = 256) -> dict:
    if not 1 <= max_cells <= 1024:
        raise ValueError("V38 audit max_cells must be within 1..1024")
    mp = pressure_operator(mp_path, max_cells)
    kp = pressure_operator(kp_path, max_cells)
    fp = pressure_operator(fp_path, max_cells)
    n = mp.shape[0]
    if n > max_cells or any(a.shape != (n, n) for a in (kp, fp)):
        raise ValueError("PCD Mp/Kp/Fp pressure dimensions mismatch")
    if not np.allclose(mp, np.diag(np.diag(mp)), atol=0.0, rtol=0.0):
        raise ValueError("V38 mass action requires exact diagonal pressure mass")
    mass_diag = np.diag(mp)
    if not np.all(np.isfinite(mass_diag)) or np.any(mass_diag <= 0):
        raise ValueError("PCD mass must be positive and finite")
    A = matrix_market(coupled_path)
    b = matrix_market(coupled_rhs_path)
    if A.shape != (4 * n, 4 * n) or b.shape != (4 * n, 1):
        raise ValueError("V38 expects exact 3-component U + 1 pressure block A and b")
    a_uu = A[:3*n, :3*n]
    G = A[:3*n, 3*n:]
    D = A[3*n:, :3*n]
    C = A[3*n:, 3*n:]
    # This is the TRUE frozen saddle Schur complement at the captured state.
    # It is not a reassembled nonlinear Jacobian or the FV physical residual.
    with warnings.catch_warnings():
        warnings.simplefilter("error", LinAlgWarning)
        lu, piv = lu_factor(a_uu, check_finite=True)
        schur = C - D @ lu_solve((lu, piv), G)
        schur_rhs = b[3*n:, 0] - D @ lu_solve((lu, piv), b[:3*n, 0])
    if not np.isfinite(schur).all() or not np.isfinite(schur_rhs).all():
        raise ValueError("nonfinite coupled Schur complement")
    fp_inverse = np.linalg.pinv(fp, rcond=1e-12)
    kp_inverse = np.linalg.pinv(kp, rcond=1e-12)
    schur_inverse = np.linalg.pinv(schur, rcond=1e-12)
    # ACTUAL CFDX PcdSchurApproximation::apply formulae.
    legacy = -fp_inverse
    commutator = -(fp @ kp_inverse) / mass_diag[:, None]
    # Small-system ORACLE ONLY: exact diagonal of the true saddle S. The
    # production SIMPLE/BlockLocal methods do not have this exact diagonal.
    diag_s = np.diag(schur)
    valid_diag = bool(np.min(np.abs(diag_s)) >
                      1e-14 * max(1.0, np.max(np.abs(diag_s))))
    oracle_diag = np.diag(1.0 / diag_s) if valid_diag else None
    tests = {"actual_coupled_schur_rhs": schur_rhs}
    for idx in sorted({0, n // 4, n // 2, n - 1}):
        vec = np.zeros(n)
        vec[idx] = 1.0
        tests[f"unit_cell_{idx}"] = vec
    # Fixed, deterministic probing vector, independent of the RHS.
    tests["seeded_random_20261111"] = np.random.default_rng(20261111).standard_normal(n)
    results = {}
    for name, r in tests.items():
        results[name] = {
            "legacy_negative_fp_inverse": normalized_gap(r, legacy @ r, schur),
            "pcd_lumped_commutator": normalized_gap(r, commutator @ r, schur),
            "exact_schur_pseudoinverse_reference": normalized_gap(r, schur_inverse @ r, schur),
        }
        if oracle_diag is not None:
            results[name]["oracle_only_exact_schur_diagonal"] = normalized_gap(
                r, oracle_diag @ r, schur)
    eye = np.eye(n)
    overall = {
        "legacy_frobenius_action_error_over_sqrt_n": float(
            np.linalg.norm(schur @ legacy - eye, ord="fro") / np.sqrt(n)),
        "commutator_frobenius_action_error_over_sqrt_n": float(
            np.linalg.norm(schur @ commutator - eye, ord="fro") / np.sqrt(n)),
        "exact_pseudoinverse_frobenius_action_error_over_sqrt_n": float(
            np.linalg.norm(schur @ schur_inverse - eye, ord="fro") / np.sqrt(n)),
        "oracle_only_exact_schur_diagonal_available": valid_diag,
    }
    if oracle_diag is not None:
        overall["oracle_only_exact_schur_diagonal_frobenius_action_error_over_sqrt_n"] = float(
            np.linalg.norm(schur @ oracle_diag - eye, ord="fro") / np.sqrt(n))
    report = {
        "schema": SCHEMA,
        "inputs": {k: {"path": str(v), "sha256": digest(v)}
                   for k, v in (("coupled_A", coupled_path), ("coupled_b", coupled_rhs_path),
                                ("Mp", mp_path), ("Kp", kp_path), ("Fp", fp_path))},
        "pressure_cells": n,
        "saddle_matrix_n": int(A.shape[0]),
        "theory": {
            "exact_schur": "C - D Auu^{-1} G",
            "actual_rhs": "bp - D Auu^{-1} bu",
            "historical_cfdx_pcd_action": "-Fp_pseudoinverse * rhs",
            "optional_cfdx_pcd_commutator": "-diag(Mp)^{-1} Fp Kp_pseudoinverse * rhs",
            "best_scalar": "argmin_alpha ||r - alpha S (approximate action r)||_2",
            "oracle_only_exact_schur_diagonal": "diag(C-D Auu^{-1} G)^{-1}; NOT production SIMPLE",
            "warning": "pseudoinverse simulates ideal inversion, not nested Krylov; NOT physical validation",
        },
        "spectral": {"Fp": spectral_info(fp), "Kp": spectral_info(kp),
                     "true_saddle_S": spectral_info(schur)},
        "schur_action_error": overall,
        "probes": results,
        "verified_physical_convergence": False,
        "operator_is_full_nonlinear_jacobian": False,
    }
    json.dumps(report, allow_nan=False)
    return report


def qualification_gates(report: dict, *, max_action_error: float | None = None,
                        require_fp_nullity_match: bool = False) -> list[str]:
    """Optional fail-closed Schur-model gates. They are NEVER physics gates."""
    failures = []
    if max_action_error is not None:
        if not np.isfinite(max_action_error) or max_action_error <= 0:
            raise ValueError("max_action_error must be positive and finite")
        for model in ("legacy", "commutator"):
            value = report["schur_action_error"][
                model + "_frobenius_action_error_over_sqrt_n"]
            if not np.isfinite(value) or value > max_action_error:
                failures.append(f"{model} action error {value:.7g} > {max_action_error:.7g}")
    if require_fp_nullity_match:
        fp = report["spectral"]["Fp"]["numerical_nullity"]
        exact = report["spectral"]["true_saddle_S"]["numerical_nullity"]
        if fp != exact:
            failures.append(f"Fp nullity {fp} differs from exact Schur nullity {exact}")
    return failures



def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coupled-A", type=Path, required=True)
    parser.add_argument("--coupled-b", type=Path, required=True)
    parser.add_argument("--mp", type=Path, required=True)
    parser.add_argument("--kp", type=Path, required=True)
    parser.add_argument("--fp", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-cells", type=int, default=256)
    parser.add_argument("--max-action-error", type=float,
                        help="optional strict algebraic model gate, NOT a convergence tolerance")
    parser.add_argument("--require-fp-nullity-match", action="store_true",
                        help="reject PCD Fp rank deficiency absent in the true pressure Schur")
    args = parser.parse_args(argv)
    try:
        report = analyze(args.coupled_A, args.coupled_b, args.mp, args.kp, args.fp,
                         max_cells=args.max_cells)
        failures = qualification_gates(
            report, max_action_error=args.max_action_error,
            require_fp_nullity_match=args.require_fp_nullity_match)
        report["algebraic_qualification_gate"] = {
            "requested": args.max_action_error is not None or args.require_fp_nullity_match,
            "pass": len(failures) == 0,
            "failure_reasons": failures,
            "is_physical_qualification": False,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.output.with_name(args.output.name + ".tmp")
        tmp.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        os.replace(tmp, args.output)
        print("PCD_SCHUR_V38 pressure_cells={} legacy_error={:.6e} "
              "commutator_error={:.6e} algebraic_gate={} physical_qualified=0".format(
                  report["pressure_cells"],
                  report["schur_action_error"]["legacy_frobenius_action_error_over_sqrt_n"],
                  report["schur_action_error"]["commutator_frobenius_action_error_over_sqrt_n"],
                  "FAIL" if failures else "PASS"))
        for failure in failures:
            print("PCD_SCHUR_V38_GATE_FAIL: " + failure, file=__import__("sys").stderr)
        return 2 if failures else 0
    except (ValueError, OSError, RuntimeError, FloatingPointError, np.linalg.LinAlgError) as exc:
        print(f"PCD_SCHUR_V38_ERROR: {exc}", file=__import__("sys").stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
