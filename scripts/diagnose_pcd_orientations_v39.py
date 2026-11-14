"""V39 forensic PCD composition A/B on one production-captured saddle matrix.

The two noncommuting orders and both signs are compared against the same
true frozen Schur. This does not select or validate a production solver.
A SIMPLE diagonal-velocity and cell-local 3x3 surrogate are oracle-only
comparators (their exact dense inverse is NOT a production preconditioner).
"""
from __future__ import annotations

import argparse
import json
import os
import warnings
from pathlib import Path

import numpy as np
from scipy.linalg import LinAlgWarning, lu_factor, lu_solve

from diagnose_pcd_schur_v38 import analyze as v38_analyze
from diagnose_pcd_schur_v38 import matrix_market, normalized_gap, pressure_operator

SCHEMA = "cfdx-pcd-orientation-audit/v39"


def exact_schur(a: np.ndarray, b: np.ndarray, n: int) -> tuple[np.ndarray, np.ndarray]:
    if a.shape != (4*n, 4*n) or b.shape != (4*n, 1):
        raise ValueError("V39 coupled A/b dimensions inconsistent with pressure size")
    u, g = a[:3*n, :3*n], a[:3*n, 3*n:]
    d, c = a[3*n:, :3*n], a[3*n:, 3*n:]
    with warnings.catch_warnings():
        warnings.simplefilter("error", LinAlgWarning)
        lu, piv = lu_factor(u)
        s = c - d @ lu_solve((lu, piv), g)
        r = b[3*n:, 0] - d @ lu_solve((lu, piv), b[:3*n, 0])
    if not np.isfinite(s).all() or not np.isfinite(r).all():
        raise ValueError("nonfinite exact Schur matrix or RHS")
    return s, r


def surrogate_schur(a: np.ndarray, n: int, mode: str) -> np.ndarray:
    """Construct a diagnostic C-D Hinv G from diag or local 3x3 velocity H."""
    u, g = a[:3*n, :3*n], a[:3*n, 3*n:]
    d, c = a[3*n:, :3*n], a[3*n:, 3*n:]
    if mode == "scalar_velocity_diagonal":
        diag = np.diag(u)
        if np.any(~np.isfinite(diag)) or np.any(diag == 0):
            raise ValueError("zero/nonfinite velocity diagonal in surrogate")
        h_inv_g = g / diag[:, None]
    elif mode == "cell_local_3x3":
        h_inv_g = np.zeros_like(g)
        for cell in range(n):
            ids = [cell + comp*n for comp in range(3)]
            block = u[np.ix_(ids, ids)]
            with warnings.catch_warnings():
                warnings.simplefilter("error", LinAlgWarning)
                lu, piv = lu_factor(block)
                h_inv_g[ids, :] = lu_solve((lu, piv), g[ids, :])
    else:
        raise ValueError("unknown local velocity surrogate")
    approx = c - d @ h_inv_g
    if not np.isfinite(approx).all():
        raise ValueError("nonfinite local Schur surrogate")
    return approx


def certify_action(s: np.ndarray, rhs: np.ndarray, action: np.ndarray) -> dict:
    if not np.isfinite(action).all():
        raise ValueError("nonfinite Schur pressure action")
    gap = normalized_gap(rhs, action @ rhs, s)
    err = float(np.linalg.norm(s @ action - np.eye(s.shape[0]), ord="fro")
                / np.sqrt(s.shape[0]))
    if not np.isfinite(err):
        raise ValueError("nonfinite Frobenius action defect")
    return {"frobenius_defect_over_sqrt_n": err,
            "actual_schur_rhs": gap}


def hotspot_coordinates(mesh_file: Path, top_rows: list[dict], n: int) -> list[dict]:
    """Use the native HDF5 cell-to-face-to-vertex connectivity, if requested."""
    import h5py  # Optional for the CLI; not required for the matrix-only audit.
    with h5py.File(mesh_file, "r") as h5:
        required = ("points", "cell_offsets", "cell_faces",
                    "face_offsets", "face_vertices")
        if any(key not in h5 for key in required):
            raise ValueError("V39 mesh lacks native CSR geometry datasets")
        points = np.asarray(h5["points"][:], dtype=np.float64)
        co, cf = h5["cell_offsets"][:], h5["cell_faces"][:]
        fo, fv = h5["face_offsets"][:], h5["face_vertices"][:]
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all() or \
       len(co) != n + 1 or co[0] != 0 or co[-1] != len(cf) or \
       len(fo) < 2 or fo[0] != 0 or fo[-1] != len(fv) or np.any(np.diff(co) <= 0) or \
       np.any(np.diff(fo) <= 0) or np.any(cf >= len(fo) - 1) or \
       np.any(fv >= len(points)):
        raise ValueError("V39 mesh CSR geometry is inconsistent with pressure cells")
    result = []
    for probe in top_rows:
        i = probe["cell"]
        faces = cf[co[i]:co[i+1]]
        vertices = np.unique(np.concatenate([fv[fo[f]:fo[f+1]] for f in faces]))
        if len(vertices) < 4:
            raise ValueError("V39 mesh cell contains fewer than four vertices")
        center = np.mean(points[vertices, :], axis=0)
        result.append({**probe, "center_m": [float(v) for v in center]})
    return result


def analyze(coupled_A: Path, coupled_b: Path,
            mp_file: Path, kp_file: Path, fp_file: Path,
            max_cells: int = 256, mesh_file: Path | None = None) -> dict:
    parent = v38_analyze(coupled_A, coupled_b, mp_file, kp_file, fp_file,
                         max_cells=max_cells)
    n = parent["pressure_cells"]
    a = matrix_market(coupled_A)
    b = matrix_market(coupled_b)
    mp = pressure_operator(mp_file, max_cells)
    kp = pressure_operator(kp_file, max_cells)
    fp = pressure_operator(fp_file, max_cells)
    s, r = exact_schur(a, b, n)
    minv = np.diag(1.0 / np.diag(mp))
    kp_inv = np.linalg.pinv(kp, rcond=1e-12)
    fp_inv = np.linalg.pinv(fp, rcond=1e-12)
    inverse_actions = {
        "legacy_negative_fp_inverse": -fp_inv,
        "old_negative_Minv_Fp_Kinv": -minv @ fp @ kp_inv,
        "v39_negative_Kinv_Fp_Minv": -kp_inv @ fp @ minv,
        "v39_positive_Kinv_Fp_Minv_oracle_only": kp_inv @ fp @ minv,
    }
    for kind in ("scalar_velocity_diagonal", "cell_local_3x3"):
        approx_s = surrogate_schur(a, n, kind)
        inverse_actions[kind + "_exact_dense_inverse_oracle_only"] = np.linalg.pinv(
            approx_s, rcond=1e-12)
    report_actions = {key: certify_action(s, r, val)
                      for key, val in inverse_actions.items()}
    # Norm mismatch is an operator diagnostic, not a statement of the PCD
    # boundary condition's physical correctness on the whole domain.
    sr = np.linalg.norm(s, axis=1)
    fr = np.linalg.norm(fp, axis=1)
    idx = np.argsort(sr / np.maximum(fr, 1e-300))[-5:][::-1]
    sample = [{"cell": int(i), "schur_row_l2": float(sr[i]),
               "Fp_row_l2": float(fr[i]),
               "ratio_schur_over_fp": float(sr[i] / max(fr[i], 1e-300))}
              for i in idx]
    best = min(report_actions, key=lambda key:
               report_actions[key]["frobenius_defect_over_sqrt_n"])
    model_only = ["legacy_negative_fp_inverse", "old_negative_Minv_Fp_Kinv",
                  "v39_negative_Kinv_Fp_Minv"]
    best_negative = min(model_only, key=lambda key:
                        report_actions[key]["frobenius_defect_over_sqrt_n"])
    if mesh_file is not None:
        sample = hotspot_coordinates(mesh_file, sample, n)
    report = {
        "schema": SCHEMA,
        "input_sha256": {k: v["sha256"] for k,v in parent["inputs"].items()},
        "cells": n,
        "spectral": parent["spectral"],
        "pressure_schur_frozen_only": True,
        "full_nonlinear_U_p_k_omega_jacobian": False,
        "pressure_nullity_mismatch": parent["spectral"]["Fp"]["numerical_nullity"] !=
                                     parent["spectral"]["true_saddle_S"]["numerical_nullity"],
        "actions": report_actions,
        "largest_S_vs_Fp_row_norm_ratios": sample,
        "mesh_geometry_used": mesh_file is not None,
        "best_observed_frobenius_model_not_recommended": best,
        "best_tested_negative_pcd_model_not_recommended": best_negative,
        "selected_for_production": False,
        "physical_converged": False,
        "constraints": {
            "both_PCD_operator_orders_noncommute": True,
            "true_S_computed_dense_in_a_bounded_diagnostic": True,
            "no_silent_solver_or_sign_choice": True,
            "no_relaxed_linear_or_physical_tolerances": True,
        },
    }
    json.dumps(report, allow_nan=False)
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coupled-A", type=Path, required=True)
    parser.add_argument("--coupled-b", type=Path, required=True)
    parser.add_argument("--mp", type=Path, required=True)
    parser.add_argument("--kp", type=Path, required=True)
    parser.add_argument("--fp", type=Path, required=True)
    parser.add_argument("--max-cells", type=int, default=256)
    parser.add_argument("--mesh", type=Path,
                        help="optional exact native HDF5 mesh for Schur-row cell centers")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--maximum-v39-action-defect", type=float,
                        help="optional strict ALGORITHMIC gate; never a physical gate")
    args = parser.parse_args(argv)
    try:
        report = analyze(args.coupled_A, args.coupled_b,
                         args.mp, args.kp, args.fp, max_cells=args.max_cells,
                         mesh_file=args.mesh)
        v39err = report["actions"]["v39_negative_Kinv_Fp_Minv"]["frobenius_defect_over_sqrt_n"]
        limit = args.maximum_v39_action_defect
        if limit is not None and (not np.isfinite(limit) or limit <= 0):
            raise ValueError("maximum-v39-action-defect must be positive and finite")
        failed = limit is not None and v39err > limit
        report["algebraic_gate"] = {
            "requested": limit is not None,
            "maximum": limit,
            "pass": not failed,
            "physical_gate": False,
        }
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temp = args.output.with_name(args.output.name + ".tmp")
        temp.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        os.replace(temp, args.output)
        print("PCD_V39_ORIENTATION cells={} old={:.6g} reverse={:.6g} "
              "gate={} physical_qualified=0".format(
                  report["cells"],
                  report["actions"]["old_negative_Minv_Fp_Kinv"]["frobenius_defect_over_sqrt_n"],
                  v39err, "FAIL" if failed else "PASS"))
        return 2 if failed else 0
    except (ValueError, OSError, RuntimeError, FloatingPointError,
            np.linalg.LinAlgError) as error:
        print("PCD_V39_ORIENTATION_ERROR: " + str(error),
              file=__import__("sys").stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
