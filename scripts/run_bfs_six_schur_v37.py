"""Run six BFS Schur approximations on the SAME native HDF5 physical case.

V37 diagnostic: same mesh, BCs, ICs, physical parameters, spatial scheme,
FGMRES and coupled_block_schur preconditioner family. Only selected Schur
approximation differs. This is NOT Driver-Seegmiller physical validation:
no NASA comparison or completed nonlinear convergence is implied by a run.

Example (small, explicitly bounded debug experiment):
  python scripts/run_bfs_six_schur_v37.py \
    --solver build/cfdx_production_solver --out build/bfs_v37 \
    --iterations 1 --upstream 4 --downstream 12 --ny 10 \
    --step-cells 2 --wall-timeout 15

Real reference input:
  python scripts/run_bfs_six_schur_v37.py \
    --solver build/cfdx_production_solver --out build/bfs_reference_v37 \
    --reference-case bfsv4.cfdx.h5 --iterations 1500
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import math
import re
from pathlib import Path
import subprocess
import sys
import time

import h5py
import numpy as np

from generate_th03_bfs_case import (
    SCHUR_MODELS, reconfigure_bfs_reference_case, write_case,
)
from diagnose_bfs_failure_v10 import classify_log
from bfs_schur_fingerprint_v40 import parse_fingerprint,compare_six_fingerprints
from qualify_bfs_schur_evidence import validate_bfs_evidence

PRECONDITIONER = "coupled_block_schur"
SCHEMA = "CFDX-BFS-SIX-SCHUR-COMMON-INPUT/V37"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def physical_input_identity(path: Path) -> dict:
    """Hash the physical input including actual datasets, not just mesh attrs.

    Numerical selections of Schur/PC are deliberately excluded from the
    common physical signature. Everything else (BCs, fields, numerics and
    raw settings) must match. Return selected Schur and preconditioner too,
    so a modified or silently substituted case is rejected.
    """
    path = Path(path)
    with h5py.File(path, "r") as h5:
        setup = json.loads(h5.attrs["case_setup_json"])
        selected = {
            x["family"]: x["configuration_key"]
            for x in setup["numerics"]["selection"]["entries"]
        }
        entries = setup["numerics"]["selection"]["entries"]
        if len(selected) != len(entries):
            raise ValueError("duplicate numerical family in BFS Schur case")
        expected_pc = f"preconditioner.{PRECONDITIONER}"
        if (selected.get("pressure_velocity") != "pressure_velocity.coupled" or
            selected.get("preconditioner") != expected_pc or
            selected.get("linear_solver") != "linear.fgmres"):
            raise ValueError("BFS six-Schur case violates the fixed solver policy")
        schur = selected.get("schur", "").removeprefix("schur.")
        if schur not in SCHUR_MODELS or selected.get("schur") != f"schur.{schur}":
            raise ValueError("BFS six-Schur case does not select an explicit Schur")
        physical_setup = copy.deepcopy(setup)
        physical_setup["numerics"]["selection"]["entries"] = [
            x for x in entries if x["family"] != "schur"
        ]
        # Exclude only the Schur field from the signature. The PC is identical
        # in all variants and stays in the signed JSON for verification.
        dataset_hashes = {}
        for name, obj in h5.items():
            if not isinstance(obj, h5py.Dataset):
                raise ValueError("unexpected HDF5 group in BFS case")
            arr = np.asarray(obj[()])
            data = hashlib.sha256()
            data.update(str(arr.dtype).encode("ascii"))
            data.update(json.dumps(list(arr.shape)).encode("ascii"))
            data.update(arr.tobytes(order="C"))
            dataset_hashes[name] = data.hexdigest()
        attrs = {}
        for key, value in h5.attrs.items():
            if key == "case_setup_json":
                continue
            if isinstance(value, np.generic):
                value = value.item()
            if not isinstance(value, (str, int, float, bool)):
                raise ValueError(f"unsupported HDF5 attribute type: {key}")
            attrs[str(key)] = value
        common = {
            "attributes": attrs,
            "datasets": dataset_hashes,
            "case_setup_except_schur": physical_setup,
        }
        data = json.dumps(common, sort_keys=True, allow_nan=False,
                          separators=(",", ":")).encode("utf-8")
        return {
            "schur": schur,
            "preconditioner": PRECONDITIONER,
            "physical_sha256": hashlib.sha256(data).hexdigest(),
            "case_sha256": sha256(path),
            "n_cells": int(h5.attrs["n_cells"]),
            "n_faces": int(h5.attrs["n_faces"]),
            "topology_hash": str(h5.attrs["topology_hash"]),
            "geometry_hash": str(h5.attrs["geometry_hash"]),
        }


def prepare_cases(out: Path, *, iterations: int, reference_case: Path | None,
                  upstream: int, downstream: int, ny: int,
                  step_cells: int, wall_stretch: float,
                  capture_operator: bool = False,
                  fingerprint_v40: bool = False) -> dict:
    """Create six independent cases; abort on any physical-input mismatch."""
    out = Path(out)
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("V37 refuses a nonempty output directory (stale evidence)")
    out.mkdir(parents=True, exist_ok=True)
    if reference_case is None:
        reference_case = out / "reference_debug.cfdx.h5"
        write_case(reference_case, n_up=upstream, n_down=downstream, ny=ny,
                   step_cells=step_cells, wall_stretch=wall_stretch,
                   algorithm="simple", iterations=iterations)
        source_kind = "GENERATED_UNIFORM_INLET_DEBUG_NOT_NASA"
    else:
        reference_case = Path(reference_case).resolve()
        source_kind = "EXTERNAL_REFERENCE_UNQUALIFIED"
    reference_sha = sha256(reference_case)
    identities = {}
    # Identical diagnostic capture path for all six variants, so the *same*
    # original saddle A,b,x0 can be hashed and compared after every run.
    # This forensic export is explicitly limited by V16 to 6000 unknowns.
    capture_prefix = out / "coupled_A_v16" if capture_operator else None
    for schur in SCHUR_MODELS:
        case = out / "cases" / f"bfs_{schur}.cfdx.h5"
        reconfigure_bfs_reference_case(
            reference_case, case, algorithm="coupled", schur=schur,
            preconditioner=PRECONDITIONER, iterations=iterations,
        )
        if capture_operator or fingerprint_v40:
            with h5py.File(case, "r+") as h5:
                if capture_operator and int(h5.attrs["n_cells"]) > 1500:
                    raise ValueError("V16 forensic capture supports at most 1500 BFS cells")
                setup = json.loads(h5.attrs["case_setup_json"])
                if capture_operator:
                    setup["numerics"]["raw_settings"]["coupled_diagnostic_capture_prefix_v16"] = str(capture_prefix)
                if fingerprint_v40:
                    setup["numerics"]["raw_settings"]["coupled_matrix_fingerprint_v40"] = "true"
                h5.attrs["case_setup_json"] = json.dumps(setup, allow_nan=False)
        identity = physical_input_identity(case)
        if identity["schur"] != schur:
            raise ValueError(f"requested Schur not present in case: {schur}")
        identities[schur] = {"case": str(case), **identity}
    if len({i["physical_sha256"] for i in identities.values()}) != 1:
        raise ValueError("BFS Schur variants DO NOT share the same physical input")
    return {
        "format": SCHEMA,
        "reference_case": str(reference_case),
        "reference_sha256": reference_sha,
        "source_kind": source_kind,
        "iterations": iterations,
        "preconditioner_policy": PRECONDITIONER,
        "operator_capture_requested": capture_operator,
        "operator_fingerprint_v40_requested": fingerprint_v40,
        "operator_capture_prefix": str(capture_prefix) if capture_prefix else None,
        "identical_physical_inputs": True,
        "common_physical_sha256": next(iter(identities.values()))["physical_sha256"],
        "cases": identities,
        "physical_validation": "NOT_EVALUATED_NO_NASA_COMPARISON",
    }


def summarize_pcd_nested_fp(text: str) -> dict:
    """Describe inexact PCD preconditioner actions without calling them physics.

    Fp is a *preconditioner*; MAX_ITER may be usable for FGMRES while
    its original, unswapped true residual remains an indispensable diagnostic.
    """
    values = []
    max_iteration = 0
    accepted = 0
    nonfinite = 0
    for line in text.splitlines():
        if not line.startswith("COUPLED_NESTED_PCD_FP "):
            continue
        pairs = dict(re.findall(r"([a-zA-Z_]+)=([^ ]+)", line))
        try:
            residual = float(pairs["original_true_relative_rhs"])
        except (KeyError, ValueError):
            nonfinite += 1
            continue
        if not math.isfinite(residual) or residual < 0:
            nonfinite += 1
            continue
        values.append(residual)
        max_iteration += int(pairs.get("status", "-1") == "1")
        accepted += int(pairs.get("accepted", "0") == "1")
    values.sort()
    if not values:
        return {"calls": 0, "malformed_or_nonfinite": nonfinite,
                "status_max_iter": 0, "accepted_actions": 0,
                "original_true_relative_p50": None,
                "original_true_relative_p90": None,
                "original_true_relative_max": None,
                "original_true_relative_over_1e3": 0,
                "original_true_relative_under_1e8": 0,
                "physical_qualified": False}
    return {
        "calls": len(values),
        "malformed_or_nonfinite": nonfinite,
        "status_max_iter": max_iteration,
        "accepted_actions": accepted,
        "original_true_relative_p50": values[len(values) // 2],
        "original_true_relative_p90": values[int(0.9 * (len(values) - 1))],
        "original_true_relative_max": values[-1],
        "original_true_relative_over_1e3": sum(x > 1e-3 for x in values),
        "original_true_relative_under_1e8": sum(x < 1e-8 for x in values),
        "physical_qualified": False,
    }


def run_one(solver: Path, out: Path, schur: str,
            case: dict, wall_timeout: float | None,
            capture_prefix: Path | None = None,
            pcd_inner_trace: bool = False) -> dict:
    case_path = Path(case["case"])
    run_dir = out / "runs" / schur
    run_dir.mkdir(parents=True, exist_ok=False)
    log_file = out / "logs" / f"bfs_{schur}.log"
    log_file.parent.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    for key in tuple(env):
        if key.startswith("CFDX_NACA_") or key.startswith("CFDX_PCD_FP_CAPTURE_"):
            env.pop(key)
    env["CFDX_BFS_STAGE_TRACE"] = "1"
    if schur == "pcd" and pcd_inner_trace:
        env["CFDX_BFS_PCD_INNER_TRACE"] = "1"
    else:
        env.pop("CFDX_BFS_PCD_INNER_TRACE", None)
    cmd = [str(solver), "--mesh", str(case_path),
           "--output-dir", str(run_dir), "--iterations", str(case["iterations"]),
           "--experimental-sst", "--benchmark-export"]
    result = {"schur": schur, "requested_pc": PRECONDITIONER,
              "exit_code": None, "execution": "NOT_STARTED",
              "command": cmd, "log": str(log_file),
              "run_dir": str(run_dir), "physical_validation": "NOT_QUALIFIED"}
    start = time.monotonic()
    with log_file.open("w", encoding="utf-8") as log:
        try:
            completed = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT,
                                       timeout=wall_timeout, env=env, check=False)
            result["exit_code"] = completed.returncode
            result["execution"] = "EXIT_ZERO" if completed.returncode == 0 else "FAILED"
        except subprocess.TimeoutExpired:
            result["execution"] = "TIMEOUT_INCOMPLETE"
        except OSError as error:
            result["execution"] = "FAILED_TO_START"
            result["error"] = str(error)
    result["wall_seconds"] = round(time.monotonic() - start, 3)
    result["log_sha256"] = sha256(log_file)
    if capture_prefix is not None:
        capture_dir = out / "operator_captures"
        capture_dir.mkdir(parents=True, exist_ok=True)
        captured = {}
        for key in ("A", "b", "x0"):
            source = Path(f"{capture_prefix}_{key}.mtx")
            if source.is_file():
                target = capture_dir / f"{schur}_{key}.mtx"
                source.rename(target)
                captured[key] = {"path": str(target), "sha256": sha256(target)}
        result["operator_capture"] = captured
        result["full_operator_captured"] = len(captured) == 3
    log_text = log_file.read_text(encoding="utf-8", errors="replace")
    result["operator_fingerprint_v40"] = parse_fingerprint(log_text)
    result["first_blocker"] = classify_log(
        log_text, execution=result["execution"], exit_code=result["exit_code"])
    if schur == "pcd" and pcd_inner_trace:
        result["nested_fp"] = summarize_pcd_nested_fp(log_text)
    result["v8_evidence"] = validate_bfs_evidence(
        run_dir, result["exit_code"], algorithm="coupled", schur=schur,
        preconditioner=PRECONDITIONER, case_sha256=case["case_sha256"],
        case_path=case_path)
    result["nonlinear_convergence_proven"] = result["v8_evidence"]["nonlinear_v8_proven"]
    return result


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--solver", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--reference-case", type=Path)
    parser.add_argument("--iterations", type=int, default=1500)
    parser.add_argument("--wall-timeout", type=float,
                        help="explicit per-run resource cap; incomplete runs remain unqualified")
    parser.add_argument("--upstream", type=int, default=4)
    parser.add_argument("--downstream", type=int, default=16)
    parser.add_argument("--ny", type=int, default=10)
    parser.add_argument("--step-cells", type=int, default=2)
    parser.add_argument("--wall-stretch", type=float, default=0.0)
    parser.add_argument("--pcd-inner-trace", action="store_true",
                        help="trace PCD Fp approximate actions; diagnostic only")
    parser.add_argument("--fingerprint-operator-v40", action="store_true",
                        help="low-memory FNV64 A,b,x0 fingerprint on all six; required for operator-equality claim")
    parser.add_argument("--capture-original-operator", action="store_true",
                        help="forensic V16 exact A,b,x0 capture; <=1500 cells; default OFF")
    args = parser.parse_args(argv)
    if not args.solver.is_file() or not os.access(args.solver, os.X_OK):
        parser.error("production solver binary not found or not executable")
    if not (1 <= args.iterations <= 100000):
        parser.error("iteration count must be in [1,100000]")
    if args.wall_timeout is not None and args.wall_timeout <= 0:
        parser.error("wall timeout must be positive")
    solver = args.solver.resolve()
    out = args.out.resolve()
    campaign = prepare_cases(
        out, iterations=args.iterations,
        reference_case=args.reference_case, upstream=args.upstream,
        downstream=args.downstream, ny=args.ny, step_cells=args.step_cells,
        wall_stretch=args.wall_stretch,
        capture_operator=args.capture_original_operator,
        fingerprint_v40=args.fingerprint_operator_v40)
    campaign["solver_sha256"] = sha256(solver)
    campaign["wall_timeout"] = args.wall_timeout
    campaign["results"] = []
    target = out / "six_schur_v37.json"
    target.write_text(json.dumps(campaign, sort_keys=True, indent=2,
                                 allow_nan=False) + "\n", encoding="utf-8")
    for schur in SCHUR_MODELS:
        print(f"V37_BFS_SCHUR_BEGIN {schur}", flush=True)
        result = run_one(solver, out, schur,
                         {**campaign["cases"][schur], "iterations": args.iterations},
                         args.wall_timeout,
                         Path(campaign["operator_capture_prefix"])
                         if args.capture_original_operator else None,
                         pcd_inner_trace=args.pcd_inner_trace)
        campaign["results"].append(result)
        target.write_text(json.dumps(campaign, sort_keys=True, indent=2,
                                     allow_nan=False) + "\n", encoding="utf-8")
        print(f"V37_BFS_SCHUR_END {schur} {result['execution']} "
              f"{result['first_blocker']['first_blocker']}", flush=True)
    proved = all(r["nonlinear_convergence_proven"] for r in campaign["results"])
    campaign["all_six_nonlinear_converged"] = proved
    if args.capture_original_operator:
        captures = [r.get("operator_capture", {}) for r in campaign["results"]]
        campaign["all_six_original_operators_captured"] = all(
            {"A", "b", "x0"} <= set(x) for x in captures)
        campaign["all_six_identical_original_operator"] = (
            campaign["all_six_original_operators_captured"] and
            all(len({c[k]["sha256"] for c in captures}) == 1
                for k in ("A", "b", "x0")))
    else:
        campaign["all_six_identical_original_operator"] = None
    if args.fingerprint_operator_v40:
        campaign["operator_fingerprint_v40_comparison"] = compare_six_fingerprints(
            campaign["results"], SCHUR_MODELS)
    target.write_text(json.dumps(campaign, sort_keys=True, indent=2,
                                 allow_nan=False) + "\n", encoding="utf-8")
    return 0 if (not args.fingerprint_operator_v40 or
                 campaign["operator_fingerprint_v40_comparison"]["all_six_same_A_b_x0_fingerprint"]) and proved and (not args.capture_original_operator or
                           campaign["all_six_identical_original_operator"]) else 2


if __name__ == "__main__":
    raise SystemExit(main())
