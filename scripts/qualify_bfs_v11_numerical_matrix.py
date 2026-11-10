"""BFS V11 controlled numerical A/B (NOT an experimental NASA certification).

Every variant is cloned from the *same* explicit second-order CFDX HDF5 case.
Geometry, initial conditions, material, turbulence BCs and source are kept
byte-identical. Changes concern only selected solver, pressure precision and
SST nonlinear method. All targets and V8/U,p,k,omega gates remain unchanged.

The parent process kills an entire process group at a diagnostic deadline;
no dangling solve from a timed-out long SST/PCD case can pollute later runs.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

import h5py

from diagnose_bfs_failure_v10 import classify_log
from generate_th03_bfs_case import reconfigure_bfs_reference_case
from qualify_bfs_schur_evidence import validate_bfs_evidence

PROFILES = {
    "simple_baseline": ("simple", None, None, {}),
    "simple_full_face_hilo": (
        "simple",
        None,
        None,
        {
            "v18_sparse_flux_pressure_jacobian": "true",
            "v19_full_pressure_flux_step": "true",
            "pressure_high_low_production": "true",
        },
    ),
    "simple_full_face_hilo_sst_newton": (
        "simple",
        None,
        None,
        {
            "v18_sparse_flux_pressure_jacobian": "true",
            "v19_full_pressure_flux_step": "true",
            "pressure_high_low_production": "true",
            "sst_omega_newton_sink": "true",
            "sst_sequential_omega_first": "true",
        },
    ),
    "block_local_baseline": ("coupled", "block_local", "coupled_block_schur", {}),
    "block_local_sst_newton": (
        "coupled",
        "block_local",
        "coupled_block_schur",
        {"sst_omega_newton_sink": "true", "sst_sequential_omega_first": "true"},
    ),
    "pcd_inner_trace": ("coupled", "pcd", "auto", {}),
    "pcd_fp_row_equilibrated": (
        "coupled",
        "pcd",
        "auto",
        {"pcd_fp_row_equilibration_experimental": "true"},
    ),
    "pcd_fp_quality_guard": (
        "coupled",
        "pcd",
        "auto",
        {"pcd_fp_quality_guard_experimental": "true"},
    ),
    "pcd_fp_second_pass": (
        "coupled",
        "pcd",
        "auto",
        {"pcd_fp_second_pass_experimental": "true"},
    ),
    "pcd_fp_recirculation_anchor": (
        "coupled",
        "pcd",
        "auto",
        {"pcd_fp_recirculation_anchor_experimental": "true"},
    ),
    "pcd_fp_inflow_and_recirculation_anchor": (
        "coupled",
        "pcd",
        "auto",
        {
            "pcd_fp_recirculation_anchor_experimental": "true",
            "pcd_pressure_inflow_dirichlet_experimental": "1",
        },
    ),
    "pcd_lumped_mass_commutator": (
        "coupled",
        "pcd",
        "auto",
        {"pcd_lumped_mass_commutator_experimental": "true"},
    ),
}


def sha(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda: f.read(2**20), b""):
            h.update(b)
    return h.hexdigest()


def case_physics_signature(path: Path) -> str:
    """Hash physical setup and all mesh datasets, excluding numerical controls."""
    dig = hashlib.sha256()
    with h5py.File(path) as f:
        setup = json.loads(f.attrs["case_setup_json"])
        setup.pop("numerics", None)
        dig.update(json.dumps(setup, sort_keys=True, allow_nan=False).encode())
        for key in sorted(f.keys()):
            d = f[key]
            dig.update(key.encode())
            dig.update(str(d.dtype).encode())
            dig.update(str(d.shape).encode())
            dig.update(d[()].tobytes())
        for key in ("topology_hash", "geometry_hash", "n_cells", "n_faces"):
            dig.update(str(f.attrs[key]).encode())
    return dig.hexdigest()


def prepare_case(
    reference: Path, case: Path, profile: str, iterations: int, *, sst_inner_budget: int
) -> dict:
    if profile not in PROFILES:
        raise ValueError(f"unknown profile {profile!r}")
    if sst_inner_budget < 1:
        raise ValueError("positive SST budget required")
    algorithm, schur, pc, overrides = PROFILES[profile]
    record = reconfigure_bfs_reference_case(
        reference,
        case,
        algorithm=algorithm,
        schur=schur,
        preconditioner=pc,
        iterations=iterations,
    )
    with h5py.File(case, "r+") as f:
        setup = json.loads(f.attrs["case_setup_json"])
        numerics = setup["numerics"]
        raw = numerics.setdefault("raw_settings", {})
        raw["sst_transport_max_iterations"] = str(sst_inner_budget)
        raw.update(overrides)
        f.attrs["case_setup_json"] = json.dumps(setup, allow_nan=False)
    reference_signature = case_physics_signature(reference)
    output_signature = case_physics_signature(case)
    if output_signature != reference_signature:
        raise RuntimeError("the numerical profile changed geometry or physical setup")
    record.update(
        {
            "profile": profile,
            "case_sha256": sha(case),
            "physics_signature_sha256": output_signature,
            "sst_inner_budget": sst_inner_budget,
            "numerical_overrides": overrides,
        }
    )
    return record


def nested_pcd_summary(log: str) -> dict:
    """Aggregate *original-operator* nested residuals, not Krylov reports."""
    values = {"COUPLED_NESTED_PCD_FP": [], "COUPLED_NESTED_PCD_KP": []}
    anchored = []
    for line in log.splitlines():
        prefix = line.split(" ", 1)[0]
        if prefix not in values:
            continue
        fields = {}
        for piece in line.split()[1:]:
            if "=" in piece:
                k, v = piece.split("=", 1)
                fields[k] = v
        key = (
            "original_true_relative_rhs"
            if prefix.endswith("_FP")
            else "independent_true_relative_rhs"
        )
        try:
            val = float(fields[key])
            accepted = int(fields["accepted"])
        except (KeyError, ValueError):
            raise ValueError("incomplete PCD inner true-residual trace")
        if not (0 <= val < float("inf") and accepted in (0, 1)):
            raise ValueError("nonfinite/invalid PCD inner true-residual trace")
        if prefix == "COUPLED_NESTED_PCD_FP" and "auxiliary_anchored" in fields:
            try:
                enabled = int(fields["auxiliary_anchored"])
                auxiliary = float(fields["auxiliary_true_relative_rhs"])
                anchors = int(fields["auxiliary_anchors"])
            except (KeyError, ValueError) as e:
                raise ValueError("incomplete experimental Fp auxiliary trace") from e
            if (
                enabled not in (0, 1)
                or not (0 <= auxiliary < float("inf"))
                or anchors < 0
            ):
                raise ValueError("nonfinite/invalid experimental Fp auxiliary trace")
            if enabled and (
                anchors == 0 or (accepted and (auxiliary > 1e-8 or val >= 0.999))
            ):
                raise ValueError(
                    "accepted anchored Fp action violated its independent residual gates"
                )
            anchored.append((enabled, auxiliary, anchors))
        values[prefix].append((val, accepted))
    summary = {
        mode: {
            "count": len(rows),
            "max_true_relative_rhs": max((v for v, _ in rows), default=None),
            "rejected_actions": sum(1 for _, accepted in rows if not accepted),
            "above_one_e_minus_three": sum(v > 1e-3 for v, _ in rows),
        }
        for mode, rows in values.items()
    }
    summary["fp_auxiliary_anchored"] = {
        "observations": len(anchored),
        "anchored_actions": sum(int(enabled) for enabled, _, _ in anchored),
        "max_auxiliary_true_relative_rhs": max(
            (v for enabled, v, _ in anchored if enabled), default=None
        ),
        "max_anchors": max((k for _, _, k in anchored), default=0),
        "original_residual_is_not_an_auxiliary_convergence_certificate": True,
    }
    return summary


def joint_returned_residual_summary(log: str) -> dict:
    """Summarize *returned* U,p,k,omega, never promote an inner SST/pressure pass.

    This is diagnostics, not a convergence or experimental certification.
    A truncated log retains incomplete status even when some local solves pass.
    """
    entries: list[dict] = []
    inner_sst_passes = 0
    inner_sst_failures = 0
    pressure_face_converged = 0
    for line in log.splitlines():
        prefix = line.split(" ", 1)[0]
        if prefix not in (
            "GLOBAL_NONLINEAR_V8",
            "BFS_PHASE_V10",
            "PRESSURE_HILO_PRODUCTION",
        ):
            continue
        fields = {}
        for piece in line.split()[1:]:
            if "=" in piece:
                k, v = piece.split("=", 1)
                fields[k] = v
        if prefix == "BFS_PHASE_V10" and fields.get("phase") == "SST_TRANSPORT_END":
            try:
                completed = int(fields["committed"]) and int(fields["converged"])
                if any(
                    int(fields[x]) not in (0, 1) for x in ("committed", "converged")
                ):
                    raise ValueError("invalid flags")
            except (KeyError, ValueError) as exc:
                raise ValueError("invalid SST transport phase evidence") from exc
            if completed:
                inner_sst_passes += 1
            else:
                inner_sst_failures += 1
        elif prefix == "PRESSURE_HILO_PRODUCTION":
            try:
                converged = int(fields["true_face_converged"])
                final = float(fields["final_face_relative"])
                tolerance = float(fields["requested_tolerance"])
            except (KeyError, ValueError) as exc:
                raise ValueError("incomplete true-face pressure evidence") from exc
            if (
                converged not in (0, 1)
                or not all(0 <= v < float("inf") for v in (final, tolerance))
                or (converged and final > tolerance)
            ):
                raise ValueError("invalid true-face pressure evidence")
            pressure_face_converged += int(converged == 1)
        elif prefix == "GLOBAL_NONLINEAR_V8":
            try:
                values = {
                    k: float(fields[k])
                    for k in ("U", "p_continuity", "k", "omega", "nut_closure", "merit")
                }
                flags = {k: int(fields[k]) for k in ("accepted", "qualified")}
                iteration = int(fields["iteration"])
            except (KeyError, ValueError) as exc:
                raise ValueError("incomplete independent U,p,k,omega evidence") from exc
            if iteration < 1 or any(not 0 <= v < float("inf") for v in values.values()):
                raise ValueError("invalid independent nonlinear residual")
            if any(v not in (0, 1) for v in flags.values()) or (
                flags["qualified"] and not flags["accepted"]
            ):
                raise ValueError("inconsistent independent nonlinear acceptance")
            # Merit is a maximum gate; allow tiny roundoff in display formatting.
            if values["merit"] + 1e-12 < max(
                values[k] for k in ("U", "p_continuity", "k", "omega", "nut_closure")
            ) * (1 - 1e-8):
                raise ValueError(
                    "nonlinear merit hides a larger physical field residual"
                )
            worst = max(
                ("U", "p_continuity", "k", "omega", "nut_closure"),
                key=lambda k: values[k],
            )
            entries.append(
                {
                    "iteration": iteration,
                    "fields": values,
                    "accepted": bool(flags["accepted"]),
                    "qualified": bool(flags["qualified"]),
                    "worst_field": worst,
                }
            )
    return {
        "independent_global_observations": len(entries),
        "first_returned_state": entries[0] if entries else None,
        "last_returned_state": entries[-1] if entries else None,
        "best_observed_merit": min(
            (e["fields"]["merit"] for e in entries), default=None
        ),
        "global_accepted_observations": sum(e["accepted"] for e in entries),
        "global_qualified_observations": sum(e["qualified"] for e in entries),
        "sst_inner_completed_observations": inner_sst_passes,
        "sst_inner_failed_observations": inner_sst_failures,
        "true_face_pressure_converged_observations": pressure_face_converged,
        "no_full_nonlinear_certificate_from_partial_local_solves": True,
    }


def run_one(
    executable: Path,
    mesh: Path,
    output: Path,
    iterations: int,
    timeout_seconds: float,
    *,
    inner_trace: bool,
) -> dict:
    output.mkdir(parents=True, exist_ok=True)
    logfile = output / "production.log"
    cmd = [
        str(executable),
        "--mesh",
        str(mesh),
        "--output-dir",
        str(output),
        "--iterations",
        str(iterations),
        "--experimental-sst",
        "--benchmark-export",
    ]
    env = os.environ.copy()
    # Campaigns must never inherit stale numerical experiments from an
    # unrelated NACA/BFS shell. The case HDF5 is the authoritative setup.
    for name in tuple(env):
        if name.startswith("CFDX_"):
            del env[name]
    env["CFDX_BFS_STAGE_TRACE"] = "1"
    if inner_trace:
        env["CFDX_BFS_PCD_INNER_TRACE"] = "1"
    # This trace exposes the distinct update/returned-equation/linear gates.
    env["CFDX_SST_DIAGNOSTICS"] = "1"
    start = time.monotonic()
    status = "NOT_RUN"
    code = None
    with logfile.open("w", encoding="utf-8") as stream:
        proc = subprocess.Popen(
            cmd,
            env=env,
            stdout=stream,
            stderr=subprocess.STDOUT,
            start_new_session=True,
        )
        try:
            code = proc.wait(timeout=timeout_seconds)
            status = "EXIT_ZERO" if code == 0 else "FAILED"
        except subprocess.TimeoutExpired:
            os.killpg(proc.pid, signal.SIGKILL)
            proc.wait()
            status = "TIMEOUT_INCOMPLETE"
    body = logfile.read_text(encoding="utf-8", errors="replace")
    try:
        pcd_metrics = nested_pcd_summary(body)
        pcd_parse_error = None
    except ValueError as exc:
        # A deadline can stop the process mid-output. Preserve the failed
        # experiment and its raw log, rather than losing the entire matrix.
        pcd_metrics = None
        pcd_parse_error = str(exc)
    try:
        joint_summary = joint_returned_residual_summary(body)
        joint_parse_error = None
    except ValueError as exc:
        joint_summary = None
        joint_parse_error = str(exc)
    return {
        "execution": status,
        "exit_code": code,
        "wall_seconds": round(time.monotonic() - start, 3),
        "command": cmd,
        "log": str(logfile),
        "log_sha256": sha(logfile),
        "phase_observation": classify_log(body, execution=status, exit_code=code),
        "nested_pcd_true_residuals": pcd_metrics,
        "returned_equations": joint_summary,
        "returned_equations_parse_error": joint_parse_error,
        "nested_pcd_trace_error": pcd_parse_error,
        "pcd_fp_nested_records": body.count("COUPLED_NESTED_PCD_FP"),
        "pcd_kp_nested_records": body.count("COUPLED_NESTED_PCD_KP"),
        "pressure_true_face_records": body.count("PRESSURE_HILO_PRODUCTION"),
        "sst_acceptance_gate_records": body.count("SST_ACCEPTANCE_GATES"),
        "status_note": (
            "No physical qualification from a single step or timeout. "
            "Use independent V8 certificate for completed runs."
        ),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-case", type=Path, required=True)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--profiles", nargs="+", choices=tuple(PROFILES), default=list(PROFILES)
    )
    parser.add_argument("--outer-iterations", type=int, default=1)
    parser.add_argument("--sst-inner-budget", type=int, default=800)
    parser.add_argument("--wall-timeout", type=float, default=60.0)
    parser.add_argument("--prepare-only", action="store_true")
    args = parser.parse_args(argv)
    if args.outer_iterations < 1 or args.sst_inner_budget < 1 or args.wall_timeout <= 0:
        parser.error("positive numerical and wall-time budgets required")
    if len(args.profiles) != len(set(args.profiles)):
        parser.error("duplicate profiles would overwrite their evidence")
    if not args.reference_case.is_file():
        parser.error("reference HDF5 case missing")
    if not args.prepare_only and not args.binary.is_file():
        parser.error("production binary missing")
    args.out.mkdir(parents=True, exist_ok=True)
    physical_reference = case_physics_signature(args.reference_case)
    results = []
    for profile in args.profiles:
        algo, schur, pc, _ = PROFILES[profile]
        mesh = args.out / "cases" / (profile + ".cfdx.h5")
        prepared = prepare_case(
            args.reference_case,
            mesh,
            profile,
            args.outer_iterations,
            sst_inner_budget=args.sst_inner_budget,
        )
        result = {
            "profile": profile,
            "prepared": prepared,
            "physical_reference_sha256": physical_reference,
            "physical_qualification": "NOT_QUALIFIED",
        }
        if args.prepare_only:
            result["execution"] = "PREPARED_NOT_RUN"
        else:
            output = args.out / "runs" / profile
            result.update(
                run_one(
                    args.binary,
                    mesh,
                    output,
                    args.outer_iterations,
                    args.wall_timeout,
                    inner_trace=profile.startswith("pcd_"),
                )
            )
            result["v8_nonlinear_gate"] = validate_bfs_evidence(
                output,
                result["exit_code"],
                algorithm=algo,
                schur=schur,
                preconditioner=pc,
                case_sha256=prepared["case_sha256"],
                case_path=mesh,
            )
            result["physical_qualification"] = (
                "NOT_QUALIFIED_EXPERIMENTAL_VALIDATION_PENDING"
            )
        results.append(result)
        (args.out / "matrix_v11.json").write_text(
            json.dumps(
                {
                    "schema": "CFDX-BFS-V11-NUMERICAL-AB/v1",
                    "reference_case_sha256": sha(args.reference_case),
                    "reference_physical_signature_sha256": physical_reference,
                    "binary_sha256": sha(args.binary)
                    if args.binary.is_file()
                    else None,
                    "all_tolerances_unchanged": True,
                    "physical_qualification": "NOT_QUALIFIED",
                    "results": results,
                },
                indent=2,
                allow_nan=False,
            )
            + "\n"
        )
        print(f"BFS_V11 {profile}: {result.get('execution', 'NOT_RUN')}", flush=True)
    # Even a completed algebraic experiment is not a NASA comparison.
    return (
        0
        if args.prepare_only
        else (
            0
            if all(
                v.get("v8_nonlinear_gate", {}).get("nonlinear_v8_proven") is True
                for v in results
            )
            else 2
        )
    )


if __name__ == "__main__":
    sys.exit(main())
