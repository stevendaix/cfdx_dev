"""Run bounded, traceable COUPLED comparisons. NEVER promote a short smoke test.

The N9 executable itself runs Couette, Poiseuille, Ghia and skew Couette in
sequence for each Schur model; a timeout is a failed/incomplete whole shard.
The production NACA driver is the actual finite-volume solver, not the small
Newton-Krylov algebraic unit test. No silent fallback, relaxed tolerance or
fabricated physical qualification is permitted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import re
import subprocess
import time
from pathlib import Path

NACA_VARIANTS = {
    "pcd_ilu0": {"CFDX_NACA_COUPLED_SCHUR": "pcd"},
    "pcd_ilu0_equilibrated": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd",
        "CFDX_NACA_PCD_FP_EQUILIBRATE": "1",
    },
    "pcd_ilu0_guarded": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd",
        "CFDX_NACA_PCD_FP_QUALITY_GUARD": "1",
    },
    "pcd_ilu0_refined": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd",
        "CFDX_NACA_PCD_FP_SECOND_PASS": "1",
    },
    "pcd_amg": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd",
        "CFDX_NACA_PCD_FP_PRECONDITIONER": "amg",
    },
    "pcd_lumped": {"CFDX_NACA_COUPLED_SCHUR": "pcd", "CFDX_NACA_PCD_LUMPED_MASS": "1"},
    "block_local": {"CFDX_NACA_COUPLED_SCHUR": "block_local"},
    "lsc": {"CFDX_NACA_COUPLED_SCHUR": "lsc"},
    "bfbt": {"CFDX_NACA_COUPLED_SCHUR": "bfbt"},
    "simple_schur": {"CFDX_NACA_COUPLED_SCHUR": "simple"},
    "simplec_schur": {"CFDX_NACA_COUPLED_SCHUR": "simplec"},
}
N9_MODELS = ("pcd", "lsc", "bfbt", "simple", "simplec", "block_local")
N9_STAGES = ("Couette", "Poiseuille", "Ghia_Re100", "skew_Couette")
RE_TRUE_FP = re.compile(r"\boriginal_true_relative_rhs=([0-9.eE+-]+)")
RE_ACCEPT = re.compile(r"\baccepted=(0|1)\b")
RE_ARNOLDI = re.compile(r"KRYLOV_(?:UNHAPPY_ARNOLDI_BREAKDOWN|ARNOLDI_BREAKDOWN)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def classify_naca(returncode: int | None, content: str) -> str:
    if returncode is None:
        return "NOT_RUN_MISSING_BINARY_OR_MESH"
    if returncode == 124:
        return "TIMEOUT_NOT_CONVERGED"
    if returncode == 0 and re.search(r"NACA_SOLVE_WALL[^\n]*\bconverged=1\b", content):
        return "SOLVER_CONVERGED_ONE_GRID_ONLY"
    if returncode == 0:
        return "UNVERIFIED_SUCCESS_NO_PHYSICAL_MARKER"
    return "FAILED_OR_SHORT_DIAGNOSTIC_NOT_CONVERGED"


def summarize_fp(content: str) -> dict:
    values = [float(x) for x in RE_TRUE_FP.findall(content)]
    values = [x for x in values if math.isfinite(x)]
    accepts = RE_ACCEPT.findall(content)
    return {
        "observed_applications": len(values),
        "true_relative_rhs_min": min(values) if values else None,
        "true_relative_rhs_max": max(values) if values else None,
        "true_relative_rhs_above_1e_3": sum(x > 1e-3 for x in values),
        "quality_rejected": sum(x == "0" for x in accepts),
        "arnoldi_breakdown_markers": len(RE_ARNOLDI.findall(content)),
    }


def summarize_naca_history(path: Path) -> dict:
    """Summarize the ACTUAL accepted FV states, never the linear log alone."""
    if not path.is_file():
        return {"status": "MISSING", "last_accepted_iteration": None}
    required = (
        "iteration",
        "algorithm",
        "convection_scheme",
        "high_order_fraction",
        "momentum_rhs_relative",
        "momentum_component_max",
        "face_momentum_x_relative",
        "face_momentum_y_relative",
        "continuity_normalized",
        "mass_local_linf",
        "pressure_residual",
        "nonfinite_mass_faces",
        "nonfinite_velocities",
    )
    try:
        with path.open(newline="", encoding="utf-8") as stream:
            reader = csv.DictReader(stream)
            if not reader.fieldnames or not set(required).issubset(reader.fieldnames):
                raise ValueError("missing independent physical history columns")
            count = 0
            final = None
            for row in reader:
                count += 1
                if int(row["iteration"]) != count:
                    raise ValueError("nonconsecutive accepted FV states")
                final = row
            if not count or final is None:
                raise ValueError("empty accepted FV history")
        metrics = {}
        for field in required:
            if field in ("iteration", "algorithm", "convection_scheme"):
                continue
            value = float(final[field])
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(f"nonfinite or invalid {field}")
            metrics[field] = value
        if metrics["high_order_fraction"] > 1.0:
            raise ValueError("high-order fraction above unity")
        return {
            "status": "OBSERVED_NOT_INDEPENDENTLY_CERTIFIED",
            "last_accepted_iteration": count,
            "algorithm": final["algorithm"],
            "final_convection_scheme": final["convection_scheme"],
            "metrics": metrics,
        }
    except (OSError, ValueError, KeyError, TypeError) as exc:
        return {
            "status": "INVALID",
            "reason": str(exc),
            "last_accepted_iteration": None,
        }


def final_naca_muscl_contract_valid(history: dict) -> bool:
    """Fail closed on the returned FULL-MUSCL FV state, never a linear status."""
    if history.get("status") != "OBSERVED_NOT_INDEPENDENTLY_CERTIFIED":
        return False
    if (
        history.get("algorithm") != "COUPLED"
        or history.get("final_convection_scheme") != "MUSCL_CELL_LIMITED"
    ):
        return False
    metrics = history.get("metrics", {})
    if metrics.get("high_order_fraction") != 1.0:
        return False
    required_gates = (
        "momentum_rhs_relative",
        "momentum_component_max",
        "face_momentum_x_relative",
        "face_momentum_y_relative",
        "continuity_normalized",
        "mass_local_linf",
    )
    return (
        all(
            isinstance(metrics.get(name), (float, int))
            and math.isfinite(metrics[name])
            and 0.0 <= metrics[name] <= 1e-8
            for name in required_gates
        )
        and metrics.get("nonfinite_mass_faces") == 0.0
        and metrics.get("nonfinite_velocities") == 0.0
    )


def summarize_n9_physical_shard(
    content: str, selected_model: str, exit_code: int
) -> dict:
    """Only post-oracle N9 COMPLETED records count; STARTING is never PASS.

    This rejects a stale binary, a different Schur selection, duplicate rows,
    missing later cases and nonfinite/residual-exceeding measurements.
    """
    prefix = f"COMPLETED COUPLED {selected_model} "
    fields = (
        "momentum_rhs_relative",
        "continuity_normalized",
        "mass_local_linf",
        "face_momentum_x_relative",
        "face_momentum_y_relative",
        "iterations",
    )
    records: dict[str, dict] = {}
    errors: list[str] = []
    requested_marker = f"STARTING COUPLED {selected_model}"
    if requested_marker not in content.splitlines():
        errors.append("missing exact requested-model STARTING identity")
    for line in content.splitlines():
        if not line.startswith("COMPLETED COUPLED "):
            continue
        if not line.startswith(prefix):
            errors.append("different Schur model in completion")
            continue
        rest = line[len(prefix) :].split()
        if not rest:
            errors.append("completion without case")
            continue
        case = rest.pop(0)
        if case not in N9_STAGES or case in records:
            errors.append(f"unknown or duplicate case {case}")
            continue
        try:
            parsed = [item.split("=", 1) for item in rest]
            if any(len(pair) != 2 for pair in parsed):
                raise ValueError("malformed physical metric")
            names = [pair[0] for pair in parsed]
            if len(names) != len(set(names)):
                raise ValueError("duplicated physical metric field")
            raw = dict(parsed)
            if set(raw) != set(fields):
                raise ValueError("missing or unexpected physical metric fields")
            data = {k: float(v) for k, v in raw.items()}
            if any(not math.isfinite(v) or v < 0 for v in data.values()):
                raise ValueError("nonfinite or negative physical metric")
            if data["iterations"] < 1 or not data["iterations"].is_integer():
                raise ValueError("invalid iteration count")
            # N9's own strict physical gates use 1e-7 on continuity and
            # momentum, but we do not hardcode a weaker substitute for that
            # code-side policy: this is a fail-closed extra report check.
            records[case] = data
            if (
                data["momentum_rhs_relative"] > 1e-8
                or data["continuity_normalized"] > 1e-8
            ):
                errors.append(
                    f"{case}: returned FV momentum/continuity exceeds 1e-8 qualification gate"
                )
            if (
                data["face_momentum_x_relative"] > 1e-8
                or data["face_momentum_y_relative"] > 1e-8
            ):
                errors.append(
                    f"{case}: independent face momentum exceeds 1e-8 qualification gate"
                )
            if data["mass_local_linf"] > 1e-8:
                errors.append(
                    f"{case}: local mass imbalance exceeds 1e-8 qualification gate"
                )
        except (ValueError, TypeError) as exc:
            errors.append(f"{case}: {exc}")
    completed = [n for n in N9_STAGES if n in records]
    complete = exit_code == 0 and len(records) == len(N9_STAGES) and not errors
    return {
        "status": "FOUR_PHYSICAL_CASES_VERIFIED"
        if complete
        else "INCOMPLETE_OR_UNVERIFIED",
        "selected_schur_model": selected_model,
        "completed_stages": completed,
        "case_metrics": records,
        "issues": errors,
        "physical_shard_gate_complete": complete,
    }


def run_case(
    *,
    label: str,
    exe: Path | None,
    args: list[str],
    env: dict[str, str],
    output: Path,
    timeout: int,
    naca: bool,
) -> dict:
    row = {
        "name": label,
        "status": "NOT_RUN_MISSING_BINARY_OR_MESH",
        "source": "CURRENT_WORKTREE_REBUILD_REQUIRED",
        "binary": str(exe) if exe else None,
        "binary_sha256": sha256(exe) if exe and exe.is_file() else None,
        "argv": [str(exe), *args] if exe else None,
        "environment_overrides": dict(sorted(env.items())),
        "elapsed_s": None,
        "exit_code": None,
        "physical_qualification": "NOT_QUALIFIED_BY_THIS_RUNNER",
    }
    if exe is None or not exe.is_file():
        return row
    output.parent.mkdir(parents=True, exist_ok=True)
    started = time.monotonic()
    combined_env = os.environ.copy()
    # Always clear potentially inherited CFDX_NACA experiment flags so each
    # variant runs exactly the requested configuration, without hidden options.
    for key in list(combined_env):
        if key.startswith("CFDX_NACA_") or key == "CFDX_KRYLOV_DIAGNOSTICS":
            combined_env.pop(key)
    combined_env.update(env)
    if naca:
        history = output.with_suffix(".history.csv")
        force = output.with_suffix(".forces.csv")
        raw_faces = output.with_suffix(".raw_wall_faces.csv")
        # Stale evidence from a previous run must never be promoted.
        for path in (history, force, raw_faces):
            path.unlink(missing_ok=True)
        combined_env.update(
            {
                "CFDX_NACA_PHYSICAL_HISTORY_CSV": str(history),
                "CFDX_NACA_CONVERGED_FORCE_CSV": str(force),
                "CFDX_NACA_RAW_WALL_FACES_CSV": str(raw_faces),
                "CFDX_NACA_WALL_TRACTION_FACE": "1",
                "CFDX_NACA_WALL_GRADIENT_WLS": "1",
            }
        )
    try:
        process = subprocess.run(
            [str(exe), *args],
            env=combined_env,
            timeout=timeout if timeout else None,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
            text=True,
            errors="replace",
        )
        content, code = process.stdout, process.returncode
    except subprocess.TimeoutExpired as exc:
        payload = exc.stdout or b""
        content = (
            payload.decode("utf-8", "replace")
            if isinstance(payload, bytes)
            else payload
        )
        content += "\nDIAGNOSTIC_TIME_BUDGET_EXHAUSTED\n"
        code = 124
    row["elapsed_s"] = round(time.monotonic() - started, 3)
    row["exit_code"] = code
    row["log_file"] = str(output)
    output.write_text(content, encoding="utf-8")
    if naca:
        row["status"] = classify_naca(code, content)
        row["fp"] = summarize_fp(content)
        row["numeric_scheme"] = "MUSCL_CELL_LIMITED"
        row["pressure_farfield"] = "NEUMANN"
        row["linear_tolerance"] = 1e-12
        row["physical_history"] = summarize_naca_history(history)
        row["history_sha256"] = sha256(history) if history.is_file() else None
        row["force_csv_sha256"] = sha256(force) if force.is_file() else None
        row["raw_wall_faces_sha256"] = (
            sha256(raw_faces) if raw_faces.is_file() else None
        )
        row["independent_force_oracle_status"] = "NOT_RUN"
        # No emission of force evidence is a valid outcome of a short or
        # failed run. Never substitute previous files or a CTest PASS.
        row["independent_physical_evidence_complete"] = False
        # A log line alone is insufficient: never mark the one-grid solver as
        # converged if the *new* final FV history or surface data are missing.
        if row["status"] == "SOLVER_CONVERGED_ONE_GRID_ONLY" and (
            not final_naca_muscl_contract_valid(row["physical_history"])
            or not force.is_file()
            or not raw_faces.is_file()
        ):
            row["status"] = "INCOMPLETE_FV_EVIDENCE_NOT_QUALIFIED"
    else:
        requested_model = (
            args[1] if len(args) >= 2 and args[0] == "--schur-model" else "UNSELECTED"
        )
        phys = summarize_n9_physical_shard(content, requested_model, code)
        row["n9_physical_cases"] = phys
        row["status"] = (
            "N9_PHYSICAL_SHARD_TEST_PASSED"
            if phys["physical_shard_gate_complete"]
            else "TIMEOUT_N9_SHARD_INCOMPLETE"
            if code == 124
            else "N9_PHYSICAL_SHARD_FAILED"
            if code != 0
            else "N9_UNVERIFIED_ZERO_EXIT"
        )
        row["cases_in_shard"] = list(N9_STAGES)
        row["last_stage_detected"] = (
            phys["completed_stages"][-1] if phys["completed_stages"] else None
        )
    return row


def build_report(args: argparse.Namespace) -> dict:
    naca_exe, mesh, n9_exe = args.naca_bin, args.naca_mesh, args.n9_bin
    report = {
        "schema": "cfdx-coupled-comparative-matrix/v2",
        "source_revision": args.source_revision,
        "mesh_sha256": sha256(mesh) if mesh and mesh.is_file() else None,
        "naca_max_nonlinear_iterations": args.naca_iterations,
        "timeout_seconds_per_case": args.timeout or None,
        "shared_FV_contract": {
            "scheme": "MUSCL_CELL_LIMITED",
            "pressure": "NEUMANN",
            "linear_tolerance": 1e-12,
            "alpha_u": 0.6,
            "alpha_p": 0.25,
            "wall_forces": "FACE_WLS_ONLY_ON_CONVERGED_SOLUTIONS",
        },
        "not_physical_qualification": True,
        "cases": [],
    }
    mesh_set = [(64, mesh), (128, args.naca_n128_mesh), (256, args.naca_n256_mesh)]
    for name in args.naca_variants:
        flags = {
            "CFDX_NACA_ALGORITHM": "COUPLED",
            "CFDX_NACA_CONVECTION_CELL_LIMITED": "1",
            "CFDX_NACA_LINEAR_TOLERANCE": "1e-12",
            "CFDX_NACA_FARFIELD_PRESSURE_NEUMANN": "1",
            "CFDX_NACA_REQUIRE_KINEMATIC_FACE_FLUX": "1",
            "CFDX_NACA_COUPLED_MATRIX_SUMMARY": "1",
            **NACA_VARIANTS[name],
        }
        for level, grid in mesh_set:
            if level != 64 and grid is None:
                continue  # only run additional levels when explicitly requested
            label = f"NACA_N{level}_{name}"
            report["cases"].append(
                run_case(
                    label=label,
                    exe=naca_exe if grid and grid.is_file() else None,
                    args=[
                        "--single",
                        str(grid),
                        "0.6",
                        "0.25",
                        str(args.naca_iterations),
                    ]
                    if grid
                    else [],
                    env=flags,
                    output=args.output_dir / f"naca_n{level}_{name}.log",
                    timeout=args.timeout,
                    naca=True,
                )
            )
    for name in args.n9_models:
        report["cases"].append(
            run_case(
                label="N9_Couette_Poiseuille_Ghia_skew_" + name,
                exe=n9_exe,
                args=["--schur-model", name, "--require-independent-face"],
                env={},
                output=args.output_dir / f"n9_{name}.log",
                timeout=args.timeout,
                naca=False,
            )
        )
    report["summary"] = {
        "requested": len(report["cases"]),
        "naca_converged_one_grid": sum(
            r["status"] == "SOLVER_CONVERGED_ONE_GRID_ONLY" for r in report["cases"]
        ),
        "n9_shards_pass": sum(
            r["status"] == "N9_PHYSICAL_SHARD_TEST_PASSED" for r in report["cases"]
        ),
        "incomplete_or_failed": sum(
            r["status"]
            not in ("SOLVER_CONVERGED_ONE_GRID_ONLY", "N9_PHYSICAL_SHARD_TEST_PASSED")
            for r in report["cases"]
        ),
        "full_3grid_NACA_GCI": "NOT_RUN",
        "independent_physical_reference_validation": "NOT_RUN",
    }
    return report


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--naca-bin", type=Path)
    p.add_argument("--naca-mesh", type=Path, help="N64 mesh")
    p.add_argument("--naca-n128-mesh", type=Path)
    p.add_argument("--naca-n256-mesh", type=Path)
    p.add_argument("--n9-bin", type=Path)
    p.add_argument("--source-revision", default="UNATTESTED_ZIP")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument(
        "--timeout",
        type=int,
        default=0,
        help="wall limit per method; 0=unlimited; never changes iteration gates",
    )
    p.add_argument("--naca-iterations", type=int, default=1500)
    p.add_argument(
        "--naca-variants",
        nargs="*",
        choices=tuple(NACA_VARIANTS),
        default=["pcd_ilu0", "pcd_ilu0_equilibrated", "pcd_amg", "block_local"],
    )
    p.add_argument(
        "--n9-models",
        nargs="*",
        choices=N9_MODELS,
        default=["pcd", "lsc", "bfbt", "simple", "simplec", "block_local"],
    )
    args = p.parse_args()
    if args.timeout < 0 or args.naca_iterations < 1:
        p.error("timeout must be nonnegative and iteration budget must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = build_report(args)
    target = args.output_dir / "coupled_matrix.json"
    target.write_text(
        json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    )
    print(f"CFDX_COUPLED_MATRIX: {result['summary']} report={target}")
    # Keep diagnostic runs scriptable even when most cases fail or time out;
    # summary never labels partial work as physically qualified.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
