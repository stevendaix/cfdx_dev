"""Run bounded, traceable COUPLED comparisons. NEVER promote a short smoke test.

The N9 executable itself runs Couette, Poiseuille, Ghia and skew Couette in
sequence for each Schur model; a timeout is a failed/incomplete whole shard.
The production NACA driver is the actual finite-volume solver, not the small
Newton-Krylov algebraic unit test. No silent fallback, relaxed tolerance or
fabricated physical qualification is permitted.
"""

from __future__ import annotations

import argparse
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
        "CFDX_NACA_COUPLED_SCHUR": "pcd", "CFDX_NACA_PCD_FP_EQUILIBRATE": "1"},
    "pcd_ilu0_guarded": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd", "CFDX_NACA_PCD_FP_QUALITY_GUARD": "1"},
    "pcd_ilu0_refined": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd", "CFDX_NACA_PCD_FP_SECOND_PASS": "1"},
    "pcd_amg": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd", "CFDX_NACA_PCD_FP_PRECONDITIONER": "amg"},
    "pcd_lumped": {
        "CFDX_NACA_COUPLED_SCHUR": "pcd", "CFDX_NACA_PCD_LUMPED_MASS": "1"},
    "block_local": {"CFDX_NACA_COUPLED_SCHUR": "block_local"},
    "lsc": {"CFDX_NACA_COUPLED_SCHUR": "lsc"},
    "bfbt": {"CFDX_NACA_COUPLED_SCHUR": "bfbt"},
    "simple_schur": {"CFDX_NACA_COUPLED_SCHUR": "simple"},
    "simplec_schur": {"CFDX_NACA_COUPLED_SCHUR": "simplec"},
}
N9_MODELS = ("pcd", "lsc", "bfbt", "simple", "simplec")
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
    if returncode == 0 and re.search(
        r"NACA_SOLVE_WALL[^\n]*\bconverged=1\b", content
    ):
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


def run_case(
    *, label: str, exe: Path | None, args: list[str], env: dict[str, str],
    output: Path, timeout: int, naca: bool
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
    try:
        process = subprocess.run(
            [str(exe), *args], env=combined_env,
            timeout=timeout, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            check=False, text=True, errors="replace",
        )
        content, code = process.stdout, process.returncode
    except subprocess.TimeoutExpired as exc:
        payload = exc.stdout or b""
        content = payload.decode("utf-8", "replace") if isinstance(payload, bytes) else payload
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
    else:
        row["status"] = (
            "N9_PHYSICAL_SHARD_TEST_PASSED" if code == 0 else
            "TIMEOUT_N9_SHARD_INCOMPLETE" if code == 124 else
            "N9_PHYSICAL_SHARD_FAILED"
        )
        row["cases_in_shard"] = ["Couette", "Poiseuille", "Ghia_Re100", "skew_Couette"]
        row["last_stage_detected"] = (
            next((name for name in reversed(row["cases_in_shard"])
                  if name.replace("Ghia_Re100", "Ghia Re=100").replace("skew_Couette", "skew Couette")
                  in content), None)
        )
    return row


def build_report(args: argparse.Namespace) -> dict:
    naca_exe, mesh, n9_exe = args.naca_bin, args.naca_mesh, args.n9_bin
    report = {
        "schema": "cfdx-coupled-comparative-matrix/v1",
        "source_revision": args.source_revision,
        "mesh_sha256": sha256(mesh) if mesh and mesh.is_file() else None,
        "naca_max_nonlinear_iterations": args.naca_iterations,
        "timeout_seconds_per_case": args.timeout,
        "not_physical_qualification": True,
        "cases": [],
    }
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
        if mesh and mesh.is_file():
            report["cases"].append(run_case(
                label="NACA_N64_" + name, exe=naca_exe,
                args=["--single", str(mesh), "0.6", "0.25", str(args.naca_iterations)],
                env=flags, output=args.output_dir / f"naca_{name}.log",
                timeout=args.timeout, naca=True,
            ))
        else:
            report["cases"].append(run_case(
                label="NACA_N64_" + name, exe=None, args=[], env=flags,
                output=args.output_dir / f"naca_{name}.log",
                timeout=args.timeout, naca=True,
            ))
    for name in args.n9_models:
        report["cases"].append(run_case(
            label="N9_Couette_Poiseuille_Ghia_skew_" + name,
            exe=n9_exe, args=["--schur-model", name], env={},
            output=args.output_dir / f"n9_{name}.log",
            timeout=args.timeout, naca=False,
        ))
    report["summary"] = {
        "requested": len(report["cases"]),
        "naca_converged_one_grid": sum(
            r["status"] == "SOLVER_CONVERGED_ONE_GRID_ONLY" for r in report["cases"]),
        "n9_shards_pass": sum(r["status"] == "N9_PHYSICAL_SHARD_TEST_PASSED" for r in report["cases"]),
        "incomplete_or_failed": sum(
            r["status"] not in ("SOLVER_CONVERGED_ONE_GRID_ONLY", "N9_PHYSICAL_SHARD_TEST_PASSED")
            for r in report["cases"]),
        "full_3grid_NACA_GCI": "NOT_RUN",
        "independent_physical_reference_validation": "NOT_RUN",
    }
    return report


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--naca-bin", type=Path)
    p.add_argument("--naca-mesh", type=Path)
    p.add_argument("--n9-bin", type=Path)
    p.add_argument("--source-revision", default="UNATTESTED_ZIP")
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--timeout", type=int, default=30)
    p.add_argument("--naca-iterations", type=int, default=2)
    p.add_argument("--naca-variants", nargs="*", choices=tuple(NACA_VARIANTS),
                   default=["pcd_ilu0", "pcd_ilu0_equilibrated", "pcd_amg", "block_local"])
    p.add_argument("--n9-models", nargs="*", choices=N9_MODELS,
                   default=["pcd", "lsc", "bfbt", "simple", "simplec"])
    args = p.parse_args()
    if args.timeout < 1 or args.naca_iterations < 1:
        p.error("timeout and iteration budget must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    result = build_report(args)
    target = args.output_dir / "coupled_matrix.json"
    target.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(f"CFDX_COUPLED_MATRIX: {result['summary']} report={target}")
    # Keep diagnostic runs scriptable even when most cases fail or time out;
    # summary never labels partial work as physically qualified.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
