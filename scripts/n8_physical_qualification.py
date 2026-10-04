#!/usr/bin/env python3
"""Run the reproducible N8 physical/production qualification campaign.

The campaign deliberately reuses existing CTest validation executables. It does
not duplicate solver implementations or alter their tolerances. The script
collects their stdout/stderr and emits a machine-readable JSON report.

Usage:
  python3 scripts/n8_physical_qualification.py --build-dir build
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path


REQUIRED_TESTS = (
    # Physical production paths: use the full pressure-velocity matrix, not
    # only the short PR smoke test, so every exposed coupling algorithm is
    # exercised (SIMPLE, SIMPLEC, PISO, PIMPLE, fractional-step, COUPLED).
    "test_n8_pressure_velocity_matrix",
    "test_poiseuille_quick",
    "test_ghia_cavity_quick",
    "test_nonorthogonal_skew_campaign",
    # Linear solver families used by N8.
    "test_cg_solver",
    "test_bicgstab_solver",
    "test_gmres_solver",
    "test_krylov_preconditioning",
    # N8 preconditioner / Schur families.
    "test_amg_preconditioner_qualification",
    "test_advanced_preconditioners",
    "test_schur_preconditioner",
    "test_schur_infrastructure",
    "test_exact_schur",
    "test_simplerc_schur",
    "test_lsc_bfbt_schur",
    "test_lsc_bfbt_schur_null_space",
    "test_schur_approximation_comparison",
    "test_schur_quantitative_qualification",
    "test_mgr_preconditioner",
    "test_coupled_block_schur_amg",
    "test_n8_schur_production_benchmark",
)


def ctest(build_dir: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        ["ctest", *args],
        cwd=build_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )


def discover_tests(build_dir: Path) -> set[str]:
    result = ctest(build_dir, "-N")
    if result.returncode != 0:
        raise RuntimeError(
            f"ctest -N failed with exit code {result.returncode}\n{result.stdout}"
        )

    names: set[str] = set()
    for line in result.stdout.splitlines():
        line = line.strip()
        if not line.startswith("Test #"):
            continue
        try:
            # ctest -N emits e.g. "Test #81: test_gmres_solver".
            _, indexed_name = line.split("#", 1)
            _, name = indexed_name.split(":", 1)
        except ValueError:
            continue
        names.add(name.strip())
    return names


def run_test(build_dir: Path, name: str) -> dict[str, object]:
    start = time.monotonic()
    result = ctest(build_dir, "--output-on-failure", "--tests-regex", f"^{name}$")
    elapsed = time.monotonic() - start
    return {
        "name": name,
        "returncode": result.returncode,
        "status": "PASS" if result.returncode == 0 else "FAIL",
        "elapsed_s": round(elapsed, 3),
        "output": result.stdout,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="JSON report path (default: <build-dir>/n8_physical_qualification.json)",
    )
    args = parser.parse_args()

    build_dir = args.build_dir.resolve()
    if not (build_dir / "CTestTestfile.cmake").exists():
        print(f"error: {build_dir} is not a configured CMake build directory", file=sys.stderr)
        return 2

    report_path = (
        args.report.resolve()
        if args.report is not None
        else build_dir / "n8_physical_qualification.json"
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)

    available = discover_tests(build_dir)
    missing = [name for name in REQUIRED_TESTS if name not in available]
    if missing:
        report = {
            "campaign": "N8 complete solver/preconditioner qualification",
            "status": "INCOMPLETE",
            "required_tests": list(REQUIRED_TESTS),
            "missing_tests": missing,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 2

    results = []
    for name in REQUIRED_TESTS:
        print(f"\n=== N8 qualification: {name} ===", flush=True)
        result = run_test(build_dir, name)
        results.append(result)
        print(result["output"], end="", flush=True)
        if result["returncode"] != 0:
            # Stop at the first failed gate. The complete output is retained in
            # the JSON report for diagnosis; later tests are not silently run.
            break

    failed = [r["name"] for r in results if r["status"] == "FAIL"]
    completed = len(results)
    report = {
        "campaign": "N8 complete solver/preconditioner qualification",
        "status": "PASS" if not failed and completed == len(REQUIRED_TESTS) else "FAIL",
        "required_tests": list(REQUIRED_TESTS),
        "completed_tests": completed,
        "failed_tests": failed,
        "results": results,
        "coverage": {
            "pressure_velocity": [
                "SIMPLE", "SIMPLEC", "PISO", "PIMPLE",
                "FRACTIONAL_STEP", "COUPLED"
            ],
            "linear_solvers": ["CG", "BiCGStab", "GMRES", "FGMRES"],
            "preconditioners": [
                "Native AMG", "Smoothed Aggregation AMG",
                "Native FieldSplit", "Coupled Block Schur", "MGR"
            ],
            "schur": [
                "Exact", "SIMPLE/SIMPLEC", "LSC", "BFBt",
                "null-space", "quantitative comparison"
            ],
            "physical_cases": ["Couette", "Poiseuille", "Ghia Re=100", "skew/non-orthogonal"]
        },
        "policy": {
            "reuses_existing_ctest_tests": True,
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_fallbacks": False,
            "independent_true_residuals": "provided by existing tests where applicable",
            "solver_vs_discretisation_separation": True,
        },
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nN8 qualification report: {report_path}")
    print(f"N8 qualification status: {report['status']}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
