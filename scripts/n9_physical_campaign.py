#!/usr/bin/env python3
"""Run the reproducible N9 pressure-velocity physical qualification campaign.

The campaign reuses existing executable CFDX tests. It does not duplicate
solver or physics implementations and never changes tolerances. Each test is
run verbosely and its complete output is retained in a machine-readable JSON
report, together with the exact source SHA the evidence was produced on.

Gate policy (fail-closed):
- PASS requires every required test to be discovered, executed and passed,
  and the campaign source SHA to match --expected-git-sha.
- The optional MPI equivalence test is required whenever the build exposes
  it (no silent MPI skip).
- Case/algorithm coverage markers found in the physical matrix output are
  recorded for traceability only: they are descriptive, never a gate. The
  strict physical gates (analytic L2, Ghia centreline oracle, cross-
  algorithm equivalence, true residuals) live inside the tests themselves.
- --with-pytest additionally runs the production body-force Python test
  (opt-in, requires a Python-enabled build); without the flag the report
  records it as not requested.

Usage:
  python3 scripts/n9_physical_campaign.py --build-dir build
  python3 scripts/n9_physical_campaign.py --build-dir build --with-pytest
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time
from pathlib import Path

CAMPAIGN_NAME = "n9-physical"

REQUIRED_TESTS = (
    "test_n9_physical_matrix",
    "test_n9_negative_diagnostics",
    "test_n9_s9_rhie_chow_consistency",
    "test_n9_s10_pressure_reference_policies",
    "test_n9_common_schur_qualification",
    "test_poiseuille_diagnostics",
    "test_ghia_cavity",
    "test_phase9_acceptance",
    "test_steady_incompressible_solver",
)
OPTIONAL_MPI_TESTS = ("test_n13_mpi_equivalence",)
PYTEST_TARGET = "tests/python/test_n9_poiseuille_body_force_production.py"
PYTEST_TEST_NAME = "pytest:" + PYTEST_TARGET

CATEGORIES = {
    "physical_matrix": ["test_n9_physical_matrix"],
    "negative_diagnostics": ["test_n9_negative_diagnostics"],
    "rhie_chow_consistency": ["test_n9_s9_rhie_chow_consistency"],
    "pressure_reference": ["test_n9_s10_pressure_reference_policies"],
    "schur_qualification": ["test_n9_common_schur_qualification"],
    "production_diagnostics": ["test_poiseuille_diagnostics", "test_ghia_cavity"],
    "steady_production": ["test_phase9_acceptance", "test_steady_incompressible_solver"],
    "mpi_equivalence": ["test_n13_mpi_equivalence"],
    "python_body_force": [PYTEST_TEST_NAME],
}

# Descriptive markers emitted by tests/validation/test_n9_physical_matrix.cpp
# on success (stderr markers: "ALG i/6 <ALGORITHM> <Case>",
# "STARTING COUPLED <model>", "COUPLED <model> <Case>",
# "EQUIVALENCE alg=<alg> couette=... poiseuille=... cavity=...").
# Reference counts derive from reading the test source at campaign-authoring
# time: 5 non-COUPLED algorithms enter the ALG loop (COUPLED is exercised
# separately through the 5 production Schur models PCD, LSC, BFBT, SIMPLE,
# SIMPLEC), and 4 cross-algorithm EQUIVALENCE rows are printed. The skew
# Couette campaign prints no per-case success marker: its evidence is the
# test PASS itself. These counts are informational traceability, never a gate.
DESCRIPTIVE_MARKERS = {
    "ALG": {"reference": 15, "what": "per-algorithm case entries"},
    "Couette": {"reference": 10, "what": "ALG + COUPLED Couette markers"},
    "Poiseuille": {"reference": 10, "what": "ALG + COUPLED Poiseuille markers"},
    "Cavity": {"reference": 5, "what": "ALG Cavity markers (Ghia Re=100)"},
    "STARTING COUPLED": {"reference": 5, "what": "per-Schur-model starts"},
    "EQUIVALENCE": {"reference": 4, "what": "cross-algorithm equivalence rows"},
}

POLICY = {
    "changes_numerical_tolerances": False,
    "disables_validation": False,
    "requires_exact_head_evidence": True,
    "stops_on_first_failed_gate": True,
    "descriptive_coverage_is_not_a_gate": True,
    "serial_evidence_is_mpi_qualification": False,
    "mpi_evidence_required_when_mpi_test_is_available": True,
    "pytest_body_force_is_opt_in": True,
}


def source_git_sha(source_dir: Path) -> str:
    result = subprocess.run(
        ["git", "-C", str(source_dir), "rev-parse", "HEAD"],
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    if result.returncode != 0:
        raise RuntimeError(
            f"Unable to determine campaign source revision: {result.stdout}"
        )
    sha = result.stdout.strip()
    if not re.fullmatch(r"[0-9a-f]{40}", sha):
        raise RuntimeError(f"Git returned an invalid full commit SHA: {sha!r}")
    return sha


def expected_sha_matches(actual_sha: str, expected_sha: str | None) -> bool:
    return not expected_sha or actual_sha == expected_sha


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
    pattern = re.compile(r"^\s*Test\s+#\s*\d+\s*:\s*(\S+)")
    return {
        match.group(1)
        for line in result.stdout.splitlines()
        if (match := pattern.match(line))
    }


def run_test(build_dir: Path, name: str) -> dict[str, object]:
    start = time.monotonic()
    result = ctest(
        build_dir,
        "-V",
        "--output-on-failure",
        "--no-tests=error",
        "--tests-regex",
        f"^{re.escape(name)}$",
    )
    elapsed = time.monotonic() - start
    no_test = "No tests were found" in result.stdout
    status = "PASS" if result.returncode == 0 and not no_test else "FAIL"
    return {
        "name": name,
        "status": status,
        "returncode": 1 if no_test else result.returncode,
        "elapsed_s": round(elapsed, 3),
        "output": result.stdout,
    }


def run_pytest(source_dir: Path) -> dict[str, object]:
    start = time.monotonic()
    result = subprocess.run(
        [sys.executable, "-m", "pytest", PYTEST_TARGET],
        cwd=source_dir,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    elapsed = time.monotonic() - start
    status = "PASS" if result.returncode == 0 else "FAIL"
    return {
        "name": PYTEST_TEST_NAME,
        "status": status,
        "returncode": result.returncode,
        "elapsed_s": round(elapsed, 3),
        "output": result.stdout,
    }


def build_category_results(
    results: list[dict[str, object]], required_tests: list[str]
) -> dict[str, dict[str, object]]:
    by_name = {str(item["name"]): item for item in results}
    required_set = set(required_tests)
    categories: dict[str, dict[str, object]] = {}
    for category, names in CATEGORIES.items():
        active_names = [name for name in names if name in required_set]
        if not active_names:
            categories[category] = {
                "required_tests": [],
                "completed": 0,
                "missing": [],
                "failed": [],
                "status": "SKIPPED",
            }
            continue
        missing = [name for name in active_names if name not in by_name]
        failed = [
            name
            for name in active_names
            if name in by_name and by_name[name]["status"] != "PASS"
        ]
        categories[category] = {
            "required_tests": active_names,
            "completed": len(active_names) - len(missing),
            "missing": missing,
            "failed": failed,
            "status": "PASS" if not missing and not failed else "FAIL",
        }
    return categories


def build_descriptive_coverage(
    results: list[dict[str, object]],
) -> dict[str, dict[str, object]]:
    matrix_output = next(
        (
            str(item["output"])
            for item in results
            if str(item["name"]) == "test_n9_physical_matrix"
        ),
        "",
    )
    coverage: dict[str, dict[str, object]] = {}
    for marker, meta in DESCRIPTIVE_MARKERS.items():
        observed = sum(
            1 for line in matrix_output.splitlines() if marker in line
        )
        coverage[marker] = {
            "what": meta["what"],
            "reference": meta["reference"],
            "observed": observed,
            "gating": False,
            "note": (
                "informational traceability only; the strict physical gates "
                "live inside the test and are never re-implemented here"
            ),
        }
    return coverage


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument(
        "--expected-git-sha",
        default=os.environ.get("GITHUB_SHA"),
        help="Expected full source SHA; defaults to GITHUB_SHA when set",
    )
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="JSON report path (default: <build-dir>/n9_physical_campaign.json)",
    )
    parser.add_argument(
        "--with-pytest",
        action="store_true",
        help=(
            "Additionally run the production body-force Python test "
            f"({PYTEST_TARGET}); requires a Python-enabled build"
        ),
    )
    args = parser.parse_args()

    build_dir = args.build_dir.resolve()
    if not (build_dir / "CTestTestfile.cmake").exists():
        print(
            f"error: {build_dir} is not a configured CMake build directory",
            file=sys.stderr,
        )
        return 2

    report_path = (
        args.report.resolve()
        if args.report is not None
        else build_dir / "n9_physical_campaign.json"
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)

    source_dir = Path(__file__).resolve().parents[1]
    try:
        git_sha = source_git_sha(source_dir)
    except RuntimeError as exc:
        report = {
            "campaign": CAMPAIGN_NAME,
            "status": "INCOMPLETE",
            "git_sha": None,
            "integrity_error": str(exc),
            "policy": POLICY,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2), file=sys.stderr)
        return 2

    expected_sha = args.expected_git_sha
    if not expected_sha_matches(git_sha, expected_sha):
        report = {
            "campaign": CAMPAIGN_NAME,
            "status": "INCOMPLETE",
            "git_sha": git_sha,
            "expected_git_sha": expected_sha,
            "integrity_error": "Campaign source SHA does not match expected GitHub SHA",
            "policy": POLICY,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2), file=sys.stderr)
        return 2

    available = discover_tests(build_dir)
    required_tests = list(REQUIRED_TESTS)
    mpi_tests = [name for name in OPTIONAL_MPI_TESTS if name in available]
    required_tests.extend(mpi_tests)
    if args.with_pytest:
        required_tests.append(PYTEST_TEST_NAME)
    missing = [name for name in required_tests if name not in available]
    if missing:
        report = {
            "campaign": CAMPAIGN_NAME,
            "status": "INCOMPLETE",
            "git_sha": git_sha,
            "expected_git_sha": expected_sha,
            "required_tests": required_tests,
            "missing_tests": missing,
            "policy": POLICY,
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 2

    results: list[dict[str, object]] = []
    for name in required_tests:
        if name == PYTEST_TEST_NAME:
            print(f"\n=== N9 campaign (pytest): {PYTEST_TARGET} ===", flush=True)
            result = run_pytest(source_dir)
        else:
            print(f"\n=== N9 pressure-velocity campaign: {name} ===", flush=True)
            result = run_test(build_dir, name)
        results.append(result)
        print(result["output"], end="", flush=True)
        if result["status"] != "PASS":
            break

    failed = [str(item["name"]) for item in results if item["status"] != "PASS"]
    categories = build_category_results(results, required_tests)
    descriptive_coverage = build_descriptive_coverage(results)
    pytest_status = (
        "not_requested"
        if not args.with_pytest
        else next(
            str(item["status"])
            for item in results
            if str(item["name"]) == PYTEST_TEST_NAME
        )
    )
    complete_execution = len(results) == len(required_tests) and not failed

    report = {
        "campaign": CAMPAIGN_NAME,
        "status": "PASS" if complete_execution else "FAIL",
        "git_sha": git_sha,
        "expected_git_sha": expected_sha,
        "required_tests": required_tests,
        "completed_tests": len(results),
        "failed_tests": failed,
        "results": results,
        "categories": categories,
        "descriptive_case_coverage": descriptive_coverage,
        "pytest_body_force": {
            "target": PYTEST_TARGET,
            "status": pytest_status,
            "opt_in": True,
        },
        "qualification_boundary": (
            "PASS means every declared campaign test executed and passed on "
            "the exact HEAD. It does not by itself close issue #601: closure "
            "additionally requires this JSON evidence attached to the issue "
            "and reviewed against the fiches methode 03 (Rhie-Chow) and 04 "
            "(pressure-velocity coupling)."
        ),
        "policy": POLICY,
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nN9 campaign report: {report_path}")
    print(f"N9 campaign status: {report['status']}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
