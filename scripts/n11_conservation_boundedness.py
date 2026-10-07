#!/usr/bin/env python3
"""Run the reproducible N11 conservation/boundedness evidence campaign.

The campaign reuses existing executable CFDX tests. It does not duplicate
solver or physics implementations and never changes tolerances. Each test is
run verbosely and its complete output is retained in a machine-readable JSON
report.

Usage:
  python3 scripts/n11_conservation_boundedness.py --build-dir build
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from pathlib import Path

REQUIRED_TESTS = (
    "test_conservation_boundedness",
    "test_transport_conservation",
    "test_n11_scalar_conservation",
    "test_conservation_assembly",
    "test_phase9_acceptance",
    "test_poiseuille_diagnostics",
    "test_ghia_cavity",
    "test_mms_scalar_diffusion",
    "test_cht_validation",
    "test_nonorthogonal_skew_campaign",
    "test_convection_polyhedral_campaign",
)

CATEGORIES = {
    "conservation_contract": [
        "test_conservation_boundedness",
        "test_transport_conservation",
        "test_conservation_assembly",
    ],
    "coupled_production": [
        "test_phase9_acceptance",
        "test_poiseuille_diagnostics",
        "test_ghia_cavity",
    ],
    "scalar_mms": ["test_mms_scalar_diffusion", "test_n11_scalar_conservation"],
    "energy": ["test_cht_validation"],
    "mesh_robustness": ["test_nonorthogonal_skew_campaign"],
    "scheme_boundedness": ["test_convection_polyhedral_campaign"],
}

QUANTITATIVE_EVIDENCE = {
    "test_conservation_boundedness": [r"(?mi)^.*(?:residual|balance|conservation).*$"],
    "test_transport_conservation": [r"(?mi)^.*(?:transport|conservation|balance).*$"],
    "test_conservation_assembly": [r"(?mi)^.*(?:assembly|conservation|balance).*$"],
    "test_phase9_acceptance": [r"(?m)^MODEL_SUMMARY\\b.*\\bsuccessful=\\d+", r"(?m)^ALGORITHM_INVARIANCE\\b.*"],
    "test_poiseuille_diagnostics": [r"(?m)^--- refinement level N=\\d+ ---$", r"(?m)^=== Refinement diagnostics ===$"],
    "test_ghia_cavity": [r"(?m)^GHIA Re=100 observed_order\\b.*", r"(?m)^GHIA_CAVITY_VALIDATION: PASS$"],
    "test_mms_scalar_diffusion": [r"(?mi)^.*(?:MMS|observed.order|L2|Linf).*$"],
    "test_n11_scalar_conservation": [r"(?m)^N11_SCALAR_QUALIFICATION: PASS\\b.*"],
    "test_cht_validation": [r"(?mi)^.*(?:CHT|heat|energy).*(?:PASS|residual|balance|flux).*$"],
    "test_nonorthogonal_skew_campaign": [r"(?m)^PHASE3_6\\s+skew=.*corrected_conservation=.*$"],
    "test_convection_polyhedral_campaign": [r"(?m)^N4_POLY_ORDER\\b.*order=.*$", r"(?m)^N4_POLY_BOUNDED\\b.*worst=.*$"],
}


KNOWN_SCOPE_GAPS = [
    {
        "id": "mpi-production-conservation",
        "status": "missing",
        "description": (
            "No dedicated production-physics N11 test currently executes the "
            "independent conservation audit under MPI; serial evidence must not "
            "be presented as MPI qualification."
        ),
    },
]


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
    quantitative = extract_quantitative_evidence(name, result.stdout)
    if status == "PASS" and not quantitative["complete"]:
        status = "FAIL"
    return {
        "name": name,
        "status": status,
        "returncode": 1 if no_test else result.returncode,
        "elapsed_s": round(elapsed, 3),
        "output": result.stdout,
        "quantitative_evidence": quantitative,
    }


def extract_quantitative_evidence(name: str, output: str) -> dict[str, object]:
    patterns = QUANTITATIVE_EVIDENCE.get(name, [])
    records: list[str] = []
    missing: list[str] = []
    for pattern in patterns:
        matches = re.findall(pattern, output)
        if matches:
            records.extend(matches)
        else:
            missing.append(pattern)
    return {
        "required_patterns": patterns,
        "records": records,
        "missing_patterns": missing,
        "complete": not missing,
    }


def build_category_results(results: list[dict[str, object]]) -> dict[str, dict[str, object]]:
    by_name = {str(item["name"]): item for item in results}
    categories: dict[str, dict[str, object]] = {}
    for category, names in CATEGORIES.items():
        missing = [name for name in names if name not in by_name]
        failed = [
            name for name in names
            if name in by_name and by_name[name]["status"] != "PASS"
        ]
        categories[category] = {
            "required_tests": names,
            "completed": len(names) - len(missing),
            "missing": missing,
            "failed": failed,
            "status": "PASS" if not missing and not failed else "FAIL",
        }
    return categories


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--build-dir", required=True, type=Path)
    parser.add_argument(
        "--report",
        type=Path,
        default=None,
        help="JSON report path (default: <build-dir>/n11_conservation_boundedness.json)",
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
        else build_dir / "n11_conservation_boundedness.json"
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)

    available = discover_tests(build_dir)
    missing = [name for name in REQUIRED_TESTS if name not in available]
    if missing:
        report = {
            "campaign": "N11 conservation and boundedness evidence",
            "status": "INCOMPLETE",
            "required_tests": list(REQUIRED_TESTS),
            "missing_tests": missing,
            "scope_gaps": KNOWN_SCOPE_GAPS,
            "policy": {
                "changes_numerical_tolerances": False,
                "disables_validation": False,
                "silent_clipping_or_repair": False,
                "serial_evidence_is_mpi_qualification": False,
            },
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 2

    results: list[dict[str, object]] = []
    for name in REQUIRED_TESTS:
        print(f"\n=== N11 conservation/boundedness: {name} ===", flush=True)
        result = run_test(build_dir, name)
        results.append(result)
        print(result["output"], end="", flush=True)
        if result["status"] != "PASS":
            break

    failed = [str(item["name"]) for item in results if item["status"] != "PASS"]
    categories = build_category_results(results)
    complete_execution = len(results) == len(REQUIRED_TESTS) and not failed

    report = {
        "campaign": "N11 conservation and boundedness evidence",
        "status": "PASS" if complete_execution else "FAIL",
        "required_tests": list(REQUIRED_TESTS),
        "completed_tests": len(results),
        "failed_tests": failed,
        "results": results,
        "categories": categories,
        "scope_gaps": KNOWN_SCOPE_GAPS,
        "qualification_boundary": (
            "PASS means every declared serial campaign test executed and passed, with its required quantitative evidence records present. "
            "It does not close N11 or imply MPI qualification. N11 closure still "
            "requires the complete declared conservation/boundedness population "
            "and its independent evidence."
        ),
        "policy": {
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_clipping_or_repair": False,
            "independent_conservation_over_linear_residual": True,
            "serial_evidence_is_mpi_qualification": False,
            "stops_on_first_failed_gate": True,
        },
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nN11 report: {report_path}")
    print(f"N11 campaign status: {report['status']}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
