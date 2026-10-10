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
import os
import re
import subprocess
import sys
import time
from pathlib import Path

QUANTITATIVE_RECORD_PREFIXES = (
    "MODEL_RESULT ",
    "GHIA ",
    "POISEUILLE_RESULT ",
    "N11_ENERGY_RESULT ",
    "N11_MPI_RESULT ",
)


def parse_quantitative_records(output: str) -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for line in output.splitlines():
        prefix = next(
            (p for p in QUANTITATIVE_RECORD_PREFIXES if line.startswith(p)),
            None,
        )
        if prefix is None:
            continue
        record: dict[str, object] = {"record_type": prefix.strip()}
        for token in line[len(prefix) :].split():
            if "=" not in token:
                continue
            key, value = token.split("=", 1)
            try:
                record[key] = float(value)
            except ValueError:
                record[key] = value
        records.append(record)
    return records


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
OPTIONAL_MPI_TESTS = ("test_n13_mpi_equivalence",)

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
    "mpi_conservation": ["test_n13_mpi_equivalence"],
}

KNOWN_SCOPE_GAPS = [
    {
        "id": "mpi-production-conservation",
        "status": "missing",
        "description": (
            "The active distributed MPI Poisson path now has an independent "
            "conservation gate, but the production incompressible/thermal "
            "physics solvers are not yet distributed under MPI; this PR must "
            "not be presented as production-physics MPI qualification."
        ),
    },
]


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
        help=(
            "JSON report path (default: <build-dir>/n11_conservation_boundedness.json)"
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
        else build_dir / "n11_conservation_boundedness.json"
    )
    report_path.parent.mkdir(parents=True, exist_ok=True)

    source_dir = Path(__file__).resolve().parents[1]
    try:
        git_sha = source_git_sha(source_dir)
    except RuntimeError as exc:
        report = {
            "campaign": "N11 conservation and boundedness evidence",
            "status": "INCOMPLETE",
            "git_sha": None,
            "integrity_error": str(exc),
            "policy": {
                "requires_exact_head_evidence": True,
                "changes_numerical_tolerances": False,
                "disables_validation": False,
            },
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2), file=sys.stderr)
        return 2

    expected_sha = args.expected_git_sha
    if not expected_sha_matches(git_sha, expected_sha):
        report = {
            "campaign": "N11 conservation and boundedness evidence",
            "status": "INCOMPLETE",
            "git_sha": git_sha,
            "expected_git_sha": expected_sha,
            "integrity_error": "Campaign source SHA does not match expected GitHub SHA",
            "policy": {
                "requires_exact_head_evidence": True,
                "changes_numerical_tolerances": False,
                "disables_validation": False,
            },
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2), file=sys.stderr)
        return 2

    available = discover_tests(build_dir)
    required_tests = list(REQUIRED_TESTS)
    mpi_tests = [name for name in OPTIONAL_MPI_TESTS if name in available]
    required_tests.extend(mpi_tests)
    missing = [name for name in required_tests if name not in available]
    if missing:
        report = {
            "campaign": "N11 conservation and boundedness evidence",
            "status": "INCOMPLETE",
            "git_sha": git_sha,
            "expected_git_sha": expected_sha,
            "required_tests": required_tests,
            "missing_tests": missing,
            "scope_gaps": KNOWN_SCOPE_GAPS,
            "policy": {
                "changes_numerical_tolerances": False,
                "disables_validation": False,
                "silent_clipping_or_repair": False,
                "serial_evidence_is_mpi_qualification": False,
                "mpi_evidence_required_when_mpi_test_is_available": True,
            },
        }
        report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(json.dumps(report, indent=2))
        return 2

    results: list[dict[str, object]] = []
    for name in required_tests:
        print(f"\n=== N11 conservation/boundedness: {name} ===", flush=True)
        result = run_test(build_dir, name)
        results.append(result)
        print(result["output"], end="", flush=True)
        if result["status"] != "PASS":
            break

    failed = [str(item["name"]) for item in results if item["status"] != "PASS"]
    categories = build_category_results(results, required_tests)
    quantitative_records = [
        record
        for item in results
        for record in parse_quantitative_records(str(item["output"]))
    ]
    quantitative_counts = {
        "MODEL_RESULT": sum(
            r["record_type"] == "MODEL_RESULT" for r in quantitative_records
        ),
        "GHIA": sum(r["record_type"] == "GHIA" for r in quantitative_records),
        "POISEUILLE_RESULT": sum(
            r["record_type"] == "POISEUILLE_RESULT" for r in quantitative_records
        ),
        "N11_ENERGY_RESULT": sum(
            r["record_type"] == "N11_ENERGY_RESULT" for r in quantitative_records
        ),
        "N11_MPI_RESULT": sum(
            r["record_type"] == "N11_MPI_RESULT" for r in quantitative_records
        ),
    }
    quantitative_requirements = {
        "MODEL_RESULT": 1,
        "GHIA": 4,
        "POISEUILLE_RESULT": 4,
        "N11_ENERGY_RESULT": 1,
        "N11_MPI_RESULT": 1 if "test_n13_mpi_equivalence" in required_tests else 0,
    }
    quantitative_coverage = {
        key: {
            "required": required,
            "observed": quantitative_counts[key],
            "complete": quantitative_counts[key] >= required,
        }
        for key, required in quantitative_requirements.items()
    }
    quantitative_complete = all(
        item["complete"] for item in quantitative_coverage.values()
    )
    complete_execution = (
        len(results) == len(required_tests) and not failed and quantitative_complete
    )

    report = {
        "campaign": "N11 conservation and boundedness evidence",
        "status": "PASS" if complete_execution else "FAIL",
        "git_sha": git_sha,
        "expected_git_sha": expected_sha,
        "required_tests": required_tests,
        "completed_tests": len(results),
        "failed_tests": failed,
        "results": results,
        "quantitative_records": quantitative_records,
        "quantitative_coverage": quantitative_coverage,
        "categories": categories,
        "scope_gaps": KNOWN_SCOPE_GAPS,
        "qualification_boundary": (
            "PASS means every declared serial campaign test executed and passed. "
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
            "mpi_evidence_required_when_mpi_test_is_available": True,
            "stops_on_first_failed_gate": True,
            "requires_exact_head_evidence": True,
        },
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nN11 report: {report_path}")
    print(f"N11 campaign status: {report['status']}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
