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
import re
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
    "test_fgmres_solver",
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
    # ctest right-aligns the test number in a variable-width field, so the
    # spacing before the '#' is not fixed: "Test   #1: name" for a single digit
    # and "Test #81: name" for two. Matching the literal prefix "Test #"
    # therefore silently skipped every test whose index had fewer digits than
    # the widest one, and the campaign reported those tests as missing rather
    # than as unparsed. Anchor on the structure instead.
    pattern = re.compile(r"^\s*Test\s+#\s*\d+\s*:\s*(\S+)")
    for line in result.stdout.splitlines():
        match = pattern.match(line)
        if match:
            names.add(match.group(1))
    return names


def run_test(build_dir: Path, name: str) -> dict[str, object]:
    start = time.monotonic()
    # -V is required, not --output-on-failure: CTest only echoes the stdout of
    # failing tests, so a passing campaign would capture no records at all and the
    # structured evidence section would be empty for every successful run.
    result = ctest(build_dir, "-V", "--output-on-failure", "--tests-regex", f"^{name}$")
    elapsed = time.monotonic() - start
    return {
        "name": name,
        "returncode": result.returncode,
        "status": "PASS" if result.returncode == 0 else "FAIL",
        "elapsed_s": round(elapsed, 3),
        "output": result.stdout,
    }


_KEY = re.compile(r"[\w.|]+")
_CTEST_LINE_PREFIX = re.compile(r"^\d+: ?")


def _parse_value(value: str) -> object:
    if value in {"true", "false"}:
        return value == "true"
    try:
        if any(ch in value for ch in ".eE"):
            return float(value)
        return int(value)
    except ValueError:
        return value


def parse_key_value_records(output: str, prefix: str) -> list[dict[str, object]]:
    """Extract structured key=value records emitted by an existing validation test.

    The first token is recorded as ``model`` when it carries no ``=``. Remaining
    bare tokens are recorded as ``gates`` so that a failure record keeps the gate
    names the test reported instead of only the model that failed.
    """
    records: list[dict[str, object]] = []
    pattern = re.compile(r"^" + re.escape(prefix) + r"\s+(.*)$")
    for raw_line in output.splitlines():
        # CTest prefixes every captured line with the test index ("21: ...").
        line = _CTEST_LINE_PREFIX.sub("", raw_line, count=1)
        match = pattern.match(line)
        if not match:
            continue
        tokens = match.group(1).split()
        record: dict[str, object] = {}
        if tokens and "=" not in tokens[0]:
            record["model"] = tokens[0]
            tokens = tokens[1:]
        gates: list[str] = []
        position = 0
        while position < len(tokens):
            token = tokens[position]
            key, separator, value = token.partition("=")
            position += 1
            if not separator or not _KEY.fullmatch(key):
                gates.append(token)
                continue
            parsed = _parse_value(value)
            trailing = tokens[position:]
            if isinstance(parsed, str) and trailing and not any("=" in item for item in trailing):
                parsed = " ".join([value, *trailing])
                position = len(tokens)
            record[key] = parsed
        if gates:
            record["gates"] = gates
        if record:
            records.append(record)
    return records


def extract_n8_evidence(results: list[dict[str, object]]) -> dict[str, object]:
    """Normalize structured evidence already emitted by N8 validation tests.

    The source tests remain the numerical oracles. This function only parses
    their existing stdout; it does not invent gates, alter tolerances, or
    reinterpret failures as passes.
    """
    evidence: dict[str, object] = {
        "physical_model_results": [],
        "physical_model_failures": [],
        "schur_quantitative": [],
        "schur_production": [],
    }

    for result in results:
        output = str(result["output"])
        evidence["physical_model_results"].extend(
            parse_key_value_records(output, "MODEL_RESULT")
        )
        evidence["physical_model_failures"].extend(
            parse_key_value_records(output, "MODEL_FAILURES")
        )
        evidence["schur_quantitative"].extend(
            parse_key_value_records(output, "N8_SCHUR")
        )
        evidence["schur_production"].extend(
            parse_key_value_records(output, "n8_schur_benchmark")
        )
        evidence["schur_production"].extend(
            parse_key_value_records(output, "n8_schur_benchmark_lifecycle")
        )

    return evidence


def audit_evidence_coverage(evidence: dict[str, object]) -> dict[str, object]:
    """Report whether structured N8 evidence contains the required fields."""
    requirements = {
        "physical_model_results": (
            "model", "solver_converged", "iterations", "gates_failed"
        ),
        "schur_quantitative": (
            "case", "cond_inf_Auu", "exact_solve_backward_error", "machine_epsilon"
        ),
        "schur_production": (
            "cells", "unknowns", "nnz", "schur_nnz", "true_residual",
            "iterations", "setup_us", "solve_us", "pressure_coarse_size",
            "hierarchy_builds", "numeric_updates",
        ),
        "schur_lifecycle": (
            "cells", "coefficient_changed", "hierarchy_builds_before",
            "hierarchy_builds_after", "numeric_updates",
            "updated_iterations", "updated_true_residual",
            "graph_change_rebuild",
        ),
    }
    missing: dict[str, list[dict[str, object]]] = {}
    checked = complete = 0

    for category, fields in requirements.items():
        if category == "schur_lifecycle":
            records = [
                record for record in evidence.get("schur_production", [])
                if isinstance(record, dict) and "coefficient_changed" in record
            ]
        elif category == "schur_production":
            records = [
                record for record in evidence.get("schur_production", [])
                if isinstance(record, dict) and "coefficient_changed" not in record
            ]
        else:
            records = evidence.get(category, [])
        if not isinstance(records, list):
            records = []

        category_missing: list[dict[str, object]] = []
        for index, record in enumerate(records):
            if not isinstance(record, dict):
                category_missing.append({"index": index, "fields": list(fields)})
                continue
            checked += 1
            absent = [field for field in fields if field not in record]
            if absent:
                category_missing.append({"index": index, "fields": absent})
            else:
                complete += 1
        if category_missing:
            missing[category] = category_missing

    return {
        "status": "COMPLETE" if not missing else "INCOMPLETE",
        "records_checked": checked,
        "records_complete": complete,
        "missing_fields": missing,
        "policy": "diagnostic_only",
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
    evidence = extract_n8_evidence(results)
    evidence_coverage = audit_evidence_coverage(evidence)
    report = {
        "campaign": "N8 complete solver/preconditioner qualification",
        "status": "PASS" if not failed and completed == len(REQUIRED_TESTS) else "FAIL",
        "required_tests": list(REQUIRED_TESTS),
        "completed_tests": completed,
        "failed_tests": failed,
        "results": results,
        "evidence": evidence,
        "evidence_coverage": evidence_coverage,
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
