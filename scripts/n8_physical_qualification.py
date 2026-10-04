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
        "run_configuration": [],
        "model_summary": [],
        "resolved_plans": [],
        "physical_model_results": [],
        "physical_model_failures": [],
        "schur_quantitative": [],
        "schur_production": [],
    }

    for result in results:
        output = str(result["output"])
        evidence["run_configuration"].extend(
            parse_key_value_records(output, "MODEL_CONFIG")
        )
        evidence["model_summary"].extend(
            parse_key_value_records(output, "MODEL_SUMMARY")
        )
        evidence["resolved_plans"].extend(
            parse_key_value_records(output, "MODEL_PLAN")
        )
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


# Required evidence fields, as (category, source list, discriminator field,
# whether the discriminator is expected, required fields). The quantitative
# Schur test and the production benchmark each print two record shapes, so a
# discriminator keeps the requirement set aligned with the emitted layout.
_EVIDENCE_REQUIREMENTS: tuple[tuple[str, str, str | None, bool | None, tuple[str, ...]], ...] = (
    (
        "run_configuration",
        "run_configuration",
        None,
        None,
        (
            "algorithm",
            "nx",
            "ny",
            "bounded",
            "coupled_requested_krylov",
            "preconditioner",
            "preconditioner_id",
            "pressure_requested_krylov",
            "pressure_requested_preconditioner",
            "alpha_u",
            "alpha_p",
            "pressure_correctors",
            "fractional_steps",
        ),
    ),
    (
        "resolved_plans",
        "resolved_plans",
        None,
        None,
        (
            "algorithm",
            "coupled_resolved",
            "coupled_krylov",
            "coupled_preconditioner",
            "pressure_resolved",
            "pressure_krylov",
            "pressure_preconditioner",
        ),
    ),
    (
        "physical_model_results",
        "physical_model_results",
        None,
        None,
        ("model", "solver_converged", "iterations", "gates_failed"),
    ),
    (
        "schur_quantitative_oracle",
        "schur_quantitative",
        "method",
        False,
        (
            "case",
            "cond_inf_Auu",
            "exact_oracle_discrepancy",
            "exact_solve_backward_error",
            "machine_epsilon",
        ),
    ),
    (
        "schur_quantitative_method",
        "schur_quantitative",
        "method",
        True,
        ("case", "method", "algebra_error", "exact_schur_error"),
    ),
    (
        "schur_production",
        "schur_production",
        "coefficient_changed",
        False,
        (
            "cells",
            "unknowns",
            "nnz",
            "schur_nnz",
            "true_residual",
            "iterations",
            "setup_us",
            "solve_us",
            "pressure_coarse_size",
            "hierarchy_builds",
            "numeric_updates",
        ),
    ),
    (
        "schur_lifecycle",
        "schur_production",
        "coefficient_changed",
        True,
        (
            "cells",
            "coefficient_changed",
            "hierarchy_builds_before",
            "hierarchy_builds_after",
            "numeric_updates",
            "updated_iterations",
            "updated_true_residual",
            "graph_change_rebuild",
        ),
    ),
)


def audit_evidence_coverage(evidence: dict[str, object]) -> dict[str, object]:
    """Report whether structured N8 evidence contains the required fields.

    The audit is descriptive: it never recomputes a gate, a tolerance, or a
    verdict. One source list can hold two record shapes, so each category is
    selected by a discriminator field instead of assuming a single layout. A
    category that produced no record is reported as well, because an empty
    evidence list is not evidence of completeness.
    """
    sources: dict[str, list[object]] = {}
    malformed_records: dict[str, list[int]] = {}
    for _, source, _, _, _ in _EVIDENCE_REQUIREMENTS:
        if source in sources:
            continue
        records = evidence.get(source, [])
        records = records if isinstance(records, list) else []
        sources[source] = records
        # Two categories can share one source list, so malformed entries are
        # collected once per source instead of once per category.
        malformed = [
            index for index, record in enumerate(records) if not isinstance(record, dict)
        ]
        if malformed:
            malformed_records[source] = malformed

    missing: dict[str, list[dict[str, object]]] = {}
    categories_without_records: list[str] = []
    checked = complete = 0

    for category, source, discriminator, expected, fields in _EVIDENCE_REQUIREMENTS:
        selected: list[tuple[int, dict[str, object]]] = []
        for index, record in enumerate(sources[source]):
            if not isinstance(record, dict):
                continue
            if discriminator is not None and (discriminator in record) is not expected:
                continue
            selected.append((index, record))
        if not selected:
            categories_without_records.append(category)

        category_missing: list[dict[str, object]] = []
        for source_index, record in selected:
            checked += 1
            absent = [field for field in fields if field not in record]
            if absent:
                category_missing.append({"index": source_index, "fields": absent})
            else:
                complete += 1
        if category_missing:
            missing[category] = category_missing

    incomplete = bool(missing or categories_without_records or malformed_records)
    return {
        "status": "INCOMPLETE" if incomplete else "COMPLETE",
        "records_checked": checked,
        "records_complete": complete,
        "missing_fields": missing,
        "categories_without_records": categories_without_records,
        "malformed_records": malformed_records,
        "policy": "diagnostic_only",
    }


def _model_names(records: list[object], key: str) -> list[str]:
    names: list[str] = []
    for record in records:
        if isinstance(record, dict):
            value = record.get(key)
            if isinstance(value, str) and value not in names:
                names.append(value)
    return names


def audit_model_resolution(evidence: dict[str, object]) -> dict[str, object]:
    """Check that every configured physical model produced a verdict record.

    The acceptance test announces the solver and preconditioner configuration of
    every model before running it. Pairing those announcements with the result
    and failure records makes a dropped, replaced, or unannounced run visible
    instead of silently shortening the campaign.
    """

    def records(source: str) -> list[object]:
        value = evidence.get(source, [])
        return value if isinstance(value, list) else []

    configured = _model_names(records("run_configuration"), "algorithm")
    resolved = _model_names(records("physical_model_results"), "model")
    failed_models = _model_names(records("physical_model_failures"), "model")
    # A model that fails still emits a result record, so the two lists overlap.
    for name in failed_models:
        if name not in resolved:
            resolved.append(name)
    unresolved = [name for name in configured if name not in resolved]
    undeclared = [name for name in resolved if name not in configured]

    summaries = [record for record in records("model_summary") if isinstance(record, dict)]
    summary = summaries[0] if len(summaries) == 1 else None
    reported_successful = summary.get("successful") if summary else None
    reported_failed = summary.get("failed") if summary else None
    reported_total = None
    if isinstance(reported_successful, int) and isinstance(reported_failed, int):
        reported_total = reported_successful + reported_failed
    observed_total = len(resolved)
    observed_failed = len(failed_models)
    # The acceptance test tallies its own model runs. A model emits a result
    # record whether it passes or fails, so the total number of result records
    # and the number of failure records are what must match the printed tally;
    # anything else means the capture is not the complete set of model runs.
    if reported_total is None or reported_failed is None:
        tallies_agree = None
    else:
        tallies_agree = (
            reported_total == observed_total and reported_failed == observed_failed
        )

    complete = not unresolved and not undeclared and tallies_agree is not False
    return {
        "status": "COMPLETE" if complete else "INCOMPLETE",
        "models_configured": len(configured),
        "models_resolved": len([name for name in resolved if name in configured]),
        "unresolved_models": unresolved,
        "undeclared_models": undeclared,
        "reported_tally": (
            None
            if summary is None
            else {
                "successful": reported_successful,
                "failed": reported_failed,
                "total": reported_total,
            }
        ),
        "observed_tally": {
            "resolved": observed_total,
            "failed": observed_failed,
            "total": observed_total,
        },
        "tallies_agree": tallies_agree,
        "policy": "diagnostic_only",
    }


_METHODS = (
    # (sub-problem, method, requested field in the run configuration, resolved
    #  flag, resolved field). The request is read only from the record emitted
    # before the run; the resolved-plan record is never a source of truth for
    # it. Krylov and preconditioner are gated symmetrically: a substituted
    # Krylov method is as much a substitution as a substituted preconditioner.
    (
        "coupled",
        "krylov",
        "coupled_requested_krylov",
        "coupled_resolved",
        "coupled_krylov",
    ),
    (
        "coupled",
        "preconditioner",
        "preconditioner",
        "coupled_resolved",
        "coupled_preconditioner",
    ),
    (
        "pressure",
        "krylov",
        "pressure_requested_krylov",
        "pressure_resolved",
        "pressure_krylov",
    ),
    (
        "pressure",
        "preconditioner",
        "pressure_requested_preconditioner",
        "pressure_resolved",
        "pressure_preconditioner",
    ),
)

# A coupled algorithm solves the 4N system and never runs the segregated
# pressure sub-problem; a segmented one does the opposite. The declared
# algorithm taxonomy is what the report checks the resolved structure against.
_COUPLED_ALGORITHM_PREFIX = "COUPLED/"


def expected_resolved_structure(algorithm: str) -> dict[str, bool]:
    """Return which sub-problems a declared algorithm is expected to resolve."""
    coupled = algorithm.startswith(_COUPLED_ALGORITHM_PREFIX)
    return {"coupled": coupled, "pressure": not coupled}


def audit_linear_plan(evidence: dict[str, object]) -> dict[str, object]:
    """Compare the announced requests with the plans the dispatcher resolved.

    The requests come from the records emitted before each run and the resolved
    plans from the records emitted after it, so a resolution cannot be checked
    against a request it also reports itself. An automatic request may
    legitimately resolve to another method, and that resolution is recorded
    instead of being invisible. An explicit request must be the method that ran,
    for the Krylov method as well as for the preconditioner: the dispatcher
    rejects an incompatible explicit request, so a mismatch means the report and
    the run disagree and must not be read as a pass.
    """

    def records(source: str) -> list[object]:
        value = evidence.get(source, [])
        return value if isinstance(value, list) else []

    configurations: dict[str, dict[str, object]] = {}
    for record in records("run_configuration"):
        if isinstance(record, dict) and isinstance(record.get("algorithm"), str):
            configurations.setdefault(str(record["algorithm"]), record)

    resolved_plans: dict[str, dict[str, object]] = {}
    for record in records("resolved_plans"):
        if isinstance(record, dict) and isinstance(record.get("algorithm"), str):
            resolved_plans[str(record["algorithm"])] = record

    automatic_resolutions: list[dict[str, str]] = []
    substitutions: list[dict[str, str]] = []
    structure_mismatches: list[dict[str, object]] = []
    models_without_plan: list[str] = []
    methods_compared = 0
    checked_structures: set[tuple[str, str]] = set()

    for algorithm, configuration in configurations.items():
        plan = resolved_plans.get(algorithm)
        if plan is None:
            models_without_plan.append(algorithm)
            continue

        expected = expected_resolved_structure(algorithm)
        for subproblem, method, request_field, resolved_flag, resolved_field in _METHODS:
            resolved_here = plan.get(resolved_flag) is True
            # The resolved flag is shared by the two methods of a sub-problem,
            # so the structure is checked once per sub-problem.
            if (algorithm, subproblem) not in checked_structures:
                checked_structures.add((algorithm, subproblem))
                if resolved_here != expected[subproblem]:
                    structure_mismatches.append(
                        {
                            "algorithm": algorithm,
                            "subproblem": subproblem,
                            "expected_resolved": expected[subproblem],
                            "resolved": resolved_here,
                        }
                    )
            if not resolved_here:
                continue
            resolved = str(plan.get(resolved_field, ""))
            asked = configuration.get(request_field)
            if not isinstance(asked, str) or asked == "":
                # A resolution that cannot be attributed to an announced request
                # is not evidence of anything.
                substitutions.append(
                    {
                        "algorithm": algorithm,
                        "subproblem": subproblem,
                        "method": method,
                        "requested": "undeclared",
                        "resolved": resolved,
                    }
                )
                continue
            methods_compared += 1
            if asked == "auto":
                automatic_resolutions.append(
                    {
                        "algorithm": algorithm,
                        "subproblem": subproblem,
                        "method": method,
                        "requested": asked,
                        "resolved": resolved,
                    }
                )
            elif asked != resolved:
                substitutions.append(
                    {
                        "algorithm": algorithm,
                        "subproblem": subproblem,
                        "method": method,
                        "requested": asked,
                        "resolved": resolved,
                    }
                )

    return {
        "status": (
            "VIOLATION"
            if substitutions or models_without_plan or structure_mismatches
            else "COMPLETE"
        ),
        "methods_compared": methods_compared,
        "automatic_resolutions": automatic_resolutions,
        "substitutions": substitutions,
        "structure_mismatches": structure_mismatches,
        "models_without_plan": models_without_plan,
        "policy": "explicit_linear_request_must_be_honored",
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
    model_resolution = audit_model_resolution(evidence)
    linear_plan = audit_linear_plan(evidence)
    # The one enforced gate of this report: an explicitly requested linear
    # method must be the method that ran. It is a contract, not a tolerance, so
    # it can fail the campaign. Everything else here stays diagnostic.
    gate_violations = (
        linear_plan["substitutions"]
        or linear_plan["models_without_plan"]
        or linear_plan["structure_mismatches"]
    )
    report = {
        "campaign": "N8 complete solver/preconditioner qualification",
        "status": (
            "FAIL"
            if failed or completed != len(REQUIRED_TESTS) or gate_violations
            else "PASS"
        ),
        "required_tests": list(REQUIRED_TESTS),
        "completed_tests": completed,
        "failed_tests": failed,
        "results": results,
        "evidence": evidence,
        "evidence_coverage": evidence_coverage,
        "model_resolution": model_resolution,
        "linear_plan": linear_plan,
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
            "enforced_gates": [
                "every explicitly requested linear method is the method that ran",
            ],
        },
    }
    report_path.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(f"\nN8 qualification report: {report_path}")
    print(f"N8 qualification status: {report['status']}")
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
