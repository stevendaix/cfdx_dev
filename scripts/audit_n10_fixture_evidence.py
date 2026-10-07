#!/usr/bin/env python3
"""Audit the verified N10 public-fixture evidence chain.

The gate combines independent byte verification, production importer
qualification and the persisted production-solver smoke evidence for the
fixture that has an explicit solver-execution contract. It does not claim
physical validation or qualification.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def load(path: Path) -> dict[str, object]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SystemExit(f"error: cannot read valid JSON report {path}: {exc}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--verification", type=Path, required=True)
    parser.add_argument("--qualification", type=Path, required=True)
    parser.add_argument("--solver-evidence", type=Path, required=True)
    parser.add_argument("--numerical-evidence", type=Path, required=True)
    parser.add_argument("--applicability", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    verification = load(args.verification)
    qualification = load(args.qualification)
    solver_evidence = load(args.solver_evidence)
    numerical_evidence = load(args.numerical_evidence)
    applicability = load(args.applicability)

    verification_fixture_rows = [
        row
        for row in verification.get("fixtures", [])
        if isinstance(row, dict) and isinstance(row.get("id"), str)
    ]
    qualification_fixture_rows = [
        row
        for row in qualification.get("fixtures", [])
        if isinstance(row, dict) and isinstance(row.get("id"), str)
    ]
    verification_ids_list = [str(row["id"]) for row in verification_fixture_rows]
    qualification_ids_list = [str(row["id"]) for row in qualification_fixture_rows]
    verification_rows = {str(row["id"]): row for row in verification_fixture_rows}
    qualification_rows = {str(row["id"]): row for row in qualification_fixture_rows}

    solver_fixture = str(solver_evidence.get("fixture", ""))
    solver_status = solver_evidence.get("execution_status")
    solver_support = solver_evidence.get("numerical_support")
    solver_physical_validation = solver_evidence.get("physical_validation")
    solver_restart = solver_evidence.get("restart_artifact")
    quantitative = solver_evidence.get("quantitative_solver_evidence")
    numerical_fixture = str(numerical_evidence.get("fixture", ""))
    numerical_status = numerical_evidence.get("status")
    numerical_cells_total = numerical_evidence.get("cells_total")
    numerical_cells_checked = numerical_evidence.get("cells_checked")
    numerical_boundary_excluded = numerical_evidence.get("boundary_cells_excluded")
    numerical_constant_linf = numerical_evidence.get("constant_gradient_linf")
    numerical_linear_linf = numerical_evidence.get("linear_gradient_linf")
    numerical_boundary_method = numerical_evidence.get("linear_boundary_reconstruction")
    numerical_gradient_method = numerical_evidence.get("linear_gradient_method")
    numerical_physical_validation = numerical_evidence.get("physical_validation")
    numerical_observed_order = numerical_evidence.get("observed_order")
    numerical_tolerances_changed = numerical_evidence.get("numerical_tolerances_changed")
    applicability_matrix = applicability.get("matrix")

    expected_fixture_ids = {
        "meshio-su2-square",
        "meshio-gmsh-insulated-2-2",
        "meshio-vtk-unstructured",
        "openfoam-airfoil2d",
    }
    verification_ids = set(verification_rows)
    qualification_ids = set(qualification_rows)
    mismatches: list[str] = []
    if len(verification_ids_list) != len(verification_ids):
        mismatches.append("independent verification contains duplicate fixture ids")
    if len(qualification_ids_list) != len(qualification_ids):
        mismatches.append("production qualification contains duplicate fixture ids")
    if verification_ids != expected_fixture_ids:
        mismatches.append(
            "independent verification fixture set does not match the required N10 public fixtures"
        )
    if qualification_ids != expected_fixture_ids:
        mismatches.append(
            "production qualification fixture set does not match the required N10 public fixtures"
        )

    expected_applicability = {
        "meshio-su2-square": "NOT_APPLICABLE_TO_CURRENT_3D_SOLVER",
        "meshio-gmsh-insulated-2-2": "NUMERICAL_APPLICABILITY_NOT_CLAIMED",
        "meshio-vtk-unstructured": "APPLICABLE_FOR_3D_EXECUTION_SMOKE",
        "openfoam-airfoil2d": "NUMERICAL_APPLICABILITY_NOT_CLAIMED",
    }

    ids = sorted(verification_ids | qualification_ids)
    rows: list[dict[str, object]] = []

    for fixture_id in ids:
        verified = verification_rows.get(fixture_id)
        qualified = qualification_rows.get(fixture_id)
        if verified is None:
            mismatches.append(f"{fixture_id}: missing from independent verification")
            continue
        if qualified is None:
            mismatches.append(f"{fixture_id}: missing from production qualification")
            continue

        byte_status = verified.get("status")
        imported = qualified.get("imported")
        topology = qualified.get("topology_valid")
        geometry = qualified.get("geometry_quality_valid")
        counts = qualified.get("counts", {})

        if byte_status != "VERIFIED":
            mismatches.append(f"{fixture_id}: byte verification status={byte_status}")
        if imported is not True:
            mismatches.append(f"{fixture_id}: imported={imported}")
        if topology is not True:
            mismatches.append(f"{fixture_id}: topology_valid={topology}")
        if geometry is not True:
            mismatches.append(f"{fixture_id}: geometry_quality_valid={geometry}")
        if not isinstance(counts, dict) or any(
            not isinstance(counts.get(key), int) or counts.get(key, 0) <= 0
            for key in ("points", "faces", "cells")
        ):
            mismatches.append(f"{fixture_id}: invalid positive mesh counts")

        solver_row: dict[str, object] = {
            "execution_status": "NOT_CLAIMED",
            "numerical_support": "NOT_CLAIMED",
            "physical_validation": "NOT_CLAIMED",
            "restart_artifact": False,
        }
        if fixture_id == solver_fixture:
            solver_row = {
                "execution_status": solver_status,
                "numerical_support": solver_support,
                "physical_validation": solver_physical_validation,
                "restart_artifact": solver_restart,
                "quantitative_solver_evidence": quantitative,
            }

        rows.append(
            {
                "id": fixture_id,
                "byte_verified": byte_status == "VERIFIED",
                "imported": imported is True,
                "topology_valid": topology is True,
                "geometry_quality_valid": geometry is True,
                "counts": counts,
                "solver_evidence": solver_row,
                "numerical_applicability": expected_applicability.get(
                    fixture_id, "NOT_CLAIMED"
                ),
            }
        )

    if solver_fixture != "meshio-vtk-unstructured":
        mismatches.append(
            "solver evidence must target the verified meshio-vtk-unstructured fixture"
        )
    if solver_status != "PASS":
        mismatches.append(f"solver evidence execution_status={solver_status}")
    if solver_support != "execution_and_convergence_smoke_only":
        mismatches.append(f"solver evidence numerical_support={solver_support}")
    if solver_physical_validation != "NOT_CLAIMED":
        mismatches.append(
            f"solver evidence physical_validation={solver_physical_validation}"
        )
    if solver_restart is not True:
        mismatches.append(f"solver evidence restart_artifact={solver_restart}")

    if numerical_fixture != "meshio-vtk-unstructured":
        mismatches.append("numerical evidence must target the verified meshio-vtk-unstructured fixture")
    if numerical_status != "PASS":
        mismatches.append(f"numerical evidence status={numerical_status}")
    if not isinstance(numerical_cells_total, int) or numerical_cells_total <= 0:
        mismatches.append(f"numerical evidence cells_total={numerical_cells_total}")
    if not isinstance(numerical_cells_checked, int) or numerical_cells_checked <= 0:
        mismatches.append(f"numerical evidence cells_checked={numerical_cells_checked}")
    elif numerical_cells_total != numerical_cells_checked:
        mismatches.append(f"numerical evidence cells_checked={numerical_cells_checked} differs from cells_total={numerical_cells_total}")
    if numerical_boundary_excluded != 0:
        mismatches.append(f"numerical evidence boundary_cells_excluded={numerical_boundary_excluded}")
    if not isinstance(numerical_constant_linf, (int, float)) or not math.isfinite(float(numerical_constant_linf)):
        mismatches.append(f"numerical evidence constant_gradient_linf={numerical_constant_linf}")
    elif float(numerical_constant_linf) > 1e-12:
        mismatches.append(f"numerical evidence constant_gradient_linf={numerical_constant_linf} > 1e-12")
    if not isinstance(numerical_linear_linf, (int, float)) or not math.isfinite(float(numerical_linear_linf)):
        mismatches.append(f"numerical evidence linear_gradient_linf={numerical_linear_linf}")
    elif float(numerical_linear_linf) > 1e-9:
        mismatches.append(f"numerical evidence linear_gradient_linf={numerical_linear_linf} > 1e-9")
    if numerical_boundary_method != "explicit_neumann":
        mismatches.append(f"numerical evidence linear_boundary_reconstruction={numerical_boundary_method}")
    if numerical_gradient_method != "weighted_least_squares":
        mismatches.append(f"numerical evidence linear_gradient_method={numerical_gradient_method}")
    if numerical_physical_validation != "NOT_CLAIMED":
        mismatches.append(f"numerical evidence physical_validation={numerical_physical_validation}")
    if numerical_observed_order != "NOT_CLAIMED":
        mismatches.append(f"numerical evidence observed_order={numerical_observed_order}")
    if numerical_tolerances_changed is not False:
        mismatches.append(f"numerical evidence numerical_tolerances_changed={numerical_tolerances_changed}")
    if not isinstance(applicability_matrix, list):
        mismatches.append("applicability matrix is missing")
    else:
        applicability_ids = [
            str(row.get("id"))
            for row in applicability_matrix
            if isinstance(row, dict) and isinstance(row.get("id"), str)
        ]
        applicability_rows = {
            fixture_id: row
            for fixture_id, row in (
                (str(row.get("id")), row)
                for row in applicability_matrix
                if isinstance(row, dict) and isinstance(row.get("id"), str)
            )
        }
        if len(applicability_ids) != len(set(applicability_ids)):
            mismatches.append("applicability matrix contains duplicate fixture ids")
        if set(applicability_ids) != set(expected_applicability):
            mismatches.append(
                "applicability matrix fixture set does not match the required N10 public fixtures"
            )
        if set(applicability_ids) != verification_ids or set(applicability_ids) != qualification_ids:
            mismatches.append(
                "applicability matrix fixture set does not match verification and production qualification fixture sets"
            )
        for fixture_id, expected_status in expected_applicability.items():
            row = applicability_rows.get(fixture_id)
            if row is None:
                continue
            if row.get("numerical_applicability") != expected_status:
                mismatches.append(
                    f"{fixture_id}: numerical_applicability={row.get('numerical_applicability')}"
                )
    applicability_policy = applicability.get("policy")
    if not isinstance(applicability_policy, dict):
        mismatches.append("applicability policy is missing")
    else:
        required_policy = {
            "import_success_is_not_numerical_qualification": True,
            "physical_validation_claimed": False,
            "forces_unsupported_2d_fixtures_through_3d_solver": False,
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_fallbacks": False,
        }
        for key, expected in required_policy.items():
            if applicability_policy.get(key) is not expected:
                mismatches.append(f"applicability policy {key}={applicability_policy.get(key)!r}")

    if not isinstance(quantitative, dict):
        mismatches.append("solver evidence quantitative_solver_evidence is missing")
    else:
        trace_count = quantitative.get("trace_count")
        final_iteration = quantitative.get("final_iteration")
        iterations = quantitative.get("iterations")
        if not isinstance(trace_count, int) or trace_count <= 0:
            mismatches.append(f"solver evidence trace_count={trace_count}")
        if not isinstance(final_iteration, int) or final_iteration <= 0:
            mismatches.append(f"solver evidence final_iteration={final_iteration}")
        if not isinstance(iterations, list) or not iterations:
            mismatches.append("solver evidence iterations is empty")
        elif isinstance(trace_count, int) and trace_count != len(iterations):
            mismatches.append(
                f"solver evidence trace_count={trace_count} differs from iterations={len(iterations)}"
            )
        required_metrics = (
            "max_momentum_residual_relative",
            "max_continuity_normalized",
            "final_momentum_residual_relative",
            "final_continuity_normalized",
            "final_velocity_change_inf",
            "final_pressure_change_inf",
            "final_flux_velocity_mismatch_linf",
        )
        for metric in required_metrics:
            value = quantitative.get(metric)
            if not isinstance(value, (int, float)) or not math.isfinite(float(value)):
                mismatches.append(f"solver evidence {metric}={value}")

    report = {
        "campaign": "N10 verified public-fixture evidence gate",
        "status": "PASS" if ids and not mismatches else "FAIL",
        "fixture_count": len(ids),
        "fixtures": rows,
        "mismatches": mismatches,
        "scope": {
            "byte_integrity": "independent SHA-256 verification",
            "production_import": "import_mesh dispatcher",
            "topology": "Mesh::topo_validate",
            "geometry": "existing 3D validator or explicit 2D edge-mesh validation",
            "numerical_solver_support": (
                "production execution/convergence smoke evidence for "
                "meshio-vtk-unstructured only"
            ),
            "physical_validation": "NOT_CLAIMED",
            "numerical_applicability": "persisted fixture-specific matrix",
            "acceptance_matrix": {
                "fixture_set": sorted(expected_fixture_ids),
                "verification_fixture_set": sorted(verification_ids),
                "qualification_fixture_set": sorted(qualification_ids),
                "applicability_fixture_set": sorted(applicability_ids)
                if isinstance(applicability_matrix, list)
                else [],
            },
        },
        "policy": {
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_fallbacks": False,
            "solver_qualification_deferred": True,
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
