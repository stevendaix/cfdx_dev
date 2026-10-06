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
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()

    verification = load(args.verification)
    qualification = load(args.qualification)
    solver_evidence = load(args.solver_evidence)

    verification_rows = {
        str(row["id"]): row
        for row in verification.get("fixtures", [])
        if isinstance(row, dict) and "id" in row
    }
    qualification_rows = {
        str(row["id"]): row
        for row in qualification.get("fixtures", [])
        if isinstance(row, dict) and "id" in row
    }

    solver_fixture = str(solver_evidence.get("fixture", ""))
    solver_status = solver_evidence.get("execution_status")
    solver_support = solver_evidence.get("numerical_support")
    solver_physical_validation = solver_evidence.get("physical_validation")
    solver_restart = solver_evidence.get("restart_artifact")
    quantitative = solver_evidence.get("quantitative_solver_evidence")

    ids = sorted(set(verification_rows) | set(qualification_rows))
    mismatches: list[str] = []
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
