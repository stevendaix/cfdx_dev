#!/usr/bin/env python3
"""Build the N10 public-fixture numerical applicability matrix.

This is an applicability contract, not a numerical qualification campaign.
It consumes the measured production-import evidence and deliberately refuses
to infer solver/accuracy applicability from provenance or import success alone.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def load(path: Path) -> dict:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except Exception as exc:
        raise SystemExit(f"error: unable to read JSON report {path}: {exc}") from exc
    if not isinstance(value, dict):
        raise SystemExit(f"error: report root is not an object: {path}")
    return value


def fixture_map(report: dict) -> dict[str, dict]:
    fixtures = report.get("fixtures")
    if not isinstance(fixtures, list):
        raise SystemExit("error: production fixture qualification has no fixture list")
    result: dict[str, dict] = {}
    for item in fixtures:
        if not isinstance(item, dict) or not isinstance(item.get("id"), str):
            raise SystemExit("error: malformed fixture entry in qualification report")
        result[item["id"]] = item
    return result


def classify(item: dict) -> tuple[str, str, list[str]]:
    fixture_id = str(item["id"])
    if not all(item.get(key) is True for key in ("imported", "topology_valid", "geometry_quality_valid")):
        return (
            "BLOCKED",
            "Production import/topology/geometry evidence is not PASS; numerical applicability cannot be assessed.",
            [],
        )

    if fixture_id == "meshio-su2-square":
        return (
            "NOT_APPLICABLE_TO_CURRENT_3D_SOLVER",
            "The imported SU2 fixture is a 2-D edge/surface representation; it must not be forced through the current 3-D incompressible production solver.",
            ["2d_formulation", "fixture_specific_2d_accuracy"],
        )

    if fixture_id == "meshio-vtk-unstructured":
        return (
            "APPLICABLE_FOR_3D_EXECUTION_SMOKE",
            "The imported VTK fixture has positive 3-D cell volumes and already has production-solver execution/convergence evidence. Physical accuracy and qualification remain unclaimed.",
            ["fixture_specific_accuracy", "physical_validation"],
        )

    if fixture_id == "meshio-gmsh-insulated-2-2":
        return (
            "NUMERICAL_APPLICABILITY_NOT_CLAIMED",
            "Production import/geometry evidence is PASS, but no fixture-specific PDE, boundary-condition, reference solution, or observed-order contract has been established.",
            ["fixture_specific_pde", "fixture_specific_accuracy", "observed_order"],
        )

    if fixture_id == "openfoam-airfoil2d":
        return (
            "NUMERICAL_APPLICABILITY_NOT_CLAIMED",
            "Production import evidence is PASS, but OpenFOAM airFoil2D semantics require an explicit CFDX 2-D/empty-boundary applicability contract before numerical qualification.",
            ["2d_formulation", "fixture_specific_accuracy", "physical_validation"],
        )

    return (
        "NUMERICAL_APPLICABILITY_NOT_CLAIMED",
        "No fixture-specific numerical applicability contract exists for this public fixture.",
        ["fixture_specific_pde", "fixture_specific_accuracy"],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--qualification", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    args = parser.parse_args()

    qualification = load(args.qualification)
    fixtures = fixture_map(qualification)

    required = (
        "meshio-su2-square",
        "meshio-gmsh-insulated-2-2",
        "meshio-vtk-unstructured",
        "openfoam-airfoil2d",
    )
    missing = [fixture_id for fixture_id in required if fixture_id not in fixtures]
    if missing:
        raise SystemExit("error: qualification report is missing fixtures: " + ", ".join(missing))

    matrix = []
    for fixture_id in required:
        status, rationale, next_steps = classify(fixtures[fixture_id])
        matrix.append(
            {
                "id": fixture_id,
                "import_evidence": {
                    "imported": fixtures[fixture_id].get("imported"),
                    "topology_valid": fixtures[fixture_id].get("topology_valid"),
                    "geometry_quality_valid": fixtures[fixture_id].get("geometry_quality_valid"),
                    "validation_mode": fixtures[fixture_id].get("validation_mode", "3d"),
                },
                "numerical_applicability": status,
                "rationale": rationale,
                "required_next_evidence": next_steps,
            }
        )

    blocked = [item["id"] for item in matrix if item["numerical_applicability"] == "BLOCKED"]
    report = {
        "campaign": "N10 public-fixture numerical applicability matrix",
        "status": "PASS" if not blocked else "FAIL",
        "matrix": matrix,
        "policy": {
            "import_success_is_not_numerical_qualification": True,
            "physical_validation_claimed": False,
            "forces_unsupported_2d_fixtures_through_3d_solver": False,
            "changes_numerical_tolerances": False,
            "disables_validation": False,
            "silent_fallbacks": False,
        },
    }

    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return 0 if not blocked else 1


if __name__ == "__main__":
    raise SystemExit(main())
