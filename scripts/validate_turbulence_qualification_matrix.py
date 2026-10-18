"""Validate the canonical 13-model turbulence qualification matrix."""

from __future__ import annotations

import json
import sys

EXPECTED = {
    "LAMINAR",
    "KEPSILON",
    "RNG_KEPSILON",
    "REALIZABLE_KEPSILON",
    "KOMEGA",
    "SST",
    "SPALART_ALLMARAS",
    "SMAGORINSKY",
    "WALE",
    "DYNAMIC_KEQN",
    "DES",
    "DDES",
    "IDDES",
}
STAGES = {
    "equations",
    "constants",
    "implementation_audit",
    "term_tests",
    "invariants",
    "analytical_mms",
    "solver_verification",
    "physical_validation",
    "mesh_yplus",
    "conservation_convergence",
    "documentation",
}


def main(path: str) -> int:
    with open(path, encoding="utf-8") as stream:
        data = json.load(stream)
    models = data.get("models", [])
    keys = {model.get("key") for model in models}
    if keys != EXPECTED or len(models) != 13:
        raise SystemExit(
            f"expected 13 catalogue models, got {len(models)}: {sorted(keys)}"
        )

    for model in models:
        if model.get("family") not in {"LAMINAR", "RANS", "LES", "HYBRID"}:
            raise SystemExit(f"invalid family for {model.get('key')}")
        if model.get("status") not in {
            "SOLVER_READY",
            "COMPONENT_READY",
            "KERNEL_ONLY",
            "PLANNED",
        }:
            raise SystemExit(f"invalid status for {model.get('key')}")
        stages = model.get("stages", {})
        if set(stages) != STAGES or not all(
            isinstance(value, bool) for value in stages.values()
        ):
            raise SystemExit(f"invalid stage map for {model.get('key')}")
        if model["status"] == "PLANNED" and stages["solver_verification"]:
            raise SystemExit(
                f"planned model cannot have solver verification: {model['key']}"
            )
        if model["status"] == "COMPONENT_READY" and any(
            stages[key]
            for key in (
                "analytical_mms",
                "solver_verification",
                "conservation_convergence",
            )
        ):
            raise SystemExit(
                "component-ready model cannot claim transport MMS, production "
                f"solver verification, or coupled conservation: {model['key']}"
            )
        if model["status"] != "SOLVER_READY" and stages["physical_validation"]:
            raise SystemExit(
                f"non-solver-ready model cannot have physical validation: {model['key']}"
            )
        if (
            model["family"] == "RANS"
            and model["key"]
            in {"KEPSILON", "RNG_KEPSILON", "REALIZABLE_KEPSILON", "KOMEGA"}
            and model["wall_distance"]
        ):
            raise SystemExit(f"unexpected wall-distance dependency: {model['key']}")
        if (
            model["key"] in {"SST", "SPALART_ALLMARAS", "DES", "DDES", "IDDES"}
            and not model["wall_distance"]
        ):
            raise SystemExit(f"missing wall-distance dependency: {model['key']}")

    rules = data.get("promotion_rules", {})
    for key in (
        "implementation_is_not_verification",
        "verification_is_not_validation",
        "wall_distance_is_not_turbulence_validation",
        "physical_validation_requires_reference_and_real_solver_execution",
        "no_tolerance_relaxation",
        "no_disabled_validation",
        "production_momentum_coupling_required",
    ):
        if rules.get(key) is not True:
            raise SystemExit(f"promotion rule missing or false: {key}")

    print(f"TURBULENCE_QUALIFICATION_MATRIX: PASS ({len(models)} catalogue models)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
