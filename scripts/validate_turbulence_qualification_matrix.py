#!/usr/bin/env python3
"""Validate the canonical 13-model turbulence qualification matrix."""
from __future__ import annotations
import json
import sys

EXPECTED = {
    "LAMINAR","KEPSILON","RNG_KEPSILON","REALIZABLE_KEPSILON","KOMEGA","SST",
    "SPALART_ALLMARAS","SMAGORINSKY","WALE","DYNAMIC_KEQN","DES","DDES","IDDES",
}
STAGES = {
    "equations","constants","implementation_audit","term_tests","invariants",
    "analytical_mms","solver_verification","physical_validation","mesh_yplus",
    "conservation_convergence","documentation",
}

def main(path: str) -> int:
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    models = data.get("models", [])
    keys = {m.get("key") for m in models}
    if keys != EXPECTED or len(models) != 13:
        raise SystemExit(f"expected 13 catalogue models, got {len(models)}: {sorted(keys)}")
    for m in models:
        if m.get("family") not in {"LAMINAR","RANS","LES","HYBRID"}:
            raise SystemExit(f"invalid family for {m.get('key')}")
        if m.get("status") not in {"SOLVER_READY","KERNEL_ONLY","PLANNED"}:
            raise SystemExit(f"invalid status for {m.get('key')}")
        stages = m.get("stages", {})
        if set(stages) != STAGES or not all(isinstance(v, bool) for v in stages.values()):
            raise SystemExit(f"invalid stage map for {m.get('key')}")
        if m["status"] == "PLANNED" and stages["solver_verification"]:
            raise SystemExit(f"planned model cannot have solver verification: {m['key']}")
        if m["status"] != "SOLVER_READY" and stages["physical_validation"]:
            raise SystemExit(f"non-solver-ready model cannot have physical validation: {m['key']}")
        if m["family"] == "RANS" and m["key"] in {"KEPSILON","RNG_KEPSILON","REALIZABLE_KEPSILON","KOMEGA"} and m["wall_distance"]:
            raise SystemExit(f"unexpected wall-distance dependency: {m['key']}")
        if m["key"] in {"SST","SPALART_ALLMARAS","DES","DDES","IDDES"} and not m["wall_distance"]:
            raise SystemExit(f"missing wall-distance dependency: {m['key']}")
    rules = data.get("promotion_rules", {})
    for key in ("implementation_is_not_verification","verification_is_not_validation",
                "wall_distance_is_not_turbulence_validation","physical_validation_requires_reference_and_real_solver_execution",
                "no_tolerance_relaxation","no_disabled_validation"):
        if rules.get(key) is not True:
            raise SystemExit(f"promotion rule missing or false: {key}")
    print(f"TURBULENCE_QUALIFICATION_MATRIX: PASS ({len(models)} catalogue models)")
    return 0

if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
