"""Compare two converged NACA wall-face CSVs without claiming GCI.

This is a spatial attribution of changes in pressure/viscous forces.
A physical reference and a third converged grid are still mandatory.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

BINS = ((0, 0.05), (0.05, 0.10), (0.10, 0.25), (0.25, 0.50),
        (0.50, 0.75), (0.75, 0.90), (0.90, 1.001))
COLUMNS = ("face", "x_over_c", "y_over_c", "Cp", "Cf_tangent",
           "Cd_pressure_face", "Cd_viscous_face")


def read_profile(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        if tuple(reader.fieldnames or ()) != COLUMNS:
            raise ValueError(f"bad NACA wall-profile schema in {path}")
        rows = [{k: float(item[k]) for k in COLUMNS} for item in reader]
    if len(rows) < 8 or any(not all(math.isfinite(v) for v in r.values()) for r in rows):
        raise ValueError("insufficient or non-finite wall-face data")
    if any(not 0 <= r["x_over_c"] <= 1.001 for r in rows):
        raise ValueError("invalid x/c coordinate")
    if len(set(int(r["face"]) for r in rows)) != len(rows):
        raise ValueError("duplicate face id")
    return rows


def audit(coarse: list[dict[str, float]], fine: list[dict[str, float]]) -> dict:
    if len(fine) <= len(coarse):
        raise ValueError("fine grid must have more wall faces")
    regions = []
    for lo, hi in BINS:
        c = [r for r in coarse if lo <= r["x_over_c"] < hi]
        f = [r for r in fine if lo <= r["x_over_c"] < hi]
        if not c or not f:
            raise ValueError("missing pressure/drag data in a chordwise region")
        regions.append({
            "x_over_c": [lo, hi],
            "coarse_faces": len(c), "fine_faces": len(f),
            "cd_pressure_coarse": sum(r["Cd_pressure_face"] for r in c),
            "cd_pressure_fine": sum(r["Cd_pressure_face"] for r in f),
            "cd_viscous_coarse": sum(r["Cd_viscous_face"] for r in c),
            "cd_viscous_fine": sum(r["Cd_viscous_face"] for r in f),
        })
    for region in regions:
        region["delta_cd_pressure"] = (
            region["cd_pressure_fine"] - region["cd_pressure_coarse"])
        region["delta_cd_viscous"] = (
            region["cd_viscous_fine"] - region["cd_viscous_coarse"])
    cdp_coarse = sum(r["Cd_pressure_face"] for r in coarse)
    cdp_fine = sum(r["Cd_pressure_face"] for r in fine)
    cdv_coarse = sum(r["Cd_viscous_face"] for r in coarse)
    cdv_fine = sum(r["Cd_viscous_face"] for r in fine)
    total_delta = (cdp_fine + cdv_fine) - (cdp_coarse + cdv_coarse)
    if not math.isfinite(total_delta):
        raise ValueError("non-finite integrated force")
    return {
        "status": "TWO_GRID_DIAGNOSTIC_NOT_QUALIFIED",
        "physical_reference_validation": "NOT_PERFORMED",
        "gci": "NOT_PERFORMED_NEEDS_CONVERGED_N256",
        "coarse_wall_faces": len(coarse),
        "fine_wall_faces": len(fine),
        "cd_coarse": cdp_coarse + cdv_coarse,
        "cd_fine": cdp_fine + cdv_fine,
        "cd_pressure_delta": cdp_fine - cdp_coarse,
        "cd_viscous_delta": cdv_fine - cdv_coarse,
        "cd_total_delta": total_delta,
        "cp_max_coarse": max(r["Cp"] for r in coarse),
        "cp_max_fine": max(r["Cp"] for r in fine),
        "regions": regions,
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coarse", type=Path, required=True)
    parser.add_argument("--fine", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = audit(read_profile(args.coarse), read_profile(args.fine))
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print("NACA_GRID_PRESSURE_AUDIT: TWO_GRID_DIAGNOSTIC_NOT_QUALIFIED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
