"""Fail-closed mesh-family gate before one-parameter Richardson/GCI.

The source mesh-quality JSON uses n_surface, n_radial and first-layer
spacing. Refining both indices by two does *not* mean first-layer spacing
is refined by two: a strong wall-normal grading can invalidate the assumed
single characteristic h. A failure does not imply that a multi-parameter
mesh-convergence analysis is impossible, only that simple 1-ratio GCI must
not be promoted as physically qualified.
"""
import argparse
import json
import math
import sys
from pathlib import Path


def verify_mesh_family(rows: list[dict], tolerance: float = 0.20) -> dict:
    if not math.isfinite(tolerance) or not 0.0 < tolerance < 1.0:
        raise ValueError("invalid similarity tolerance")
    if len(rows) != 3:
        raise ValueError("exactly three mesh-quality records required")
    items = sorted(rows, key=lambda row: int(row["n_surface"]))
    diagnostic = {"status": "BLOCKED", "pairs": [], "reasons": []}
    keys = ("n_surface", "n_radial", "radius", "radial_power",
            "minimum_first_layer_spacing")
    for row in items:
        for key in keys:
            value = float(row[key])
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"invalid mesh characteristic {key}")
        if int(row["inverted_or_degenerate_cells"]) != 0:
            diagnostic["reasons"].append("inverted_or_degenerate_cells")
    strategies = {row.get("radial_refinement_strategy", "unspecified") for row in items}
    if len(strategies) != 1:
        diagnostic["reasons"].append("inconsistent_radial_refinement_strategy")
    for a, b in zip(items, items[1:]):
        tangential = float(b["n_surface"]) / float(a["n_surface"])
        radial = float(b["n_radial"]) / float(a["n_radial"])
                raise ValueError(f"invalid mesh characteristic {key}")
        if int(row["inverted_or_degenerate_cells"]) != 0:
            diagnostic["reasons"].append("inverted_or_degenerate_cells")
    for a, b in zip(items, items[1:]):
        tangential = float(b["n_surface"]) / float(a["n_surface"])
        radial = float(b["n_radial"]) / float(a["n_radial"])
        wall = (float(a["minimum_first_layer_spacing"]) /
                float(b["minimum_first_layer_spacing"]))
        mismatch = abs(wall / tangential - 1.0)
        pair = {"from": int(a["n_surface"]), "to": int(b["n_surface"]),
                "tangential_ratio": tangential, "radial_ratio": radial,
                "first_layer_ratio": wall, "relative_wall_ratio_mismatch": mismatch}
        diagnostic["pairs"].append(pair)
        if tangential <= 1 or abs(radial / tangential - 1.0) > tolerance:
            diagnostic["reasons"].append("different_tangential_and_radial_ratios")
        if mismatch > tolerance:
            diagnostic["reasons"].append("first_layer_refinement_not_similar")
        if (abs(float(b["radius"]) / float(a["radius"]) - 1) > 1e-12 or
                abs(float(b["radial_power"]) - float(a["radial_power"])) > 1e-12):
            diagnostic["reasons"].append("different_geometry_or_grading")
    if not diagnostic["reasons"]:
        diagnostic["status"] = "COMPATIBLE_WITH_SINGLE_RATIO_GCI"
    else:
        diagnostic["reasons"] = sorted(set(diagnostic["reasons"]))
    return diagnostic


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quality_json", type=Path)
    parser.add_argument("--tolerance", type=float, default=0.20)
    parser.add_argument("--report-only", action="store_true",
                        help="print diagnostic but do not certify single-ratio GCI")
    args = parser.parse_args()
    rows = json.loads(args.quality_json.read_text())
    result = verify_mesh_family(rows, args.tolerance)
    print(json.dumps(result, indent=2, sort_keys=True))
    if result["status"] != "COMPATIBLE_WITH_SINGLE_RATIO_GCI":
        if args.report_only:
            return 0
        print("GCI MESH SIMILARITY NOT QUALIFIED", file=sys.stderr)
        return 2
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError) as error:
        print(f"MESH FAMILY QUALIFICATION FAIL: {error}", file=sys.stderr)
        sys.exit(2)
