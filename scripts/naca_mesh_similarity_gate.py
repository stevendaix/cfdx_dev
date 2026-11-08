"""Fail-closed mesh-family gate before one-parameter Richardson/GCI.

The source mesh-quality JSON uses n_surface, n_radial and first-layer
spacing. Refining both indices by two does *not* mean first-layer spacing
is refined by two: a strong wall-normal grading can invalidate the assumed
single characteristic h. A failure does not imply that a multi-parameter
mesh-convergence analysis is impossible, only that simple 1-ratio GCI must
not be promoted as physically qualified.
"""

import argparse
import itertools
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
    keys = (
        "n_surface",
        "n_radial",
        "radius",
        "radial_power",
        "minimum_first_layer_spacing",
        "mean_first_layer_spacing",
        "maximum_first_layer_spacing",
        "minimum_absolute_cell_area",
        "maximum_cell_edge_aspect_ratio",
        "minimum_radial_growth_ratio",
        "maximum_radial_growth_ratio",
    )
    # A valid one-parameter family must have comparable *distributions*, not
    # merely the same N and minimum near-wall spacing. Missing quality data
    # cannot be interpreted as evidence of similarity.
    if [row.get("n_surface") for row in items] != [64, 128, 256]:
        diagnostic["reasons"].append("wrong_or_duplicate_naca_mesh_levels")
    for row in items:
        for integer_key in ("n_surface", "n_radial", "inverted_or_degenerate_cells"):
            if type(row.get(integer_key)) is not int:
                raise ValueError(f"invalid integer mesh characteristic {integer_key}")
        for key in keys:
            value = float(row[key])
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"invalid mesh characteristic {key}")
        if int(row["inverted_or_degenerate_cells"]) != 0:
            diagnostic["reasons"].append("inverted_or_degenerate_cells")
        max_angle = float(row["maximum_nonorthogonality_degrees"])
        p95 = float(row["p95_nonorthogonality_degrees"])
        if not (
            math.isfinite(max_angle)
            and math.isfinite(p95)
            and 0 <= p95 <= max_angle < 90
        ):
            raise ValueError("invalid nonorthogonality distribution")
        if max_angle >= 80:
            diagnostic["reasons"].append("near_degenerate_nonorthogonality")
        te = float(row.get("trailing_edge_regularization", 0.0))
        if not math.isfinite(te) or not 0.0 <= te <= 0.08:
            raise ValueError("invalid trailing-edge regularization")
    strategies = {row.get("radial_refinement_strategy", "unspecified") for row in items}
    if len(strategies) != 1 or "unspecified" in strategies:
        diagnostic["reasons"].append(
            "inconsistent_or_missing_radial_refinement_strategy"
        )
    for a, b in itertools.pairwise(items):
        tangential = float(b["n_surface"]) / float(a["n_surface"])
        radial = float(b["n_radial"]) / float(a["n_radial"])
        wall = float(a["minimum_first_layer_spacing"]) / float(
            b["minimum_first_layer_spacing"]
        )
        mismatch = abs(wall / tangential - 1.0)
        means = float(a["mean_first_layer_spacing"]) / float(
            b["mean_first_layer_spacing"]
        )
        maxima = float(a["maximum_first_layer_spacing"]) / float(
            b["maximum_first_layer_spacing"]
        )
        max_angle_change = abs(
            float(b["maximum_nonorthogonality_degrees"])
            - float(a["maximum_nonorthogonality_degrees"])
        )
        p95_angle_change = abs(
            float(b["p95_nonorthogonality_degrees"])
            - float(a["p95_nonorthogonality_degrees"])
        )
        pair = {
            "from": int(a["n_surface"]),
            "to": int(b["n_surface"]),
            "tangential_ratio": tangential,
            "radial_ratio": radial,
            "first_layer_ratio": wall,
            "relative_wall_ratio_mismatch": mismatch,
            "mean_first_layer_ratio": means,
            "maximum_first_layer_ratio": maxima,
            "max_nonorth_angle_change_deg": max_angle_change,
            "p95_nonorth_angle_change_deg": p95_angle_change,
            "aspect_ratio_coarse": float(a["maximum_cell_edge_aspect_ratio"]),
            "aspect_ratio_fine": float(b["maximum_cell_edge_aspect_ratio"]),
        }
        diagnostic["pairs"].append(pair)
        if tangential <= 1 or abs(radial / tangential - 1.0) > tolerance:
            diagnostic["reasons"].append("different_tangential_and_radial_ratios")
        if (
            max(mismatch, abs(means / tangential - 1.0), abs(maxima / tangential - 1.0))
            > tolerance
        ):
            diagnostic["reasons"].append("first_layer_distribution_not_similar")
        if max_angle_change > 10.0 or p95_angle_change > 10.0:
            diagnostic["reasons"].append(
                "nonorthogonality_distribution_changes_with_refinement"
            )
        if (
            abs(float(b["radius"]) / float(a["radius"]) - 1) > 1e-12
            or abs(float(b["radial_power"]) - float(a["radial_power"])) > 1e-12
        ):
            diagnostic["reasons"].append("different_geometry_or_grading")
        if (
            abs(
                float(b.get("trailing_edge_regularization", 0.0))
                - float(a.get("trailing_edge_regularization", 0.0))
            )
            > 1e-12
        ):
            diagnostic["reasons"].append("different_trailing_edge_sampling")
    if not diagnostic["reasons"]:
        diagnostic["status"] = "COMPATIBLE_WITH_SINGLE_RATIO_GCI"
    else:
        diagnostic["reasons"] = sorted(set(diagnostic["reasons"]))
    return diagnostic


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quality_json", type=Path)
    parser.add_argument("--tolerance", type=float, default=0.20)
    parser.add_argument(
        "--report-only",
        action="store_true",
        help="print diagnostic but do not certify single-ratio GCI",
    )
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
