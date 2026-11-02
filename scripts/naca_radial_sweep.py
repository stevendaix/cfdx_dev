#!/usr/bin/env python3
"""Fail-closed wall-normal geometry sweep; this is not CFD qualification."""

import argparse
import json
import math
from pathlib import Path

from generate_2d_validation_meshes import naca_o_grid_quality


def sweep(n, nr, powers, max_nonorth=65.0):
    if not math.isfinite(max_nonorth) or max_nonorth <= 0.0:
        raise ValueError("max_nonorth must be finite and positive")
    baseline = naca_o_grid_quality(n, nr, radial_power=1.8)
    records = []
    for power in powers:
        if not math.isfinite(power) or not 1.0 <= power <= 3.0:
            raise ValueError("invalid radial power")
        quality = naca_o_grid_quality(n, nr, radial_power=power)
        geometry_pass = (
            quality["inverted_or_degenerate_cells"] == 0
            and quality["maximum_nonorthogonality_degrees"] <= max_nonorth
            and quality["minimum_absolute_cell_area"] > 0.0
            and math.isfinite(quality["maximum_radial_growth_ratio"])
            and math.isfinite(quality["minimum_first_layer_spacing"])
        )
        quality["geometry_pass"] = geometry_pass
        quality["first_layer_ratio_to_baseline"] = (
            quality["minimum_first_layer_spacing"]
            / baseline["minimum_first_layer_spacing"]
        )
        records.append(quality)
    return {
        "mesh": f"N{n}",
        "baseline": baseline,
        "candidates": records,
        "warning": (
            "Geometry-only screening; wall y+, interpolation quality, "
            "volume-neighbor ratios and CFD convergence remain unqualified."
        ),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=64)
    parser.add_argument("--nr", type=int, default=24)
    parser.add_argument("--powers", type=float, nargs="+", default=[1.4, 1.6, 1.8, 2.0])
    parser.add_argument("--max-nonorth", type=float, default=65.0)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    result = sweep(args.n, args.nr, args.powers, args.max_nonorth)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    print(
        json.dumps(
            {
                "mesh": result["mesh"],
                "candidates": len(result["candidates"]),
                "passed_geometry": sum(
                    candidate["geometry_pass"] for candidate in result["candidates"]
                ),
            }
        )
    )


if __name__ == "__main__":
    main()
