"""Fail-closed NACA three-grid drag verification; requires actual converged data.

The grid levels have an approximate 2:1 linear refinement ratio. An observed
order and GCI are reported only when the data support Richardson extrapolation.
This is numerical verification, NOT validation against an experimental Cd.
"""

import argparse
import json
import math
from pathlib import Path


def verify_three_grid(summary: dict) -> dict:
    if summary.get("three_grid_solver_gate") != "PASS":
        raise ValueError("the three-grid CFD solver gate did not pass")
    rows = summary.get("results")
    if not isinstance(rows, list) or len(rows) != 3:
        raise ValueError("exactly three grid results are required")
    cds = []
    for n, row in zip((64, 128, 256), rows):
        if row.get("N") != n or row.get("status") != "PASS":
            raise ValueError(f"N{n}: missing physical PASS")
        metrics = row.get("metrics")
        if not isinstance(metrics, dict):
            raise TypeError(f"N{n}: metrics must be a mapping")
        for key in ("Cd", "Cl", "continuity", "Cd_pressure", "Cd_viscous"):
            value = metrics.get(key)
            if (
                isinstance(value, bool)
                or not isinstance(value, (int, float))
                or not math.isfinite(value)
            ):
                raise ValueError(f"N{n}: invalid {key}")
        if not (
            0 <= metrics["continuity"] < 1e-7
            and abs(metrics["Cl"]) <= 2e-2
            and metrics["Cd_pressure"] > 0
        ):
            raise ValueError(f"N{n}: physical gate failed")
        if not math.isclose(
            metrics["Cd"],
            metrics["Cd_pressure"] + metrics["Cd_viscous"],
            rel_tol=1e-7,
            abs_tol=1e-9,
        ):
            raise ValueError(f"N{n}: nonconservative drag decomposition")
        cds.append(float(metrics["Cd"]))
    d01 = cds[1] - cds[0]
    d12 = cds[2] - cds[1]
    if not (d01 * d12 > 0 and abs(d12) < abs(d01)):
        raise ValueError("drag has no monotone asymptotic refinement signature")
    order = math.log(abs(d01 / d12)) / math.log(2.0)
    if not (0 < order < 10.0 and math.isfinite(order)):
        raise ValueError("invalid observed Richardson order")
    fine = cds[2]
    if abs(fine) <= 1e-14:
        raise ValueError("zero drag prevents relative GCI")
    extrapolated = fine + (fine - cds[1]) / (2.0**order - 1.0)
    gci_value = 1.25 * abs((fine - cds[1]) / fine) / (2.0**order - 1.0)
    if not (math.isfinite(extrapolated) and math.isfinite(gci_value)):
        raise ValueError("non-finite refinement estimate")
    return {
        # This lightweight summary does NOT inspect mesh-family geometry,
        # physical histories or independent wall-traction certificates.
        # A numerical Richardson estimate is never a qualified force GCI.
        "status": "DIAGNOSTIC_ONLY",
        "grids": [64, 128, 256],
        "cd": cds,
        "observed_order": order,
        "richardson_extrapolated_cd": extrapolated,
        "fine_grid_relative_gci": gci_value,
        "reference_validation": "NOT_PERFORMED",
        "wall_shear_validation": "NOT_PERFORMED",
        "qualification_route": "scripts/naca_force_mesh_qualification.py",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--output", type=Path, default=Path("naca_gci.json"))
    args = parser.parse_args()
    try:
        report = verify_three_grid(json.loads(args.summary.read_text()))
    except (OSError, ValueError, TypeError, KeyError) as exc:
        print(f"NACA GCI FAIL: {exc}")
        return 1
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
