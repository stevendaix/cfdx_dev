"""Fail-closed 3-grid NACA verification, with independent mesh-family evidence.

Never report a single-ratio Richardson/GCI from N64/N128/N256 solely from
mesh names. The original wall-normal spacing refines roughly 3.2:1 rather
than 2:1, and an apparently monotone Cd cannot repair that inconsistency.
Requires a successful 3-grid CFD solver campaign AND a geometrically similar
mesh-quality family. Distinguishes Cd, Cd_pressure, Cd_viscous, Cl and Cm.
This is discretization verification, not NASA experimental validation.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from naca_mesh_similarity_gate import verify_mesh_family


def estimate_quantity(values: list[float], ratio: float) -> dict:
    if any(not math.isfinite(v) for v in values):
        raise ValueError("nonfinite aerodynamic coefficient")
    coarse, medium, fine = values
    d1, d2 = medium - coarse, fine - medium
    if d1 * d2 <= 0.0 or abs(d2) >= abs(d1):
        return {"status": "NOT_QUALIFIED_NONMONOTONE_OR_PLATEAU", "gci": None}
    order = math.log(abs(d1 / d2)) / math.log(ratio)
    if not math.isfinite(order) or not 0 < order < 8:
        return {"status": "NOT_QUALIFIED_ORDER", "gci": None}
    denom = ratio ** order - 1.0
    extrapolated = fine + (fine - medium) / denom
    absolute_gci = 1.25 * abs(fine - medium) / denom
    return {"status": "RICHARDSON_CANDIDATE", "order": order,
            "extrapolated": extrapolated, "gci_absolute": absolute_gci,
            "gci_relative": absolute_gci / abs(fine) if abs(fine) > 1e-14 else None}


def verify_three_grid(summary: dict, mesh_quality: list[dict] | None = None) -> dict:
    if summary.get("three_grid_solver_gate") != "PASS":
        raise ValueError("the three-grid CFD solver gate did not pass")
    quality = mesh_quality if mesh_quality is not None else summary.get("mesh_quality")
    if not isinstance(quality, list):
        raise ValueError("GCI requires independent mesh-quality records")
    mesh_report = verify_mesh_family(quality)
    if mesh_report["status"] != "COMPATIBLE_WITH_SINGLE_RATIO_GCI":
        raise ValueError("mesh family is not similar: " + ", ".join(mesh_report["reasons"]))
    rows = summary.get("results")
    if not isinstance(rows, list) or len(rows) != 3:
        raise ValueError("exactly three grid results are required")
    expected = (64, 128, 256)
    # Mesh records must correspond one-to-one to the CFD solution records.
    if sorted(int(item["n_surface"]) for item in quality) != list(expected):
        raise ValueError("mesh-quality N labels differ from CFD solution labels")
    names = ("Cd", "Cd_pressure", "Cd_viscous", "Cl", "Cm")
    quantities = {name: [] for name in names}
    for level, row in zip(expected, rows):
        if row.get("N") != level or row.get("status") != "PASS":
            raise ValueError(f"N{level}: missing physical PASS")
        metrics = row.get("metrics")
        if not isinstance(metrics, dict):
            raise ValueError(f"N{level}: metrics must be an object")
        for name in (*names, "continuity"):
            value = metrics.get(name)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
                raise ValueError(f"N{level}: invalid {name}")
        if not (metrics["continuity"] < 1e-7 and abs(metrics["Cl"]) <= 2e-2):
            raise ValueError(f"N{level}: physical conservation or symmetry failed")
        if abs(metrics["Cd"] - metrics["Cd_pressure"] - metrics["Cd_viscous"]) > 1e-8:
            raise ValueError(f"N{level}: pressure and viscous drag do not sum to Cd")
        for name in names:
            quantities[name].append(float(metrics[name]))
    r1 = mesh_report["pairs"][0]["tangential_ratio"]
    r2 = mesh_report["pairs"][1]["tangential_ratio"]
    if abs(r1-r2) / max(r1,r2) > 0.05 or r1 <= 1:
        raise ValueError("unequal nominal refinement ratios")
    ratio = math.sqrt(r1*r2)
    results = {name: estimate_quantity(values, ratio) for name, values in quantities.items()}
    if results["Cd"]["status"] != "RICHARDSON_CANDIDATE":
        raise ValueError("Cd lacks an asymptotic Richardson signature")
    return {"grids": list(expected), "mesh_family": mesh_report,
            "refinement_ratio": ratio, "coefficients": quantities,
            "quantities": results,
            "observed_order": results["Cd"]["order"],
            "richardson_extrapolated_cd": results["Cd"]["extrapolated"],
            "fine_grid_relative_gci": results["Cd"]["gci_relative"],
            "status": "RICHARDSON_CANDIDATE_NOT_EXPERIMENTAL_VALIDATION",
            "reference_validation": "NOT_PERFORMED",
            "wall_shear_validation": "NOT_PERFORMED"}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--mesh-quality", type=Path, required=True,
                        help="independently measured N64/N128/N256 mesh-quality JSON")
    parser.add_argument("--output", type=Path, default=Path("naca_gci.json"))
    args = parser.parse_args()
    try:
        summary = json.loads(args.summary.read_text(encoding="utf-8"))
        quality = json.loads(args.mesh_quality.read_text(encoding="utf-8"))
        report = verify_three_grid(summary, quality)
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as exc:
        args.output.unlink(missing_ok=True)
        print(f"NACA GCI FAIL: {exc}")
        return 2
    temporary = args.output.with_name(args.output.name + ".tmp")
    temporary.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    temporary.replace(args.output)
    print(json.dumps(report, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
