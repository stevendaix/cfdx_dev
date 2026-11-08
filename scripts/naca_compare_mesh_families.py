"""Compare two NACA mesh families without confusing anisotropy with convergence.

A single-grid drag agreement is never grid convergence. A qualified
comparison needs two independent three-grid physical certificates, each
produced by naca_force_mesh_qualification.py, on identical physical setups.
Geometry-only analysis is available to diagnose incompatible wall resolution.
The ~0.119 reference is a published *numerical* NACA0012/Re=1000/AoA=0
comparison, not an experimental measurement or a calibration target.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path

REQUIRED_PARAMETERS = (
    "airfoil",
    "reynolds",
    "angle_deg",
    "chord",
    "farfield_radius_chords",
    "velocity_boundary",
    "pressure_boundary",
    "convection_scheme",
    "wall_shear_method",
    "solver_revision",
)
COMPONENTS = ("cd", "cdp", "cdv")


def _finite_positive(value: object, name: str) -> float:
    number = float(value)
    if not math.isfinite(number) or number <= 0:
        raise ValueError(f"{name}: expected finite, positive number")
    return number


def _quality_records(value: object) -> dict[int, dict]:
    if not isinstance(value, list) or not value:
        raise TypeError("mesh quality must be a nonempty JSON list")
    records = {}
    for row in value:
        if not isinstance(row, dict):
            raise TypeError("invalid mesh quality record")
        level = int(row["n_surface"])
        if level in records:
            raise ValueError(f"duplicate mesh N{level}")
        for key in (
            "n_radial",
            "radius",
            "minimum_first_layer_spacing",
            "mean_first_layer_spacing",
            "maximum_cell_edge_aspect_ratio",
        ):
            _finite_positive(row[key], key)
        if int(row["inverted_or_degenerate_cells"]) != 0:
            raise ValueError(f"N{level}: degenerate/inverted cells")
        records[level] = row
    return records


def _mesh_diagnostics(a: dict[int, dict], b: dict[int, dict]) -> dict:
    levels = sorted(set(a) & set(b))
    if not levels:
        raise ValueError("no common mesh refinement level")
    pairs = []
    reasons = []
    for level in levels:
        ar, br = a[level], b[level]
        if (
            int(ar["n_radial"]) != int(br["n_radial"])
            or abs(float(ar["radius"]) - float(br["radius"])) > 1e-12
        ):
            reasons.append(f"N{level}: different radial count or farfield extent")
        h_a = float(ar["minimum_first_layer_spacing"])
        h_b = float(br["minimum_first_layer_spacing"])
        layer_ratio = max(h_a, h_b) / min(h_a, h_b)
        if layer_ratio > 1.25:
            reasons.append(f"N{level}: first-wall-layer mismatch >25%")
        pairs.append(
            {
                "N": level,
                "first_wall_layer_a": h_a,
                "first_wall_layer_b": h_b,
                "first_wall_layer_ratio": layer_ratio,
                "max_edge_aspect_ratio_a": float(ar["maximum_cell_edge_aspect_ratio"]),
                "max_edge_aspect_ratio_b": float(br["maximum_cell_edge_aspect_ratio"]),
            }
        )
    return {
        "status": (
            "COMPARABLE_WALL_RESOLUTION" if not reasons else "DIFFERENT_RESOLUTION"
        ),
        "pairs": pairs,
        "reasons": reasons,
        "note": (
            "Changing radial grading is not equivalent to a controlled "
            "isotropic/anisotropic comparison at equal near-wall resolution."
        ),
    }


def _candidate(result: object, family: str) -> dict:
    if not isinstance(result, dict):
        raise TypeError(f"{family}: invalid physical report")
    if result.get("status") != "GCI_CANDIDATE":
        raise ValueError(
            f"{family}: complete three-grid nonlinear certification missing"
        )
    if (
        result.get("mesh_family", {}).get("status")
        != "COMPATIBLE_WITH_SINGLE_RATIO_GCI"
    ):
        raise ValueError(f"{family}: three-grid mesh similarity not certified")
    gci = result.get("gci", {})
    if gci.get("status") != "RICHARDSON_CANDIDATE_ONLY":
        raise ValueError(f"{family}: independent Richardson evidence missing")
    cases = result.get("cases")
    if (
        not isinstance(cases, list)
        or len(cases) != 3
        or any(case.get("status") != "CONVERGED" for case in cases)
    ):
        raise ValueError(f"{family}: not three physically converged solutions")
    return gci["result"]["quantities"]


def compare(
    quality_a: object,
    quality_b: object,
    report_a: object | None = None,
    report_b: object | None = None,
    metadata_a: object | None = None,
    metadata_b: object | None = None,
    reference_cd: float = 0.119,
    relative_tolerance: float = 0.03,
    max_relative_gci: float = 0.10,
) -> dict:
    reference_cd = _finite_positive(reference_cd, "reference_cd")
    if not math.isfinite(relative_tolerance) or relative_tolerance <= 0:
        raise ValueError("relative_tolerance must be positive")
    if not math.isfinite(max_relative_gci) or max_relative_gci <= 0:
        raise ValueError("max_relative_gci must be positive")
    a, b = _quality_records(quality_a), _quality_records(quality_b)
    geometry = _mesh_diagnostics(a, b)
    result = {
        "status": "BLOCKED",
        "geometry": geometry,
        "gci_comparison": {},
        "reference": {
            "cd": reference_cd,
            "type": "PUBLISHED_NUMERICAL_BENCHMARK_NOT_EXPERIMENTAL",
            "physical_validation": "NOT_PERFORMED",
        },
        "reasons": [],
    }
    if report_a is None or report_b is None:
        result["reasons"].append(
            "two certified three-grid physical GCI reports required"
        )
        return result
    if metadata_a is None or metadata_b is None:
        result["reasons"].append(
            "two solver/physics provenance metadata records required"
        )
        return result
    if not isinstance(metadata_a, dict) or not isinstance(metadata_b, dict):
        raise TypeError("physics metadata must be JSON objects")
    for key in REQUIRED_PARAMETERS:
        if key not in metadata_a or key not in metadata_b:
            result["reasons"].append(f"missing physical provenance: {key}")
        elif metadata_a[key] != metadata_b[key]:
            result["reasons"].append(f"incompatible physical setup: {key}")
    if result["reasons"]:
        return result
    for family, records in (("A", a), ("B", b)):
        if set(records) != {64, 128, 256}:
            result["reasons"].append(
                f"{family}: exactly N64/N128/N256 mesh-quality evidence required"
            )
    if result["reasons"]:
        return result
    try:
        values_a, values_b = _candidate(report_a, "A"), _candidate(report_b, "B")
        all_agree = True
        for name in COMPONENTS:
            aa = values_a[name]
            bb = values_b[name]
            if (
                aa["status"] != "asymptotic_candidate"
                or bb["status"] != "asymptotic_candidate"
            ):
                raise ValueError(f"{name}: non-asymptotic force component")
            xa = float(aa["extrapolated"])
            xb = float(bb["extrapolated"])
            ga = float(aa["GCI_absolute"])
            gb = float(bb["GCI_absolute"])
            if not all(math.isfinite(v) for v in (xa, xb, ga, gb)) or ga < 0 or gb < 0:
                raise ValueError(
                    f"{name}: nonfinite/negative extrapolation or uncertainty"
                )
            scale = max(abs(xa), abs(xb), 1e-12)
            discrepancy = abs(xa - xb)
            relative = discrepancy / scale
            bands_overlap = discrepancy <= ga + gb
            precise = max(ga, gb) / scale <= max_relative_gci
            agrees = bands_overlap and precise and relative <= relative_tolerance
            all_agree = all_agree and agrees
            result["gci_comparison"][name] = {
                "extrapolated_a": xa,
                "extrapolated_b": xb,
                "absolute_difference": discrepancy,
                "relative_difference": relative,
                "gci_absolute_a": ga,
                "gci_absolute_b": gb,
                "gci_bands_overlap": bands_overlap,
                "uncertainty_bounded": precise,
                "within_requested_tolerance": relative <= relative_tolerance,
                "agreement_candidate": agrees,
            }
    except (ValueError, TypeError, KeyError, OverflowError) as exc:
        result["reasons"].append(str(exc))
        return result
    result["reference"]["relative_error_a"] = (
        result["gci_comparison"]["cd"]["extrapolated_a"] - reference_cd
    ) / reference_cd
    result["reference"]["relative_error_b"] = (
        result["gci_comparison"]["cd"]["extrapolated_b"] - reference_cd
    ) / reference_cd
    if all_agree:
        result["status"] = "CROSS_MESH_GCI_AGREEMENT_CANDIDATE"
    else:
        result["reasons"].append(
            "Cd, pressure and viscous drag must each agree within "
            "bounded GCI uncertainty"
        )
    # Even a successful cross-mesh comparison does NOT verify Cp(x), Cf(x),
    # body wake, farfield sensitivity or agreement to experimental data.
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--a-quality", type=Path, required=True)
    parser.add_argument("--b-quality", type=Path, required=True)
    parser.add_argument("--a-gci", type=Path)
    parser.add_argument("--b-gci", type=Path)
    parser.add_argument("--a-metadata", type=Path)
    parser.add_argument("--b-metadata", type=Path)
    parser.add_argument("--reference-cd", type=float, default=0.119)
    parser.add_argument("--relative-tolerance", type=float, default=0.03)
    parser.add_argument("--max-relative-gci", type=float, default=0.10)
    parser.add_argument("--report-only", action="store_true")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    try:

        def load(path: Path | None) -> object | None:
            return json.loads(path.read_text(encoding="utf-8")) if path else None

        report = compare(
            load(args.a_quality),
            load(args.b_quality),
            load(args.a_gci),
            load(args.b_gci),
            load(args.a_metadata),
            load(args.b_metadata),
            args.reference_cd,
            args.relative_tolerance,
            args.max_relative_gci,
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        tmp = args.output.with_name(args.output.name + ".tmp")
        tmp.write_text(
            json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        os.replace(tmp, args.output)
        print(json.dumps(report, indent=2, sort_keys=True))
        if report["status"] == "CROSS_MESH_GCI_AGREEMENT_CANDIDATE":
            return 0
        return 0 if args.report_only else 2
    except (OSError, ValueError, TypeError, KeyError, OverflowError) as exc:
        print(f"NACA CROSS-MESH QUALIFICATION FAIL: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())