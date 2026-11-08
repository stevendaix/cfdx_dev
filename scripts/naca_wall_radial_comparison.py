"""Compare wall forces on two NACA grids with identical surface sampling.

This explicitly detects cancellation of the pressure/viscous drag errors.
A small total Cd difference is NEVER a mesh-convergence certificate and
this tool does not replace a 3-grid GCI or experimental validation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from naca_wall_pressure_grid_audit import BINS, read_profile


def read_force(path: Path) -> dict[str, float]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        rows = list(reader)
        if not reader.fieldnames or len(rows) != 1:
            raise ValueError("expected one force certificate")
    record = rows[0]
    if record.get("converged", "").lower() not in {"yes", "true", "1"}:
        raise ValueError("force certificate is not converged")
    out = {key: float(record[key]) for key in ("cd", "cdp", "cdv")}
    if not all(math.isfinite(v) for v in out.values()):
        raise ValueError("non-finite force certificate")
    if not math.isclose(out["cd"], out["cdp"] + out["cdv"], abs_tol=1e-9):
        raise ValueError("inconsistent total drag certificate")
    return out


def compare(
    a: list[dict[str, float]], b: list[dict[str, float]],
    certificate_a: dict[str, float], certificate_b: dict[str, float],
    coordinate_tolerance: float = 1e-9,
) -> dict:
    if not 0 < coordinate_tolerance < 1e-3:
        raise ValueError("invalid surface-coordinate tolerance")
    if len(a) != len(b) or len(a) < 32:
        raise ValueError("identical wall-surface face counts required")
    def sort_surface(rows: list[dict[str, float]]):
        return sorted(rows, key=lambda r: (
            0 if r["y_over_c"] >= 0 else 1, r["x_over_c"], r["y_over_c"]
        ))
    a, b = sort_surface(a), sort_surface(b)
    for ai, bi in zip(a, b):
        for coord in ("x_over_c", "y_over_c"):
            if abs(ai[coord] - bi[coord]) > coordinate_tolerance:
                raise ValueError("wall surface sampling is not geometrically comparable")
    def check_integrals(rows: list[dict[str, float]], cert: dict[str, float]):
        for comp, col in (("cdp", "Cd_pressure_face"),
                          ("cdv", "Cd_viscous_face")):
            integrated = sum(row[col] for row in rows)
            if not math.isclose(integrated, cert[comp], rel_tol=2e-7, abs_tol=2e-9):
                raise ValueError(f"force certificate does not match face sum: {comp}")
    check_integrals(a, certificate_a)
    check_integrals(b, certificate_b)
    regions = []
    for lo, hi in BINS:
        paired = [(ai, bi) for ai, bi in zip(a, b) if lo <= ai["x_over_c"] < hi]
        if not paired:
            raise ValueError("one chordwise region lacks face samples")
        regions.append({
            "x_over_c": [lo, hi],
            "face_count": len(paired),
            "delta_cd_pressure": sum(bi["Cd_pressure_face"] - ai["Cd_pressure_face"]
                                     for ai, bi in paired),
            "delta_cd_viscous": sum(bi["Cd_viscous_face"] - ai["Cd_viscous_face"]
                                    for ai, bi in paired),
            "max_abs_delta_cp": max(abs(ai["Cp"] - bi["Cp"]) for ai, bi in paired),
            "max_abs_delta_cf": max(abs(ai["Cf_tangent"] - bi["Cf_tangent"])
                                    for ai, bi in paired),
        })
    p_delta = certificate_b["cdp"] - certificate_a["cdp"]
    v_delta = certificate_b["cdv"] - certificate_a["cdv"]
    t_delta = certificate_b["cd"] - certificate_a["cd"]
    scale = max(abs(certificate_a["cd"]), abs(certificate_b["cd"]), 1e-15)
    p_scale = max(abs(certificate_a["cdp"]), abs(certificate_b["cdp"]), 1e-15)
    v_scale = max(abs(certificate_a["cdv"]), abs(certificate_b["cdv"]), 1e-15)
    if abs((p_delta + v_delta) - t_delta) > 1e-9:
        raise ValueError("pressure+viscous deltas do not reproduce total drag")
    cancellation = abs(t_delta) / max(abs(p_delta) + abs(v_delta), 1e-15)
    deceptive = (abs(t_delta) / scale <= 0.03 and
                 max(abs(p_delta) / p_scale, abs(v_delta) / v_scale) >= 0.05)
    return {
        "status": "TWO_GRID_DIAGNOSTIC_NOT_QUALIFIED",
        "false_total_cd_agreement_risk": deceptive,
        "force_cancellation_fraction": 1.0 - cancellation,
        "cd_a": certificate_a["cd"], "cd_b": certificate_b["cd"],
        "relative_total_cd_difference": abs(t_delta)/scale,
        "cd_pressure_a": certificate_a["cdp"],
        "cd_pressure_b": certificate_b["cdp"],
        "relative_pressure_difference": abs(p_delta)/p_scale,
        "cd_viscous_a": certificate_a["cdv"],
        "cd_viscous_b": certificate_b["cdv"],
        "relative_viscous_difference": abs(v_delta)/v_scale,
        "max_abs_delta_cp": max(r["max_abs_delta_cp"] for r in regions),
        "max_abs_delta_cf": max(r["max_abs_delta_cf"] for r in regions),
        "regions": regions,
        "gci": "NOT_PERFORMED_REQUIRES_THREE_CONVERGED_COMPARABLE_GRIDS",
        "experimental_validation": "NOT_PERFORMED",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ("a", "b"):
        parser.add_argument(f"--{key}-profiles", required=True, type=Path)
        parser.add_argument(f"--{key}-forces", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    try:
        result = compare(
            read_profile(args.a_profiles), read_profile(args.b_profiles),
            read_force(args.a_forces), read_force(args.b_forces),
        )
        args.output.parent.mkdir(parents=True, exist_ok=True)
        temporary = args.output.with_name(args.output.name + ".tmp")
        temporary.write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
        temporary.replace(args.output)
        print(json.dumps({k:v for k,v in result.items() if k != "regions"}, indent=2))
        return 0
    except (ValueError, TypeError, OSError, KeyError) as exc:
        print("NACA WALL RADIAL COMPARISON FAILED:", exc)
        args.output.unlink(missing_ok=True)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
