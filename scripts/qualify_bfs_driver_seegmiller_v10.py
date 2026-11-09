"""Fail-closed BFS Driver-Seegmiller postprocess of a *converged* CFDX run.

This extractor is not a physical oracle: first-order wall shear in current
benchmark_observables is diagnostic only, and the present CSV lacks the
Reynolds shear stresses required for a full NASA validation.  Refuse to
compute a reattachment length from transient/rejected/unfinished states.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
import math
from pathlib import Path

from qualify_bfs_schur_evidence import validate_bfs_evidence

EXPERIMENTAL_XR_H = 6.26
EXPERIMENTAL_XR_UNCERTAINTY = 0.10
PROFILES_X_H = (1.0, 4.0, 6.0, 10.0)


def _number(cell: dict, name: str) -> float:
    v = float(cell[name])
    if not math.isfinite(v):
        raise ValueError(f"nonfinite {name}")
    return v


def cf_and_reattachment(
    faces_csv: Path, *, rho: float, u_ref: float, step_height: float
) -> dict:
    """Return strictly diagnostic Cf and zero-crossing; never a physical PASS."""
    if not (
        math.isfinite(rho)
        and rho > 0
        and math.isfinite(u_ref)
        and u_ref > 0
        and math.isfinite(step_height)
        and step_height > 0
    ):
        raise ValueError("invalid physical reference units")
    rows: list[tuple[float, float]] = []
    methods = set()
    with faces_csv.open(newline="", encoding="utf-8") as stream:
        for item in csv.DictReader(stream):
            if item["patch"] != "lowerWall":
                continue
            x = _number(item, "x_m") / step_height
            # The step's vertical wall belongs to lowerWall but its tau_x is
            # NOT streamwise bottom-wall skin friction.  Keep only horizontal
            # downstream wall faces with a predominantly vertical normal.
            if x < 0 or abs(_number(item, "area_y_m2")) < 0.95 * _number(
                item, "area_m2"
            ):
                continue
            method = item["wall_shear_method"]
            if method == "not_available_non_dirichlet_wall":
                raise ValueError("lower wall lacks a resolved shear boundary condition")
            shear = _number(item, "wall_shear_x_Pa")
            cf = shear / (0.5 * rho * u_ref * u_ref)
            rows.append((x, cf))
            methods.add(method)
    rows.sort(key=lambda item: item[0])
    if len(rows) < 4 or len({x for x, _ in rows}) != len(rows):
        raise ValueError("insufficient distinct downstream wall faces for Cf")
    if any(x2 <= x1 for (x1, _), (x2, _) in itertools.pairwise(rows)):
        raise ValueError("nonmonotone downstream wall geometry")
    zero_crossings = []
    for (x0, y0), (x1, y1) in itertools.pairwise(rows):
        if y0 < 0.0 <= y1:
            crossing = x1 if y1 == 0.0 else x0 - y0 * (x1 - x0) / (y1 - y0)
            zero_crossings.append(crossing)
    if not zero_crossings:
        xr = None
    else:
        xr = zero_crossings[0]
    return {
        "cf_series": [{"x_over_h": x, "cf_diagnostic": y} for x, y in rows],
        "wall_shear_methods": sorted(methods),
        "reattachment_x_over_h_diagnostic": xr,
        "reattachment_error_vs_driver_seegmiller": (
            xr - EXPERIMENTAL_XR_H if xr is not None else None
        ),
        "experimental_xr_over_h": EXPERIMENTAL_XR_H,
        "experimental_uncertainty_x_over_h": EXPERIMENTAL_XR_UNCERTAINTY,
        "physical_certification": "NOT_QUALIFIED_FIRST_ORDER_SHEAR_AND_NO_GRID_GCI",
    }


def extract_velocity_profiles(
    cells_csv: Path, *, step_height: float, u_ref: float
) -> dict:
    data: list[tuple[float, float, float, float]] = []
    with cells_csv.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            data.append(
                (
                    _number(row, "x_m") / step_height,
                    _number(row, "y_m") / step_height,
                    _number(row, "u_x_m_s") / u_ref,
                    _number(row, "u_y_m_s") / u_ref,
                )
            )
    if not data:
        raise ValueError("no cellwise velocity fields")
    result = {}
    for station in PROFILES_X_H:
        # Select a unique computational x station, keeping all vertical cells;
        # the reported station offset is explicit, never masked as interpolation.
        x_station = min((x for x, _, _, _ in data), key=lambda x: abs(x - station))
        points = sorted(
            (y, ux, uy) for x, y, ux, uy in data if abs(x - x_station) < 1e-9
        )
        if len(points) < 2:
            raise ValueError("insufficient cells in profile station")
        result[f"x_over_h_{station:g}"] = {
            "requested": station,
            "actual_x_over_h": x_station,
            "samples": [
                {"y_over_h": y, "ux_over_uref": ux, "uy_over_uref": uy}
                for y, ux, uy in points
            ],
            "reference_comparison": "NOT_EVALUATED_NO_NASA_PROFILE_FILE",
        }
    return result


def make_report(
    *,
    output_dir: Path,
    case_path: Path,
    algorithm: str,
    schur: str | None,
    preconditioner: str | None,
    process_exit_code: int | None,
    rho: float,
    u_ref: float,
    step_height: float,
) -> dict:
    fingerprint = hashlib.sha256(case_path.read_bytes()).hexdigest()
    gate = validate_bfs_evidence(
        output_dir,
        process_exit_code,
        algorithm=algorithm,
        schur=schur,
        preconditioner=preconditioner,
        case_sha256=fingerprint,
        case_path=case_path,
    )
    result = {
        "schema": "CFDX-BFS-DRIVER-SEEGMILLER/v10",
        "case_sha256": fingerprint,
        "v8_nonlinear_certificate_passed": gate["nonlinear_v8_proven"],
        "gate_defects": gate["defects"],
        "experimental_source": "NASA TMR Driver-Seegmiller 2DBFS",
        "experimental_reattachment_x_over_h": EXPERIMENTAL_XR_H,
        "experimental_reattachment_uncertainty": EXPERIMENTAL_XR_UNCERTAINTY,
        "physical_validation": "NOT_QUALIFIED",
        "comparison": None,
        "reynolds_shear_stresses": "NOT_EXPORTED_FROM_PRODUCTION_SOLVER",
        "pressure_coefficient": "NOT_QUALIFIED_CELL_CENTRE_P_IS_NOT_WALL_FACE_PRESSURE",
        "grid_convergence": "NOT_EVALUATED_REQUIRES_3_NESTED_GRIDS",
    }
    if not gate["nonlinear_v8_proven"]:
        result["reason"] = "GLOBAL_U_P_K_OMEGA_NOT_CONVERGED_NO_QOI_COMPUTED"
        return result
    obs = output_dir / "benchmark_observables"
    result["comparison"] = {
        "cf_and_reattachment": cf_and_reattachment(
            obs / "face_fluxes.csv", rho=rho, u_ref=u_ref, step_height=step_height
        ),
        "velocity_profiles": extract_velocity_profiles(
            obs / "cell_fields.csv", step_height=step_height, u_ref=u_ref
        ),
    }
    result["reason"] = "DIAGNOSTIC_ONLY_NASA_PROFILES_AND_STRESSES_NOT_CHECKED"
    return result


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run", type=Path, required=True)
    p.add_argument("--case", type=Path, required=True)
    p.add_argument(
        "--algorithm", choices=("simple", "simplec", "coupled"), required=True
    )
    p.add_argument("--schur", default=None)
    p.add_argument("--preconditioner", default=None)
    p.add_argument("--exit-code", type=int, required=True)
    p.add_argument("--density", type=float, required=True)
    p.add_argument("--u-ref", type=float, required=True)
    p.add_argument("--step-height", type=float, default=0.0127)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    try:
        report = make_report(
            output_dir=args.run,
            case_path=args.case,
            algorithm=args.algorithm,
            schur=args.schur,
            preconditioner=args.preconditioner,
            process_exit_code=args.exit_code,
            rho=args.density,
            u_ref=args.u_ref,
            step_height=args.step_height,
        )
    except (ValueError, KeyError, OSError, ZeroDivisionError) as e:
        report = {
            "schema": "CFDX-BFS-DRIVER-SEEGMILLER/v10",
            "physical_validation": "NOT_QUALIFIED",
            "reason": f"invalid_evidence:{e}",
        }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(report["reason"])
    return 0 if report["physical_validation"] == "QUALIFIED" else 2


if __name__ == "__main__":
    raise SystemExit(main())
