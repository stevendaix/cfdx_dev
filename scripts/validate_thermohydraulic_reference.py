#!/usr/bin/env python3
"""Validate frozen thermohydraulic benchmark reference contracts.

This is a reference/metadata CI, not a CFDX solver qualification. Physical
validation remains BLOCKED until a production turbulence/energy path exists
and solver results are compared with the pinned reference data.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
CASE_DIR = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
CASE_IDS = ("TH-01", "TH-03", "TH-07", "TH-08")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def close(actual: float, expected: float, *, atol: float = 1.0e-6) -> bool:
    return math.isclose(actual, expected, rel_tol=0.0, abs_tol=atol)


def validate_case(case_id: str) -> dict[str, Any]:
    path = CASE_DIR / f"{case_id}.json"
    data = json.loads(path.read_text(encoding="utf-8"))
    require(data.get("schema_version") == 1, f"{case_id}: unsupported schema_version")
    require(data.get("case_id") == case_id, f"{case_id}: case_id mismatch")
    for key in ("reference", "physics", "qoi", "status"):
        require(key in data, f"{case_id}: missing top-level {key}")
    require(bool(data["reference"].get("citation")), f"{case_id}: missing citation")
    require(bool(data["reference"].get("primary_url")), f"{case_id}: missing primary source")
    require(data["status"].get("contract") == "READY_FOR_CONTRACT_CI",
            f"{case_id}: contract must be explicitly ready for contract CI")
    require(data["status"].get("physical_validation", "").startswith("BLOCKED_"),
            f"{case_id}: physical validation must remain fail-closed/blocked")

    derived: dict[str, Any] = {}
    if case_id == "TH-01":
        physics = data["physics"]
        require(physics["re_tau_targets"] == [180, 395, 590], "TH-01: Re_tau matrix mismatch")
        require(len(data["reference"]["datasets"]) == 3, "TH-01: missing DNS datasets")
        required = {"mean profiles", "Reynolds-stress profiles", "spectra", "two-point correlations"}
        require(required.issubset(set(data["reference"]["data_contents"])),
                "TH-01: incomplete DNS data inventory")
        require("constant bulk mass flux" in physics["forcing"],
                "TH-01: forcing convention must match Moser DNS")
        derived["re_tau_cases"] = physics["re_tau_targets"]
        derived["solver_validation"] = "BLOCKED: RANS/LES production coupling and exact DNS data ingestion"

    elif case_id == "TH-03":
        physics = data["physics"]
        geom = physics["geometry"]
        flow = physics["flow"]
        inlet = physics["inlet_boundary_layer_at_x_over_h_minus4"]
        require(close(geom["step_height_m"], 0.0127), "TH-03: step height mismatch")
        require(close(geom["expansion_ratio"], 1.125), "TH-03: expansion ratio mismatch")
        require(flow["reynolds_step"] == 37500, "TH-03: Re_h mismatch")
        require(close(flow["velocity_reference_m_s"], 44.2, atol=1e-9),
                "TH-03: U_ref mismatch")
        require(close(inlet["delta_m"], 0.019), "TH-03: inlet delta mismatch")
        require(close(inlet["displacement_thickness_m"], 0.0027), "TH-03: delta-star mismatch")
        require(close(inlet["momentum_thickness_m"], 0.0019), "TH-03: theta mismatch")
        require(close(data["qoi"]["uncertainty"]["mean_velocity_relative"], 0.015),
                "TH-03: mean-velocity uncertainty mismatch")
        require(close(data["qoi"]["uncertainty"]["reynolds_stress_relative"], 0.12),
                "TH-03: Reynolds-stress uncertainty mismatch")
        require(close(data["qoi"]["reference_targets"]["reattachment_x_over_h"], 6.26),
                "TH-03: reattachment reference mismatch")
        require(close(data["qoi"]["reference_targets"]["reattachment_uncertainty_absolute"], 0.10),
                "TH-03: reattachment uncertainty mismatch")
        derived["reattachment_x_over_h"] = 6.26
        derived["reattachment_uncertainty"] = 0.10
        derived["reference_reynolds"] = flow["reynolds_step"]
        derived["reattachment_gate"] = "zero crossing of tau_w/Cf from negative to positive; value loaded from CMP32 data, not guessed"

    elif case_id == "TH-07":
        dim = data["physics"]["dimensionless"]
        re_b = float(dim["bulk_reynolds"])
        gr_values = [float(v) for v in dim["grashof_cases"]]
        ri_values = [float(v) for v in dim["richardson_cases"]]
        calculated = [g / (re_b * re_b) for g in gr_values]
        require(len(calculated) == len(ri_values), "TH-07: Richardson list length mismatch")
        for index, (actual, expected) in enumerate(zip(calculated, ri_values, strict=True)):
            require(close(actual, expected, atol=5.0e-7),
                    f"TH-07: Ri mismatch at Gr index {index}: {actual} != {expected}")
        targets = data["reference_targets"]["wall_shear_stress_mean"]
        require(close(targets["Gr_6.4e5"]["heated"], 0.0047), "TH-07: heated-wall tau mismatch")
        require(close(targets["Gr_6.4e5"]["cooled"], 0.0041), "TH-07: cooled-wall tau mismatch")
        require(close(targets["Gr_1.6e6"]["heated"], 0.0059), "TH-07: high-Gr heated tau mismatch")
        require(close(targets["Gr_1.6e6"]["cooled"], 0.0039), "TH-07: high-Gr cooled tau mismatch")
        require(close(data["physics"]["geometry"]["channel_height_over_delta"], 1.0),
                "TH-07: channel gap normalization mismatch")
        derived["richardson_numbers"] = calculated
        derived["reference_bulk_reynolds"] = re_b
        derived["reference_prandtl"] = dim["prandtl"]

    elif case_id == "TH-08":
        physics = data["physics"]
        inlets = physics["inlets"]
        main = inlets["main_cold"]
        branch = inlets["branch_hot"]
        main_u = main["volumetric_flow_m3_s"] / (math.pi * main["diameter_m"] ** 2 / 4.0)
        branch_u = branch["volumetric_flow_m3_s"] / (math.pi * branch["diameter_m"] ** 2 / 4.0)
        require(close(main_u, physics["derived"]["bulk_velocity_main_m_s"], atol=2.0e-7),
                "TH-08: main bulk velocity mismatch")
        require(close(branch_u, physics["derived"]["bulk_velocity_branch_m_s"], atol=2.0e-7),
                "TH-08: branch bulk velocity mismatch")
        require(close(branch["temperature_C"] - main["temperature_C"], 17.0),
                "TH-08: inlet delta-T mismatch")
        require(close(physics["derived"]["mixed_temperature_star_fully_mixed"], 0.38),
                "TH-08: fully mixed T* target mismatch")
        require(data["qoi"]["frequency"]["reported_temperature_spectral_peak_Hz"] == 3.5,
                "TH-08: reported spectral peak mismatch")
        require(data["measurement_uncertainty"]["normalized_mean_temperature_absolute"] == 0.03,
                "TH-08: normalized temperature uncertainty mismatch")
        derived["bulk_velocity_main_m_s"] = main_u
        derived["bulk_velocity_branch_m_s"] = branch_u
        derived["fully_mixed_temperature_star"] = 0.38

    return {
        "case_id": case_id,
        "contract_status": "PASS",
        "physical_validation_status": "BLOCKED_NOT_A_SOLVER_RUN",
        "reference": data["reference"]["citation"],
        "derived_checks": derived,
        "source_file": str(path.relative_to(ROOT)),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=CASE_IDS, action="append")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    selected = args.case or list(CASE_IDS)
    results = [validate_case(case_id) for case_id in selected]
    report = {
        "schema_version": 1,
        "scope": "reference-contract integrity and dimensionless-number checks only",
        "not_a_solver_validation": True,
        "cases": results,
    }
    rendered = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
