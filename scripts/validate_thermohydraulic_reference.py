"""Validate frozen thermohydraulic benchmark reference contracts.

This is a reference/metadata CI, not a CFDX solver qualification. Physical
validation remains BLOCKED until a production turbulence/energy path exists
and solver results are compared with the pinned reference data.
"""

from __future__ import annotations

import argparse
import hashlib
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
    require(
        bool(data["reference"].get("primary_url")), f"{case_id}: missing primary source"
    )
    require(
        data["status"].get("contract") == "READY_FOR_CONTRACT_CI",
        f"{case_id}: contract must be explicitly ready for contract CI",
    )
    require(
        data["status"].get("physical_validation", "").startswith("BLOCKED_"),
        f"{case_id}: physical validation must remain fail-closed/blocked",
    )

    derived: dict[str, Any] = {}
    if case_id == "TH-01":
        physics = data["physics"]
        require(
            physics["re_tau_targets"] == [180, 395, 590],
            "TH-01: Re_tau matrix mismatch",
        )
        require(
            physics["re_tau_nominal_targets"] == [180, 395, 590],
            "TH-01: nominal Re_tau labels mismatch",
        )
        require(len(data["reference"]["datasets"]) == 3, "TH-01: missing DNS datasets")
        require(
            data["reference"]["data_url"]
            == "https://turbulence.oden.utexas.edu/data/MKM/",
            "TH-01: official DNS data portal mismatch",
        )
        require(
            data["reference"]["profile_files"]
            == [
                f"chan{re}/profiles/chan{re}.{ext}"
                for re in (180, 395, 590)
                for ext in ("means", "reystress")
            ],
            "TH-01: official MKM relative profile paths mismatch",
        )
        required = {
            "mean profiles",
            "Reynolds-stress profiles",
            "spectra",
            "two-point correlations",
        }
        require(
            required.issubset(set(data["reference"]["data_contents"])),
            "TH-01: incomplete DNS data inventory",
        )
        require(
            "constant mass flux" in physics["forcing"],
            "TH-01: forcing convention must match Moser DNS",
        )
        derived["re_tau_cases"] = physics["re_tau_targets"]
        actual = [entry["re_tau_actual"] for entry in physics["reference_cases"]]
        require(
            actual == [178.12, 392.24, 587.19],
            "TH-01: actual DNS Re_tau values mismatch",
        )
        require(
            [entry["spectral_grid"] for entry in physics["reference_cases"]]
            == ["128x129x128", "256x193x192", "384x257x384"],
            "TH-01: DNS spectral grid inventory mismatch",
        )
        derived["actual_re_tau"] = actual
        derived["profile_files"] = data["reference"]["profile_files"]
        derived["profile_note"] = (
            "chan180.means header gives 178.12; nominal label is 180"
        )
        derived["solver_validation"] = (
            "BLOCKED: RANS/LES production coupling and exact DNS data ingestion"
        )

    elif case_id == "TH-03":
        physics = data["physics"]
        geom = physics["geometry"]
        flow = physics["flow"]
        inlet = physics["inlet_boundary_layer_at_x_over_h_minus4"]
        require(close(geom["step_height_m"], 0.0127), "TH-03: step height mismatch")
        require(
            close(geom["expansion_ratio"], 1.125), "TH-03: expansion ratio mismatch"
        )
        require(flow["reynolds_step"] == 37500, "TH-03: Re_h mismatch")
        require(
            close(flow["velocity_reference_m_s"], 44.2, atol=1e-9),
            "TH-03: U_ref mismatch",
        )
        require(close(inlet["delta_m"], 0.019), "TH-03: inlet delta mismatch")
        require(
            close(inlet["displacement_thickness_m"], 0.0027),
            "TH-03: delta-star mismatch",
        )
        require(close(inlet["momentum_thickness_m"], 0.0019), "TH-03: theta mismatch")
        require(
            close(data["qoi"]["uncertainty"]["mean_velocity_relative"], 0.015),
            "TH-03: mean-velocity uncertainty mismatch",
        )
        require(
            close(data["qoi"]["uncertainty"]["reynolds_stress_relative"], 0.12),
            "TH-03: Reynolds-stress uncertainty mismatch",
        )
        require(
            close(data["qoi"]["uncertainty"]["pressure_coefficient_absolute"], 0.0002),
            "TH-03: pressure-coefficient uncertainty mismatch",
        )
        require(
            close(data["qoi"]["reference_targets"]["reattachment_x_over_h"], 6.26),
            "TH-03: reattachment reference mismatch",
        )
        require(
            close(
                data["qoi"]["reference_targets"]["reattachment_uncertainty_absolute"],
                0.10,
            ),
            "TH-03: reattachment uncertainty mismatch",
        )
        derived["reattachment_x_over_h"] = 6.26
        derived["reattachment_uncertainty"] = 0.10
        # Independently recompute Re_h from the dimensional inputs, and H1.
        computed_re = (
            flow["velocity_reference_m_s"]
            * geom["step_height_m"]
            / flow["kinematic_viscosity_m2_s"]
        )
        require(
            close(computed_re, float(flow["reynolds_step"]), atol=0.1),
            f"TH-03: inconsistent Re_h versus nu/U/h ({computed_re})",
        )
        computed_shape = (
            inlet["displacement_thickness_m"] / inlet["momentum_thickness_m"]
        )
        require(
            close(computed_shape, inlet["shape_factor"], atol=0.01),
            "TH-03: inlet H1 != delta*/theta",
        )
        derived["reference_reynolds"] = computed_re
        derived["inlet_shape_factor_from_thicknesses"] = computed_shape
        derived["implied_kinematic_viscosity_m2_s"] = flow["kinematic_viscosity_m2_s"]
        derived["reattachment_gate"] = (
            "zero crossing of tau_w/Cf from negative to positive; value loaded from CMP32 data, not guessed"
        )

    elif case_id == "TH-07":
        dim = data["physics"]["dimensionless"]
        require(dim["bulk_reynolds"] == 4328, "TH-07: bulk Reynolds number mismatch")
        require(close(dim["prandtl"], 0.71), "TH-07: Prandtl number mismatch")
        require(
            dim["grashof_cases"] == [0, 640000, 950000, 1600000],
            "TH-07: Grashof matrix mismatch",
        )
        re_b = float(dim["bulk_reynolds"])
        gr_values = [float(v) for v in dim["grashof_cases"]]
        ri_values = [float(v) for v in dim["richardson_cases"]]
        calculated = [g / (re_b * re_b) for g in gr_values]
        require(
            len(calculated) == len(ri_values), "TH-07: Richardson list length mismatch"
        )
        for index, (actual, expected) in enumerate(
            zip(calculated, ri_values, strict=True)
        ):
            require(
                close(actual, expected, atol=5.0e-7),
                f"TH-07: Ri mismatch at Gr index {index}: {actual} != {expected}",
            )
        targets = data["reference_targets"]["wall_shear_stress_mean"]
        require(
            close(targets["Gr_6.4e5"]["heated"], 0.0047),
            "TH-07: heated-wall tau mismatch",
        )
        require(
            close(targets["Gr_6.4e5"]["cooled"], 0.0041),
            "TH-07: cooled-wall tau mismatch",
        )
        require(
            close(targets["Gr_1.6e6"]["heated"], 0.0059),
            "TH-07: high-Gr heated tau mismatch",
        )
        require(
            close(targets["Gr_1.6e6"]["cooled"], 0.0039),
            "TH-07: high-Gr cooled tau mismatch",
        )
        require(
            close(data["physics"]["geometry"]["channel_height_over_delta"], 1.0),
            "TH-07: channel gap normalization mismatch",
        )
        derived["richardson_numbers"] = calculated
        derived["reference_bulk_reynolds"] = re_b
        derived["reference_prandtl"] = dim["prandtl"]
        # Check the actual DNS table rather than comparing only hard-coded
        # copies of the wall-shear values. Re_tau scales with sqrt(tau_wall).
        friction = data["reference_targets"]["friction_reynolds"]
        baseline_re = float(friction["Gr_0"]["global"])
        require(
            close(baseline_re, 284.02, atol=0.005),
            "TH-07: incorrect isothermal friction Reynolds baseline",
        )
        baseline_tau = float(targets["Gr_0"]["reference_isothermal"])
        require(baseline_tau > 0.0, "TH-07: invalid base wall shear")
        percent_changes = {}
        for name, gr in (
            ("Gr_6.4e5", 640000),
            ("Gr_9.5e5", 950000),
            ("Gr_1.6e6", 1600000),
        ):
            tabulated = friction[name]
            pair = targets[name]
            heated = float(tabulated["heated"])
            cooled = float(tabulated["cooled"])
            global_re = float(tabulated["global"])
            expected_global = math.sqrt((heated * heated + cooled * cooled) / 2)
            require(
                close(global_re, expected_global, atol=0.02),
                f"TH-07: global Re_tau is not consistent at Gr={gr}",
            )
            measured_ratio = float(pair["heated"]) / float(pair["cooled"])
            expected_ratio = (heated / cooled) ** 2
            require(
                math.isclose(measured_ratio, expected_ratio, rel_tol=0.05),
                f"TH-07: shear ratios disagree with friction Re at Gr={gr}",
            )
            percent_changes[str(gr)] = {
                "heated_friction_re_change_percent": 100 * (heated / baseline_re - 1),
                "cooled_friction_re_change_percent": 100 * (cooled / baseline_re - 1),
                "heated_shear_change_percent_from_rounded_table": 100
                * (float(pair["heated"]) / baseline_tau - 1),
                "cooled_shear_change_percent_from_rounded_table": 100
                * (float(pair["cooled"]) / baseline_tau - 1),
            }
        derived["corrected_wall_changes"] = percent_changes
        expected_width = [3.80, 3.85, 3.90, 4.10]
        require(
            data["reference_targets"]["grid_reference"]["mean_cell_width_star"]
            == expected_width,
            "TH-07: reference grid mean spacing mismatch",
        )

    elif case_id == "TH-08":
        physics = data["physics"]
        inlets = physics["inlets"]
        main = inlets["main_cold"]
        branch = inlets["branch_hot"]
        main_u = main["volumetric_flow_m3_s"] / (
            math.pi * main["diameter_m"] ** 2 / 4.0
        )
        branch_u = branch["volumetric_flow_m3_s"] / (
            math.pi * branch["diameter_m"] ** 2 / 4.0
        )
        require(
            close(main_u, physics["derived"]["bulk_velocity_main_m_s"], atol=2.0e-7),
            "TH-08: main bulk velocity mismatch",
        )
        require(
            close(
                branch_u, physics["derived"]["bulk_velocity_branch_m_s"], atol=2.0e-7
            ),
            "TH-08: branch bulk velocity mismatch",
        )
        require(
            close(branch["temperature_C"] - main["temperature_C"], 17.0),
            "TH-08: inlet delta-T mismatch",
        )
        require(
            close(physics["derived"]["mixed_temperature_star_fully_mixed"], 0.38),
            "TH-08: fully mixed T* target mismatch",
        )
        require(
            data["qoi"]["frequency"]["reported_temperature_spectral_peak_Hz"] == 3.5,
            "TH-08: reported spectral peak mismatch",
        )
        require(
            data["measurement_uncertainty"]["normalized_mean_temperature_absolute"]
            == 0.03,
            "TH-08: normalized temperature uncertainty mismatch",
        )
        require(
            data["qoi"]["velocity_stations"]["ldv_downstream_x_over_D"] == [2.6, 6.6],
            "TH-08: LDV station matrix mismatch",
        )
        require(
            data["qoi"]["velocity_stations"]["piv_downstream_x_over_D"]
            == [1.6, 2.6, 3.6, 4.6],
            "TH-08: PIV station matrix mismatch",
        )
        require(
            close(data["measurement_uncertainty"]["temperature_rms_relative"], 0.08),
            "TH-08: RMS temperature uncertainty mismatch",
        )
        derived["bulk_velocity_main_m_s"] = main_u
        derived["bulk_velocity_branch_m_s"] = branch_u
        volume_weighted = branch["volumetric_flow_m3_s"] / (
            main["volumetric_flow_m3_s"] + branch["volumetric_flow_m3_s"]
        )
        require(
            close(
                physics["derived"]["flow_ratio_main_over_branch"],
                main["volumetric_flow_m3_s"] / branch["volumetric_flow_m3_s"],
            ),
            "TH-08: main/branch volume-flow ratio mismatch",
        )
        require(
            close(
                physics["derived"][
                    "constant_property_volume_weighted_temperature_star"
                ],
                volume_weighted,
            ),
            "TH-08: constant-property volume mixing fraction mismatch",
        )
        require(
            close(
                physics["derived"]["delta_temperature_K"],
                branch["temperature_C"] - main["temperature_C"],
            ),
            "TH-08: inlet temperature difference mismatch",
        )
        derived["fully_mixed_temperature_star_published"] = 0.38
        derived["constant_property_volume_mixing_star"] = volume_weighted
        derived["reference_minus_constant_property_mixing"] = 0.38 - volume_weighted
        derived["mixing_note"] = (
            "0.38 is the value in OECD/NEA 2011 p.48, not the 0.40 "
            "equal-property volume-weighted inlet estimate; verify properties "
            "and mass/enthalpy balance before solver qualification"
        )

    return {
        "case_id": case_id,
        "contract_status": "PASS",
        "physical_validation_status": "BLOCKED_NOT_A_SOLVER_RUN",
        "reference": data["reference"]["citation"],
        "derived_checks": derived,
        "source_file": str(path.relative_to(ROOT)),
        "contract_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
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
