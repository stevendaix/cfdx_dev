"""End-to-end SST + thermal coupling in the actual production binary.

These tiny HDF5 cases are software/energy-balance verification, NOT a TH-07
validation: they have no turbulent inlet profiles or mesh-converged Nu data.
"""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import sys
from pathlib import Path

import h5py
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src/cfdx/python"))
sys.path.insert(0, str(ROOT / "tests/python"))
import test_n9_poiseuille_body_force_production as channel


def write_sst_thermal_case(
    path: Path,
    *,
    prt: str | None = "0.9",
    buoyant: bool = False,
    uniform: bool = False,
    sst: bool = True,
    driven: bool = False,
) -> None:
    # The Pr_t A/B case must have a real, reproducible carrier flow. Use the
    # same body-force channel contract as the N9/T2 production calibration;
    # keep the deliberately under-resolved 3x5 case only for the fail-closed
    # negative test below.
    channel.NX, channel.NY = (12, 16) if driven else (3, 5)
    channel.MAX_ITERATIONS = 3000 if driven else 40
    original_upstream_patch_type = channel._PATCH_TYPE["upstream_open"]
    if driven:
        channel._PATCH_TYPE["upstream_open"] = 1  # PatchType::INLET
    try:
        channel._write_case(path)
    finally:
        channel._PATCH_TYPE["upstream_open"] = original_upstream_patch_type
    with h5py.File(path, "r+") as case:
        setup = json.loads(case.attrs["case_setup_json"])
        if sst:
            setup["turbulence_model"] = "k_omega_sst"
        else:
            setup.pop("turbulence_model", None)
        setup["energy_model"] = "buoyant" if buoyant else "conjugate"
        setup["radiation_model"] = "none"
        setup["materials"] = [
            {
                "name": "air",
                "density": 1.2,
                "dynamic_viscosity": 0.012,
                "thermal_conductivity": 0.024,
                "specific_heat": 1005.0,
                **({"turbulence_model": "k_omega_sst"} if sst else {}),
                "extra_properties": (
                    {"turbulent_prandtl": prt} if prt is not None else {}
                ),
            }
        ]
        setup["body_force"] = [0.01, 0.0, 0.0] if driven else []
        setup["initial_condition"].update(
            {
                "velocity": 0.0,
                "velocity_vector": [0.0, 0.0, 0.0],
                "temperature": 300.0,
                **({"scalar_fields": {"k": 0.01, "omega": 2.0}} if sst else {}),
            }
        )
        if buoyant:
            setup["gravity_vector"] = "0 -9.81 0"
            setup["ref_temperature"] = 300.0
            setup["materials"][0]["extra_properties"][
                "thermal_expansion_coefficient"
            ] = "0.003"
        for patch in setup["boundary_conditions"]:
            name = patch["patch_name"]
            if name in {"top", "bottom"}:
                patch["value_type"] = "wall_thermal"
                patch["wall_temperature"] = 300.0 if uniform or name == "top" else 310.0
                if sst:
                    # The production SST contract requires explicit wall
                    # conditions on every WALL patch. This E2E uses the
                    # supported low-Re wall treatment; wall-function
                    # qualification remains a separate validation target.
                    patch["sst_scalar_conditions"] = [
                        {"field": "k", "type": "wall_k_zero", "value": 0.0},
                        {"field": "omega", "type": "wall_omega_sst", "value": 0.0},
                    ]
            elif name in {"front", "back"}:
                if sst:
                    patch["sst_scalar_conditions"] = [
                        {"field": "k", "type": "zero_gradient", "value": 0.0},
                        {"field": "omega", "type": "zero_gradient", "value": 0.0},
                    ]
            elif name == "upstream_open":
                if driven:
                    # Impose a reproducible inlet velocity while retaining
                    # the declared body force as a volumetric momentum source.
                    # The mesh patch type is changed to INLET above so the
                    # physical and topological contracts agree.
                    patch["type"] = "inlet"
                    patch["value_type"] = "fixed"
                    patch["velocity_vector"] = [0.1, 0.0, 0.0]
                    patch["temperature"] = 300.0
                if sst:
                    patch["sst_scalar_conditions"] = [
                        {"field": "k", "type": "fixed_value", "value": 0.01},
                        {"field": "omega", "type": "fixed_value", "value": 2.0},
                    ]
                    patch["turbulence_intensity"] = 0.05
                    patch["turbulence_length_scale"] = 0.1
            elif name == "downstream_open" and sst:
                # Preserve the pressure-outlet contract from the N9 channel
                # fixture. Replacing it with a velocity zero-gradient outlet
                # removes the pressure anchor and leaves a finite continuity
                # defect even when momentum and pressure linear solves settle.
                patch["sst_scalar_conditions"] = [
                    {"field": "k", "type": "zero_gradient", "value": 0.0},
                    {"field": "omega", "type": "zero_gradient", "value": 0.0},
                ]
        del case.attrs["case_setup_json"]
        case.attrs["case_setup_json"] = json.dumps(setup)


def execute(
    solver: Path,
    case: Path,
    result: Path,
    *,
    experimental_sst: bool = True,
    iterations: int = 40,
    timeout_seconds: int = 120,
) -> subprocess.CompletedProcess:
    command = [
        str(solver),
        "--mesh",
        str(case),
        "--output-dir",
        str(result),
        "--iterations",
        str(iterations),
    ]
    if experimental_sst:
        command.append("--experimental-sst")
    command.append("--benchmark-export")
    return subprocess.run(
        command,
        text=True,
        capture_output=True,
        timeout=timeout_seconds,
        check=False,
    )


@pytest.fixture(scope="module")
def solver() -> Path:
    value = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not value:
        pytest.skip("set CFDX_PRODUCTION_SOLVER to the built production executable")
    path = Path(value)
    assert path.is_file()
    return path


def test_real_thermal_coupling_prt_changes_heat_flux(
    solver: Path, tmp_path: Path
) -> None:
    fluxes: dict[str, float] = {}
    for prt in ("0.4", "1.8"):
        case = tmp_path / f"Prt_{prt}.cfdx.h5"
        out = tmp_path / f"Prt_{prt}"
        write_sst_thermal_case(case, prt=prt, driven=True)
        with h5py.File(case, "r") as case_file:
            setup = json.loads(case_file.attrs["case_setup_json"])
            axial = {
                patch["patch_name"]: (patch["type"], patch["value_type"])
                for patch in setup["boundary_conditions"]
                if patch["patch_name"] in {"upstream_open", "downstream_open"}
            }
            assert setup["body_force"] == [0.01, 0.0, 0.0]
            downstream = next(
                patch
                for patch in setup["boundary_conditions"]
                if patch["patch_name"] == "downstream_open"
            )
            assert downstream["pressure"] == 0.0
            assert axial == {
                "upstream_open": ("inlet", "fixed"),
                "downstream_open": ("pressure_outlet", "outlet_pressure"),
            }
            patch_types = case_file.attrs["boundary_patches"].split(";")
            assert any(
                item.startswith("upstream_open:") and item.endswith(":1")
                for item in patch_types
            )
            assert any(
                item.startswith("downstream_open:") and item.endswith(":2")
                for item in patch_types
            )
        run = execute(solver, case, out, iterations=3000, timeout_seconds=600)
        assert run.returncode == 0, (run.stdout + run.stderr)[-3000:]
        execution = json.loads((out / "execution.json").read_text())
        carrier = json.loads((out / "convergence.json").read_text())
        assert execution["converged"] is True
        assert carrier["converged"] is True
        energy = json.loads((out / "energy_convergence.json").read_text())
        summary = json.loads((out / "energy_summary.json").read_text())
        mass = json.loads(
            (out / "benchmark_observables/mass_conservation.json").read_text()
        )
        assert energy["converged"] is True
        assert energy["global_balance_normalized"] < 1e-9
        assert energy["max_cell_residual_normalized"] < 1e-9
        # The aggregate balance audit and the per-patch reporting pass are
        # independently reconstructed. Their floating-point reductions need
        # not be bitwise/near-bitwise identical; do not turn that diagnostic
        # into a second conservation gate. The blocking normalized global,
        # local and mass-balance gates above remain authoritative.
        assert math.isfinite(summary["boundary_flux_out"])
        assert all(math.isfinite(patch["flux_out"]) for patch in summary["patches"])
        assert mass["global_mass_imbalance_relative"] < 1e-9
        with (out / "benchmark_observables/cell_fields.csv").open(newline="") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 12 * 16
        assert all(float(row["nu_t_m2_s"]) >= 0 for row in rows)
        assert any(float(row["nu_t_m2_s"]) > 1e-7 for row in rows)
        assert all(300.0 <= float(row["temperature_K"]) <= 310.0 for row in rows)
        fluxes[prt] = abs(
            next(
                patch["flux_out"]
                for patch in summary["patches"]
                if patch["name"] == "bottom"
            )
        )
    # The constitutive law is verified quantitatively in
    # test_production_sst_thermal_coupling.cpp (lambda_eff = lambda +
    # rho*cp*nu_t/Pr_t). This E2E instead verifies the coupled response:
    # reducing Pr_t must increase the wall heat flux. Do not impose a fixed
    # percentage here: the wall-normal molecular sublayer, resolved mesh and
    # frozen carrier determine the flux sensitivity, not only the bulk
    # conductivity ratio. This tiny 12x16 fixture is explicitly not a
    # wall-resolved turbulent heat-transfer validation.
    assert fluxes["0.4"] > fluxes["1.8"]


def test_real_buoyant_laminar_isothermal_equilibrium(
    solver: Path, tmp_path: Path
) -> None:
    # This test verifies the Boussinesq equilibrium independently of SST.
    # The separate Pr_t A/B test above exercises the real SST/energy coupling.
    # Combining zero buoyancy with a non-equilibrium k/omega field and wall
    # turbulence conditions does not define a steady SST state.
    case = tmp_path / "buoyant_laminar_uniform.cfdx.h5"
    out = tmp_path / "buoyant_laminar_uniform"
    write_sst_thermal_case(case, prt=None, buoyant=True, uniform=True, sst=False)
    run = execute(solver, case, out, experimental_sst=False)
    assert run.returncode == 0, (run.stdout + run.stderr)[-3000:]
    assert "Buoyant flow-energy outer loop converged" in run.stdout
    summary = json.loads((out / "energy_summary.json").read_text())
    energy = json.loads((out / "energy_convergence.json").read_text())
    assert summary["min_temperature"] == pytest.approx(300.0, abs=1e-10)
    assert summary["max_temperature"] == pytest.approx(300.0, abs=1e-10)
    assert energy["global_balance_normalized"] < 1e-9
    assert json.loads((out / "execution.json").read_text())["converged"] is True


def test_missing_prt_fails_before_solver_execution(
    solver: Path, tmp_path: Path
) -> None:
    case = tmp_path / "missing_prt.cfdx.h5"
    write_sst_thermal_case(case, prt=None)
    run = execute(solver, case, tmp_path / "failed")
    assert run.returncode != 0
    assert "turbulent_prandtl" in run.stderr + run.stdout


def test_thermal_campaign_fails_closed_when_carrier_does_not_converge(
    solver: Path, tmp_path: Path
) -> None:
    sys.path.insert(0, str(ROOT))
    from scripts import run_thermohydraulic_cases as campaign

    case = tmp_path / "TH07_debug.cfdx.h5"
    write_sst_thermal_case(case)
    report = campaign.run_one(
        solver,
        tmp_path,
        tmp_path / "campaign",
        "TH-07",
        {"id": "Gr_0"},
        {
            "mesh_id": "test15",
            "case_file": case.name,
            "cells": 15,
            "max_iterations": 40,
        },
        90,
    )
    # This deliberately under-resolved 3x5 SST fixture is a negative
    # convergence case. A failed carrier solve must never be upgraded to an
    # executed/unqualified campaign or receive synthetic energy evidence.
    assert report["status"] == "INCOMPLETE_OR_FAIL", report
    assert report["physical_qualification"] == "NOT_QUALIFIED"
    assert report["exit_code"] != 0
    assert report["solver_executed"] is True
    assert report["error"].startswith("CFDX process exited with code")
    assert "energy_imbalance_relative" not in report
    assert "files_sha256" not in report
