"""End-to-end SST + thermal coupling in the actual production binary.

These tiny HDF5 cases are software/energy-balance verification, NOT a TH-07
validation: they have no turbulent inlet profiles or mesh-converged Nu data.
"""

from __future__ import annotations

import csv
import json
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
    path: Path, *, prt: str | None = "0.9", buoyant: bool = False, uniform: bool = False
) -> None:
    channel.NX, channel.NY = 3, 5
    channel.MAX_ITERATIONS = 40
    channel._write_case(path)
    with h5py.File(path, "r+") as case:
        setup = json.loads(case.attrs["case_setup_json"])
        setup["turbulence_model"] = "k_omega_sst"
        setup["energy_model"] = "buoyant" if buoyant else "conjugate"
        setup["radiation_model"] = "none"
        setup["materials"] = [
            {
                "name": "air",
                "density": 1.2,
                "dynamic_viscosity": 0.012,
                "thermal_conductivity": 0.024,
                "specific_heat": 1005.0,
                "turbulence_model": "k_omega_sst",
                "extra_properties": (
                    {"turbulent_prandtl": prt} if prt is not None else {}
                ),
            }
        ]
        setup["body_force"] = []
        setup["initial_condition"].update(
            {
                "velocity": 0.0,
                "velocity_vector": [0.0, 0.0, 0.0],
                "temperature": 300.0,
                "scalar_fields": {"k": 0.01, "omega": 2.0},
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
                patch["sst_scalar_conditions"] = [
                    {"field": "k", "type": "wall_k_zero", "value": 0.0},
                    {"field": "omega", "type": "wall_omega_sst", "value": 0.0},
                ]
            elif name in {"front", "back"}:
                patch["sst_scalar_conditions"] = [
                    {"field": "k", "type": "zero_gradient", "value": 0.0},
                    {"field": "omega", "type": "zero_gradient", "value": 0.0},
                ]
            elif name == "outlet":
                # A prescribed-velocity inlet needs a pressure reference at
                # the outlet; a zero-gradient pressure boundary leaves the
                # all-Neumann pressure system unable to balance the imposed
                # mass flux on this smoke-test channel.
                patch["type"] = "pressure_outlet"
                patch["value_type"] = "outlet_pressure"
                patch["pressure"] = 0.0
                patch["sst_scalar_conditions"] = [
                    {"field": "k", "type": "fixed_value", "value": 0.01},
                    {"field": "omega", "type": "fixed_value", "value": 2.0},
                ]
            else:
                patch["sst_scalar_conditions"] = [
                    {"field": "k", "type": "fixed_value", "value": 0.01},
                    {"field": "omega", "type": "fixed_value", "value": 2.0},
                ]
                if name == "inlet":
                    patch["value_type"] = "fixed"
                    patch["turbulence_intensity"] = 0.05
                    patch["turbulence_length_scale"] = 0.1
                    # This is an energy/SST coupling smoke test, not a flow
                    # qualification. Keep the carrier at rest so the
                    # independent local face-momentum gate is exact; the
                    # nonzero k/omega field still exercises turbulent lambda_eff.
                    patch["velocity_magnitude"] = 0.0
                    patch["velocity_vector"] = [0.0, 0.0, 0.0]
        del case.attrs["case_setup_json"]
        case.attrs["case_setup_json"] = json.dumps(setup)


def execute(solver: Path, case: Path, result: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [
            str(solver),
            "--mesh",
            str(case),
            "--output-dir",
            str(result),
            "--iterations",
            "40",
            "--experimental-sst",
            "--benchmark-export",
        ],
        text=True,
        capture_output=True,
        timeout=120,
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
        write_sst_thermal_case(case, prt=prt)
        run = execute(solver, case, out)
        assert run.returncode == 0, (run.stdout + run.stderr)[-3000:]
        energy = json.loads((out / "energy_convergence.json").read_text())
        summary = json.loads((out / "energy_summary.json").read_text())
        mass = json.loads(
            (out / "benchmark_observables/mass_conservation.json").read_text()
        )
        assert energy["converged"] is True
        assert energy["global_balance_normalized"] < 1e-9
        assert energy["max_cell_residual_normalized"] < 1e-9
        assert abs(summary["boundary_flux_out"]) < 1e-9
        assert mass["global_mass_imbalance_relative"] < 1e-9
        with (out / "benchmark_observables/cell_fields.csv").open(newline="") as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 15
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
    # With the same strictly positive nu_t, a lower Pr_t gives a larger
    # turbulent lambda_eff; proving the energy stage consumes the SST field.
    assert fluxes["0.4"] > fluxes["1.8"] * 1.15


def test_real_buoyant_sst_isothermal_equilibrium(solver: Path, tmp_path: Path) -> None:
    case = tmp_path / "buoyant_uniform.cfdx.h5"
    out = tmp_path / "buoyant_uniform"
    write_sst_thermal_case(case, buoyant=True, uniform=True)
    run = execute(solver, case, out)
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


def test_thermal_campaign_requires_native_energy_audit(
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
    assert report["status"] == "EXECUTED_UNQUALIFIED", report
    assert report["physical_qualification"] == "NOT_QUALIFIED"
    assert report["energy_imbalance_relative"] < 1e-9
    assert report["max_cell_energy_imbalance_relative"] < 1e-9
    assert "energy_convergence.json" in report["files_sha256"]
    assert "energy_summary.json" in report["files_sha256"]
