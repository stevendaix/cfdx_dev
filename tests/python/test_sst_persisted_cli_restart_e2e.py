"""Production CLI regression for persisted SST cases and DAT restart equivalence.

The fixture starts from a real converted CFDX case so the C++ case reader,
boundary metadata reader, and production executable are exercised together.
This is software-integration evidence, not a physical SST qualification.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

from cfdx.dat_io import read_dat_restart


def _prepare_persisted_sst_case(case_path: Path) -> None:
    """Declare explicit SST initial data and patch-level scalar BCs in a case."""
    with h5py.File(case_path, "r+") as h5:
        raw = h5.attrs["case_setup_json"]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        setup = json.loads(str(raw))
        setup["physics_model"] = "incompressible_laminar"
        setup["turbulence_model"] = "k_omega_sst"
        setup["energy_model"] = "none"
        setup["radiation_model"] = "none"
        setup["multiphase_model"] = "none"
        setup["transient"] = False
        setup["solver_mode"] = "steady"
        setup["materials"] = [
            {
                "name": "air",
                "density": 1.2,
                "dynamic_viscosity": 0.012,
                "turbulence_model": "k_omega_sst",
            }
        ]
        setup["initial_condition"] = {
            "velocity": 1.0,
            "pressure": 0.0,
            "velocity_vector": [1.0, 0.0, 0.0],
            "scalar_fields": {"k": 0.01, "omega": 2.0},
        }
        # Preserve the converter's resolved N1 selection; the test must not
        # bypass the production numerics gate.
        for boundary in setup.get("boundary_conditions", []):
            patch_type = str(boundary.get("type", "")).lower()
            value_type = str(boundary.get("value_type", "")).lower()
            if patch_type == "wall" or value_type in {"wall_no_slip", "wall"}:
                boundary["sst_scalar_conditions"] = [
                    {"field": "k", "type": "wall_k_zero", "value": 0.0},
                    {"field": "omega", "type": "wall_omega_sst", "value": 0.0},
                ]
            else:
                boundary["sst_scalar_conditions"] = [
                    {"field": "k", "type": "fixed_value", "value": 0.01},
                    {"field": "omega", "type": "fixed_value", "value": 2.0},
                ]
        h5.attrs.modify("case_setup_json", json.dumps(setup))


def _run(
    solver: Path,
    case: Path,
    output: Path,
    iterations: int,
    restart: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    command = [
        str(solver),
        "--mesh",
        str(case),
        "--output-dir",
        str(output),
        "--iterations",
        str(iterations),
        "--experimental-sst",
    ]
    if restart is not None:
        command.extend(["--restart", str(restart)])
    return subprocess.run(
        command, capture_output=True, text=True, timeout=120, check=False
    )


def _assert_close_field(left, right, name: str) -> None:
    assert left.dimension == right.dimension, f"{name}: dimension changed"
    a = np.asarray(left.values, dtype=float)
    b = np.asarray(right.values, dtype=float)
    assert a.shape == b.shape, f"{name}: shape changed ({a.shape} != {b.shape})"
    assert np.isfinite(a).all() and np.isfinite(b).all(), f"{name}: non-finite state"
    # The restart is a serialization boundary, not a new discretization.
    np.testing.assert_allclose(a, b, rtol=2e-10, atol=2e-12, err_msg=name)


def test_persisted_sst_cli_wall_omega_and_dat_restart_equivalence(
    tmp_path: Path,
) -> None:
    solver_text = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver_text:
        pytest.skip("CFDX_PRODUCTION_SOLVER is provided by CTest")
    solver = Path(solver_text)
    assert solver.is_file()

    root = Path(__file__).resolve().parents[2]
    source = root / "tests" / "data" / "su2" / "mesh_NACA0012_inv.su2"
    config = root / "tests" / "data" / "su2" / "inv_NACA0012_basic.cfg"
    assert source.is_file() and config.is_file()

    case_dir = tmp_path / "persisted_case"
    case_dir.mkdir()
    source_case = case_dir / source.name
    config_case = case_dir / source.with_suffix(".cfg").name
    source_case.write_bytes(source.read_bytes())
    config_case.write_bytes(config.read_bytes())
    case_path = case_dir / "case.cfdx.h5"

    from cfdx.io.converter import convert

    converted = convert(source_case, output=case_path)
    assert converted.case is not None
    _prepare_persisted_sst_case(case_path)
    # Prove that the contract is stored in the HDF5 file, not only in memory.
    with h5py.File(case_path, "r") as h5:
        raw = h5.attrs["case_setup_json"]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        setup = json.loads(str(raw))
        assert setup["turbulence_model"] == "k_omega_sst"
        assert any(
            condition.get("type") == "wall_omega_sst"
            for boundary in setup["boundary_conditions"]
            for condition in boundary.get("sst_scalar_conditions", [])
        ), "fixture must persist at least one WALL_OMEGA_SST obligation"

    continuous = _run(solver, case_path, tmp_path / "continuous", 2)
    assert continuous.returncode == 0, (continuous.stdout + continuous.stderr)[-16000:]
    continuous_dat = tmp_path / "continuous" / "restart.dat"
    assert continuous_dat.is_file()
    continuous_state = read_dat_restart(continuous_dat)
    assert {"U", "p", "k", "omega"}.issubset(continuous_state.fields)

    first = _run(solver, case_path, tmp_path / "first_segment", 1)
    assert first.returncode == 0, (first.stdout + first.stderr)[-16000:]
    first_dat = tmp_path / "first_segment" / "restart.dat"
    assert first_dat.is_file()
    first_state = read_dat_restart(first_dat)
    assert {"U", "p", "k", "omega"}.issubset(first_state.fields)

    resumed = _run(
        solver, case_path, tmp_path / "resumed_segment", 2, restart=first_dat
    )
    assert resumed.returncode == 0, (resumed.stdout + resumed.stderr)[-16000:]
    resumed_dat = tmp_path / "resumed_segment" / "restart.dat"
    assert resumed_dat.is_file()
    resumed_state = read_dat_restart(resumed_dat)
    assert resumed_state.cells == continuous_state.cells
    for name in ("U", "p", "k", "omega"):
        _assert_close_field(
            continuous_state.fields[name], resumed_state.fields[name], name
        )
