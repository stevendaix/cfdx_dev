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
        # The source SU2 fixture is an implicit-time case. Remove only its
        # temporal selection because this executable is steady; retain the
        # converter's other resolved N1 selections and their validation gate.
        for key in ("numerics", "numerical_config"):
            numerics = setup.get(key)
            if not isinstance(numerics, dict):
                continue
            selection = numerics.get("selection")
            if not isinstance(selection, dict):
                continue
            entries = selection.get("entries", [])
            selection["entries"] = [
                entry for entry in entries if entry.get("family") != "temporal"
            ]
            required = selection.get("required_families")
            if isinstance(required, list):
                selection["required_families"] = [
                    family for family in required if family != "temporal"
                ]

            # The SU2 adapter records JST as unmapped; the production driver
            # intentionally rejects that conversion gap. This integration
            # fixture explicitly chooses a supported convection method rather
            # than weakening the driver's numerical-selection contract.
            if key == "numerics":
                entries = selection.setdefault("entries", [])
                if not any(entry.get("family") == "convection" for entry in entries):
                    entries.append(
                        {
                            "family": "convection",
                            "configuration_key": "numerics.convection.upwind",
                        }
                    )
                required = selection.setdefault("required_families", [])
                if "convection" not in required:
                    required.append("convection")
        # The steady production executable rejects temporal declarations from
        # both the legacy scheme and the resolved numerical-method report.
        numerics = setup.get("numerics")
        if isinstance(numerics, dict):
            numerics["transient_scheme"] = "steady"
        numerical_config = setup.get("numerical_config")
        if isinstance(numerical_config, dict):
            numerical_config["transient_scheme"] = "steady"
            entries = numerical_config.get("entries")
            if isinstance(entries, list):
                numerical_config["entries"] = [
                    entry for entry in entries if entry.get("family") != "temporal"
                ]
        numerical_report = setup.get("numerical_report")
        if isinstance(numerical_report, dict):
            resolved = numerical_report.get("resolved")
            if isinstance(resolved, list):
                numerical_report["resolved"] = [
                    entry for entry in resolved if entry.get("family") != "temporal"
                ]
        for boundary in setup.get("boundary_conditions", []):
            patch_type = str(boundary.get("type", "")).lower()
            value_type = str(boundary.get("value_type", "")).lower()
            if patch_type == "inlet":
                # The production SST path derives inlet k/omega from these
                # explicit, dimensionally consistent inputs. Defaults of zero
                # are intentionally rejected by the case contract.
                boundary["velocity_magnitude"] = 1.0
                boundary["velocity_vector"] = [1.0, 0.0, 0.0]
                boundary["turbulence_intensity"] = 0.05
                boundary["turbulence_length_scale"] = 0.1
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
        # Recreate the attribute because the serialized JSON may grow beyond
        # the converter's original HDF5 string-attribute width.
        del h5.attrs["case_setup_json"]
        h5.attrs["case_setup_json"] = json.dumps(setup)


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

    # Keep this integration fixture deliberately tiny.  The regression is
    # about persisted SST metadata and restart determinism, not airfoil physics;
    # the NACA mesh made each first SST iteration exceed the CI test budget.
    case_dir = tmp_path / "persisted_case"
    case_dir.mkdir()
    source_case = case_dir / "sst_restart_fixture.su2"
    config_case = case_dir / "sst_restart_fixture.cfg"
    source_case.write_text(
        """NDIME= 2
NELEM= 2
5 0 1 2 0
5 0 2 3 1
NPOIN= 4
0.0 0.0 0
1.0 0.0 1
1.0 1.0 2
0.0 1.0 3
NMARK= 4
MARKER_TAG= inlet
MARKER_ELEMS= 1
3 0 3
MARKER_TAG= outlet
MARKER_ELEMS= 1
3 1 2
MARKER_TAG= wall_bottom
MARKER_ELEMS= 1
3 0 1
MARKER_TAG= wall_top
MARKER_ELEMS= 1
3 2 3
""",
        encoding="utf-8",
    )
    config_case.write_text(
        """SOLVER= INC_RANS
KIND_TURB_MODEL= SST
INC_DENSITY_INIT= 1.2
INC_VELOCITY_INIT= (1.0, 0.0, 0.0)
MU_CONSTANT= 0.012
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= NO
NUM_METHOD_GRAD= GREEN_GAUSS
TIME_DISCRE_FLOW= EULER_IMPLICIT
ITER= 2
""",
        encoding="utf-8",
    )
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

    # --iterations is a per-invocation iteration budget, not a global target:
    # the resumed segment performs one more iteration after the first segment.
    resumed = _run(
        solver, case_path, tmp_path / "resumed_segment", 1, restart=first_dat
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
