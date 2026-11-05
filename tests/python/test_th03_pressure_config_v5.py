"""V5: independent pressure settings never silently alter momentum/physics gates."""

import json
from pathlib import Path

import h5py
import pytest

from scripts.generate_th03_openfoam_blockmesh import (
    make_case_setup,
    make_structured_mesh,
    write_case,
)


def test_independent_pressure_settings_preserved(tmp_path: Path):
    output = tmp_path / "p.cfdx.h5"
    mesh = make_structured_mesh("smoke")
    summary = write_case(
        output,
        mesh,
        mesh_source="six-block smoke fixture",
        pressure_krylov="cg",
        pressure_preconditioner="native_amg",
        pressure_tolerance=1e-11,
    )
    assert summary["pressure_krylov"] == "cg"
    with h5py.File(output) as handle:
        setup = json.loads(handle.attrs["case_setup_json"])
    settings = setup["numerics"]["raw_settings"]
    assert settings["pressure_krylov"] == "cg"
    assert settings["pressure_preconditioner"] == "native_amg"
    assert settings["pressure_linear_tolerance"] == "1e-11"
    assert settings["pressure_gmres_restart"] == "160"
    assert setup["numerics"]["selection"]["entries"][-1]["family"] == "linear_solver"


@pytest.mark.parametrize("tol", [0.0, -1.0, float("nan"), float("inf")])
def test_pressure_target_rejects_invalid_values(tol):
    with pytest.raises(ValueError, match="pressure tolerance"):
        make_case_setup(2, 0.25, 0.1, pressure_tolerance=tol)


def test_three_grid_campaign_fails_closed_with_rejected_solver(tmp_path: Path):
    import sys

    root = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(root / "scripts"))
    from run_th03_pressure_campaign_v5 import campaign

    fake = tmp_path / "reject.sh"
    fake.write_text(
        '#!/bin/sh\necho "solve_steady_incompressible: pressure-correction solve did not converge (status=1, iterations=108, residual=1e-12, relative=2e-11)"\nexit 2\n'
    )
    fake.chmod(0o755)
    summary = campaign(fake, tmp_path / "result", 5.0, 1e-11, ("smoke",))
    assert summary["results"][0]["return_code"] == 2
    assert summary["results"][0]["physical_qualification"] == "NOT_QUALIFIED"
    assert summary["results"][0]["last_pressure_failure"]["relative_residual"] == 2e-11
    assert summary["spatial_convergence"] == "NOT_AVAILABLE"
