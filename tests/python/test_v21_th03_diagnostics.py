"""V21: source and numeric output must agree; no fake physical qualification."""

import os
import subprocess
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
import pytest

from scripts.generate_th03_openfoam_blockmesh import (
    make_case_setup,
    make_structured_mesh,
    write_case,
)
from scripts.qualify_th03_v21 import diagnostic, wall_sign_change


def test_v21_momentum_pressure_adaptation_requires_v19():
    with pytest.raises(ValueError, match="requires V19"):
        make_case_setup(1, 0.55, 0.3, v21_adaptive_relaxation=True)
    setup = make_case_setup(
        1,
        0.55,
        0.3,
        v18_sparse_flux_pressure_jacobian=True,
        v19_full_pressure_flux_step=True,
        v21_adaptive_relaxation=True,
    )
    assert setup.numerics.raw_settings["v21_adaptive_relaxation"] == "true"
    assert setup.numerics.under_relaxation_momentum == 0.55


def test_reattachment_proxy_detects_downstream_crossing_and_rejects_bad_samples():
    x = np.array([1.0, 2.0, 3.0, 4.0])
    u = np.array([1.0, -2.0, -1.0, 1.0])
    assert wall_sign_change(x, u) == [3.5]
    with pytest.raises(ValueError, match="increase"):
        wall_sign_change(x[::-1], u)
    with pytest.raises(ValueError, match="nonfinite"):
        wall_sign_change(x, np.array([1.0, float("nan"), 3.0, 4.0]))


def test_v21_export_of_real_sst_fields_and_qualification_flag(tmp_path: Path):
    solver = os.getenv("CFDX_PRODUCTION_SOLVER_BINARY")
    if not solver:
        pytest.skip("set CFDX_PRODUCTION_SOLVER_BINARY to test real production SST VTU")
    case = tmp_path / "th03_v21.cfdx.h5"
    write_case(
        case,
        make_structured_mesh("smoke"),
        mesh_source="V21 SST VTU gate",
        max_iterations=1,
        alpha_u=0.55,
        alpha_p=0.3,
        pressure_tolerance=1e-10,
        sst_inner_iterations=10,
        sst_newton_sink=True,
        sst_omega_first=True,
        sst_provisional=True,
        sst_relative_gate=True,
        sst_molecular_wall_diffusion=True,
        v18_sparse_flux_pressure_jacobian=True,
        v19_full_pressure_flux_step=True,
    )
    out = tmp_path / "result"
    p = subprocess.run(
        [
            solver,
            "--mesh",
            str(case),
            "--output-dir",
            str(out),
            "--experimental-sst",
            "--iterations",
            "1",
        ],
        capture_output=True,
        text=True,
        check=False,
        timeout=120,
    )
    assert p.returncode in (0, 1), (p.stdout[-2500:], p.stderr[-2500:])
    path = out / "result_1.vtu"
    assert path.exists(), (p.stdout[-1000:], p.stderr[-2000:])
    tree = ET.parse(path)
    fields = {i.get("Name") for i in tree.findall(".//CellData/DataArray")}
    assert {
        "p",
        "u_x",
        "u_y",
        "u_z",
        "k",
        "omega",
        "nu_t",
        "wall_distance",
        "sst_strain_rate",
    } <= fields
    result = diagnostic(path, require_sst=True)
    assert result["physical_qualification"] == "NOT_QUALIFIED"
    assert result["sst"]["omega_min"] > 0
    assert result["sst"]["nut_over_nu_max"] >= 0
    assert abs(result["near_wall_table"][-1]["cp_shift_proxy"]) < 1e-13
    assert result["wall_cp_shift_anchor_x_h"] > 0
    assert result["wall_cp_shift_proxy_min"] <= result["wall_cp_shift_proxy_max"]
