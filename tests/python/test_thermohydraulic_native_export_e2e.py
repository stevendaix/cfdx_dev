"""Real production-executable benchmark-output smoke, NOT a TH validation.

Uses a tiny, zero-body-force quiescent channel so the expected mass flux is
exactly zero. No claimed DNS/RANS accuracy; exercises the complete HDF5 input,
final solver state, native CSVs, and velocity-bearing VTU path.
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
sys.path.insert(0, str(ROOT / "tests/python"))
import test_n9_poiseuille_body_force_production as channel  # noqa: E402


def test_real_production_solver_thermohydraulic_export(tmp_path: Path) -> None:
    solver_name = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver_name:
        pytest.skip("CFDX_PRODUCTION_SOLVER must point to a built production binary")
    solver = Path(solver_name)
    assert solver.is_file()

    # Deliberately a numerical smoke, not the original forced Poiseuille
    # validation: remove the driving force and require exactly quiescent flow.
    channel.NX = 3
    channel.NY = 4
    channel.MAX_ITERATIONS = 5
    case_file = tmp_path / "native_thermo_smoke.cfdx.h5"
    channel._write_case(case_file)
    with h5py.File(case_file, "r+") as h5:
        setup = json.loads(h5.attrs["case_setup_json"])
        setup["body_force"] = [0.0, 0.0, 0.0]
        setup["numerics"]["under_relaxation_momentum"] = 0.2
        setup["numerics"]["under_relaxation_pressure"] = 0.1
        h5.attrs["case_setup_json"] = json.dumps(setup)

    output = tmp_path / "solver_output"
    result = subprocess.run([str(solver), "--mesh", str(case_file),
                             "--output-dir", str(output), "--iterations", "5",
                             "--benchmark-export"],
                            capture_output=True, text=True, timeout=90, check=False)
    assert result.returncode == 0, (result.stdout[-3000:], result.stderr[-3000:])
    assert 'Case under-relaxation alpha_u=0.2 alpha_p=0.1' in result.stdout
    assert 'alpha_u=0.2 alpha_p=0.1' in result.stdout
    summary = json.loads((output / "execution.json").read_text(encoding="utf-8"))
    mass = json.loads((output / "benchmark_observables/mass_conservation.json").read_text(encoding="utf-8"))
    assert summary["converged"] is True
    assert summary["process_exit_code"] == 0
    assert mass["cells"] == 12
    assert abs(mass["boundary_net_mass_flux_kg_s"]) < 1e-12
    assert mass["global_mass_imbalance_relative"] == pytest.approx(0.0, abs=1e-12)
    assert mass["maximum_cell_mass_imbalance_relative"] == pytest.approx(0.0, abs=1e-12)
    assert mass["physical_qualification"] == "NOT_QUALIFIED"
    with (output / "benchmark_observables/cell_fields.csv").open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 12
    for row in rows:
        for key in ("u_x_m_s", "u_y_m_s", "u_z_m_s", "p_static_Pa"):
            assert abs(float(row[key])) < 1e-10
    with (output / "benchmark_observables/face_fluxes.csv").open(newline="", encoding="utf-8") as stream:
        faces = list(csv.DictReader(stream))
    assert len(faces) == mass["faces"]
    assert all(abs(float(row["mass_flux_kg_s"])) < 1e-10 for row in faces)
    vtk = sorted(output.glob("result_*.vtu"))
    assert vtk
    text = vtk[-1].read_text(encoding="utf-8")
    for name in ('Name="p"', 'Name="u_x"', 'Name="u_y"', 'Name="u_z"'):
        assert name in text
