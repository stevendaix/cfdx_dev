"""Opt-in BFS pressure-seed import: data validation, not CFD convergence."""

from __future__ import annotations

import os
import shutil
import subprocess

import h5py
import numpy as np
import pytest

from scripts.generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
from scripts.initialize_th03_hybrid import initialise_case

SOLVER = os.environ.get("CFDX_HYBRID_PRODUCTION_BIN")
pytestmark = pytest.mark.skipif(
    not SOLVER, reason="production solver binary not supplied"
)


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    root = tmp_path_factory.mktemp("bfs-pressure-import")
    raw = root / "raw.cfdx.h5"
    hybrid = root / "hybrid.cfdx.h5"
    write_case(
        raw,
        make_structured_mesh("smoke"),
        mesh_source="pressure-import-test",
        max_iterations=1,
        sst_provisional=True,
        sst_relative_gate=True,
        sst_inner_iterations=10,
        sst_newton_sink=True,
        sst_omega_first=True,
        v18_sparse_flux_pressure_jacobian=True,
        v19_full_pressure_flux_step=True,
        pressure_tolerance=1e-10,
    )
    initialise_case(raw, hybrid)
    return root, hybrid


def run(root, case_path, extra=()):
    proc = subprocess.run(
        [
            str(SOLVER),
            "--mesh",
            str(case_path),
            "--output-dir",
            str(root / "output"),
            "--experimental-sst",
            "--initial-velocity-field",
            "--initial-sst-fields",
            "--initial-face-mass-flux",
            "--iterations",
            "1",
            *extra,
        ],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    return proc.stdout + proc.stderr


def variant(case, tmp_path, kind):
    root, original = case
    target = tmp_path / "seeded.cfdx.h5"
    shutil.copy2(original, target)
    with h5py.File(target, "r+") as file:
        scalar = file.require_group("fields").require_group("scalar")
        cells = int(file.attrs["n_cells"])
        if kind != "absent":
            data = np.full(cells + (1 if kind == "shape" else 0), 0.0)
            if kind == "nan":
                data[0] = np.nan
            scalar.create_dataset("p", data=data)
    return target


def test_valid_seed_is_imported_without_claiming_convergence(case, tmp_path):
    path = variant(case, tmp_path, "valid")
    log = run(tmp_path, path, ("--initial-pressure-field",))
    assert "INITIAL_PRESSURE_FIELD_IMPORTED cells=" in log
    assert "initial guess only; NOT a convergence certificate" in log


def test_missing_pressure_seed_fails_closed(case, tmp_path):
    log = run(
        tmp_path, variant(case, tmp_path, "absent"), ("--initial-pressure-field",)
    )
    assert "scalar field 'p' is absent" in log
    assert "INITIAL_PRESSURE_FIELD_IMPORTED" not in log


@pytest.mark.parametrize(
    "kind,error",
    [
        ("shape", "pressure field shape mismatch"),
        ("nan", "pressure field is not finite"),
    ],
)
def test_invalid_pressure_seed_is_rejected(case, tmp_path, kind, error):
    log = run(tmp_path, variant(case, tmp_path, kind), ("--initial-pressure-field",))
    assert error in log
    assert "INITIAL_PRESSURE_FIELD_IMPORTED" not in log


def test_pressure_seed_requires_velocity_import(case, tmp_path):
    path = variant(case, tmp_path, "valid")
    # Rejection happens before the pressure seed or any FV equations are used.
    proc = subprocess.run(
        [
            str(SOLVER),
            "--mesh",
            str(path),
            "--output-dir",
            str(tmp_path / "bad-output"),
            "--initial-pressure-field",
            "--experimental-sst",
            "--iterations",
            "1",
        ],
        capture_output=True,
        text=True,
        timeout=60,
        check=False,
    )
    assert "--initial-pressure-field requires --initial-velocity-field" in (
        proc.stdout + proc.stderr
    )
