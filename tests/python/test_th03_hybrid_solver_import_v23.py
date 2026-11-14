"""V23 hybrid runtime integration and fail-closed dataset validation."""

from __future__ import annotations

import os
import shutil
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

from scripts.generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
from scripts.initialize_th03_hybrid import initialise_case

BIN = os.environ.get("CFDX_HYBRID_PRODUCTION_BIN")
pytestmark = pytest.mark.skipif(
    not BIN, reason="set CFDX_HYBRID_PRODUCTION_BIN to compiled solver"
)


@pytest.fixture(scope="module")
def source(tmp_path_factory):
    d = tmp_path_factory.mktemp("hybrid-import")
    raw = d / "raw.cfdx.h5"
    hybrid = d / "hybrid.cfdx.h5"
    write_case(
        raw,
        make_structured_mesh("smoke"),
        mesh_source="v23-test",
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
    initialise_case(raw, hybrid, require_native_geometry=False)
    return d, hybrid


def run(binary, case: Path, out: Path, extra=()):
    p = subprocess.run(
        [
            binary,
            "--mesh",
            str(case),
            "--output-dir",
            str(out),
            "--experimental-sst",
            "--initial-velocity-field",
            "--initial-sst-fields",
            "--initial-face-mass-flux",
            "--iterations",
            "1",
            *extra,
        ],
        check=False,
        text=True,
        capture_output=True,
        timeout=45,
    )
    return p, p.stdout + p.stderr


def test_hybrid_runtime_reads_all_fields_and_accepts_projection(source):
    d, hybrid = source
    _p, logs = run(BIN, hybrid, d / "good")
    assert "HYBRID_INITIAL_SST_FIELDS imported" in logs
    assert "HYBRID_INITIAL_FACE_FLUX continuity_l1=" in logs
    assert "INCOMPRESSIBLE_ITER iter=1" in logs
    # The solver returns nonzero on an unfinished run; this is expected.
    assert "Converged NO Iterations 1" in logs
    assert "fatal" not in logs.lower()


def test_hybrid_rejects_invalid_sst_scalar(source):
    d, hybrid = source
    bad = d / "badomega.cfdx.h5"
    shutil.copy2(hybrid, bad)
    with h5py.File(bad, "r+") as f:
        f["fields/scalar/omega"][0] = -1.0
    _p, logs = run(BIN, bad, d / "badomega_out")
    assert "hybrid SST field contains an invalid or subfloor value" in logs
    assert "INCOMPRESSIBLE_ITER" not in logs


def test_hybrid_rejects_nonconservative_or_bc_flux(source):
    d, hybrid = source
    bad = d / "badflux.cfdx.h5"
    shutil.copy2(hybrid, bad)
    with h5py.File(bad, "r+") as f:
        f["fields/face/mdot"][50] += 0.005
    _p, logs = run(BIN, bad, d / "badflux_out")
    assert (
        "hybrid initial mass flux fails local continuity gate" in logs
        or "hybrid flux violates" in logs
    )
    assert "INCOMPRESSIBLE_ITER" not in logs


def test_hybrid_rejects_bad_sst_shape_before_hdf5_read(source):
    d, hybrid = source
    bad = d / "badshape.cfdx.h5"
    shutil.copy2(hybrid, bad)
    with h5py.File(bad, "r+") as f:
        data = f["fields/scalar/omega"][:]
        del f["fields/scalar/omega"]
        f["fields/scalar"].create_dataset("omega", data=np.column_stack((data, data)))
    _p, logs = run(BIN, bad, d / "badshape_out")
    assert "hybrid SST fields cannot be read from HDF5" in logs
    assert "INCOMPRESSIBLE_ITER" not in logs
