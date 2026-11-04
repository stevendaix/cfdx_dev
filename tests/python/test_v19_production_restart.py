"""Opt-in V19 production integration: replay the authoritative face flux.

Run with CFDX_PRODUCTION_SOLVER_BINARY=/absolute/path/to/cfdx_production_solver
so the test uses the actual compiled production executable. Normal unit-test
jobs without a production binary skip, instead of accepting a fake solver.
"""
import json
import os
import shutil
from pathlib import Path
import subprocess

import h5py
import pytest

from scripts.generate_th03_openfoam_blockmesh import (
    make_structured_mesh,
    write_case,
)


def _solver_binary() -> Path:
    value = os.environ.get("CFDX_PRODUCTION_SOLVER_BINARY")
    if not value:
        pytest.skip("set CFDX_PRODUCTION_SOLVER_BINARY for the full V19 restart test")
    path = Path(value)
    if not path.is_file():
        pytest.fail(f"production binary is missing: {path}")
    return path


def _run(solver: Path, case: Path, output: Path, iterations: int,
         restart: Path | None = None) -> subprocess.CompletedProcess[str]:
    cmd = [str(solver), "--mesh", str(case), "--output-dir", str(output),
           "--iterations", str(iterations)]
    if restart is not None:
        cmd.extend(["--restart", str(restart)])
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120,
                          check=False)
    assert proc.returncode in (0, 1), (
        f"production process failed: {proc.returncode}\n"
        f"stdout tail: {proc.stdout[-1500:]}\n"
        f"stderr tail: {proc.stderr[-1500:]}")
    assert (output / "restart.dat").is_file()
    return proc


def test_v19_split_pressure_flux_restart_replays_uninterrupted_mdot(tmp_path: Path):
    solver = _solver_binary()
    case = tmp_path / "v19_smoke.cfdx.h5"
    write_case(case, make_structured_mesh("smoke"),
               mesh_source="V19 replay contract fixture", max_iterations=4,
               alpha_u=0.25, alpha_p=0.1, pressure_tolerance=1e-10,
               v18_sparse_flux_pressure_jacobian=True,
               v19_full_pressure_flux_step=True)
    # Run a diagnostic LAMINAR carrier explicitly: the fixture generator's
    # RANS SST declaration is intentionally rejected without --experimental-sst.
    # This test isolates the pressure/flux state and must not fake SST fields.
    with h5py.File(case, "r+") as handle:
        setup = json.loads(handle.attrs["case_setup_json"])
        setup["physics_model"] = "incompressible"
        setup["turbulence_model"] = "laminar"
        setup["materials"][0]["turbulence_model"] = "laminar"
        handle.attrs["case_setup_json"] = json.dumps(setup)
    _run(solver, case, tmp_path / "full", 4)
    _run(solver, case, tmp_path / "split", 2)
    checkpoint = tmp_path / "split" / "restart.dat"
    _run(solver, case, tmp_path / "resumed", 2, checkpoint)

    whole = (tmp_path / "full" / "restart.dat").read_text().splitlines()
    resumed = (tmp_path / "resumed" / "restart.dat").read_text().splitlines()
    # Current DAT iteration metadata is per solver invocation. Field and
    # face-flux contents must still be byte-identical to a continuous run.
    assert whole[2] == "iteration 4"
    assert resumed[2] == "iteration 2"
    assert whole[:2] == resumed[:2]
    assert whole[3:] == resumed[3:]
    assert any(line.startswith("face_fields ") for line in whole)
    assert "field mdot 1" in whole

    # A valid physical case named like a plain mesh must not silently lose
    # its inlet/wall pressure-velocity boundary conditions and converge at 0.
    mislabeled = tmp_path / "case_with_wrong_suffix.h5"
    shutil.copyfile(case, mislabeled)
    wrong_name = subprocess.run(
        [str(solver), "--mesh", str(mislabeled), "--output-dir",
         str(tmp_path / "wrong_case_name"), "--iterations", "2"],
        capture_output=True, text=True, timeout=60, check=False)
    assert wrong_name.returncode == 2
    assert "contains case_setup_json" in (wrong_name.stdout + wrong_name.stderr)

    # Old/incomplete checkpoints are rejected: recomputing mdot from the
    # *relaxed* stored pressure would silently change the next iteration.
    truncated = tmp_path / "without_mdot.dat"
    original = checkpoint.read_text()
    assert "face_fields " in original
    truncated.write_text(original.split("face_fields ", 1)[0])
    failed = subprocess.run(
        [str(solver), "--mesh", str(case), "--output-dir",
         str(tmp_path / "missing_flux"), "--iterations", "1",
         "--restart", str(truncated)],
        capture_output=True, text=True, timeout=60, check=False)
    assert failed.returncode == 2
    assert "requested mdot but DAT has none" in (failed.stdout + failed.stderr)