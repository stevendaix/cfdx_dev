"""CI integration gate: real native CFDX solver -> real trace -> reader."""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import pytest
from cfdx.cell_convergence import read_cell_trace


def test_production_solver_writes_accepted_cell_snapshots(tmp_path):
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    mesh = os.environ.get("CFDX_PRODUCTION_MESH")
    if not solver or not mesh:
        pytest.skip("native solver and mesh supplied by the CFDX CTest environment")
    out = tmp_path / "native_cell_diagnostics"
    result = subprocess.run(
        [
            solver,
            "--mesh",
            mesh,
            "--output-dir",
            str(out),
            "--iterations",
            "5",
            "--cell-hotspots",
            "--cell-hotspots-every",
            "2",
        ],
        check=False,
        text=True,
        capture_output=True,
        timeout=60,
    )
    assert result.returncode in (0, 1), result.stderr[-4000:]
    trace = read_cell_trace(out / "cell_convergence.jsonl")
    assert trace.n_cells > 0
    assert trace.n_faces > 0
    assert trace.frames, "at least one native nonlinear iteration must be accepted"
    assert trace.complete
    assert trace.accepted_iterations >= trace.latest.iteration
    assert not trace.to_dict()["solver_convergence_certificate"]
    assert Path(trace.mesh_source).resolve() == Path(mesh).resolve()
    assert all(
        sample.cell_id < trace.n_cells
        for frame in trace.frames
        for sample in frame.cells
    )
