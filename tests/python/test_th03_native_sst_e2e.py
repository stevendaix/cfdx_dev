"""Run actual CFDX SST on a tiny BFS; software regression, not TH-03 validation.

The synthetic inlet and low-resolution mesh can yield a perfectly converged
but physically wrong reattachment. This test must never certify TH-03.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import audit_th03_wall_resolution as wall_audit
from scripts import compare_th03_wall_gradients as diag
from scripts import generate_th03_bfs_case as bfs


def test_debug_bfs_pressure_outlet_matches_physical_operator() -> None:
    setup = json.loads(bfs.make_setup(height=0.0127, reynolds=36000, velocity=44.2))
    outlet = next(
        boundary
        for boundary in setup["boundary_conditions"]
        if boundary["patch_name"] == "outlet"
    )
    assert outlet["type"] == "pressure_outlet"
    assert outlet["value_type"] == "outlet_pressure"
    assert outlet["pressure"] == 0


@pytest.mark.skipif(
    os.environ.get("CFDX_RUN_TH03_DEBUG_E2E") != "1",
    reason="opt-in actual SST BFS solve (145 nonlinear iterations)",
)
def test_real_sst_bfs_debug_is_never_physically_qualified(tmp_path: Path) -> None:
    solver = Path(
        os.environ.get(
            "CFDX_PRODUCTION_SOLVER",
            str(ROOT / "build/th-preflight/cfdx_production_solver"),
        )
    )
    assert solver.is_file(), f"real CFDX solver unavailable: {solver}"
    case = tmp_path / "debug_bfs.cfdx.h5"
    out = tmp_path / "out"
    metadata = bfs.write_case(case)
    assert metadata["physical_qualification"] == "NOT_QUALIFIED"
    args = [
        str(solver),
        "--mesh",
        str(case),
        "--output-dir",
        str(out),
        "--iterations",
        "250",
        "--experimental-sst",
        "--benchmark-export",
    ]
    with (tmp_path / "solver.log").open("w") as log:
        run = subprocess.run(
            args, stdout=log, stderr=subprocess.STDOUT, timeout=240, check=False
        )
    solver_log = (tmp_path / "solver.log").read_text()
    assert run.returncode in (0, 1), solver_log[-5000:]
    execution_path = out / "execution.json"
    assert execution_path.is_file(), (
        "production solver did not write execution.json; "
        f"returncode={run.returncode}; log tail:\\n{solver_log[-5000:]}"
    )
    execution = json.loads(execution_path.read_text())
    convergence_path = out / "convergence.json"
    assert convergence_path.is_file(), (
        "production solver did not write convergence.json; "
        f"returncode={run.returncode}; execution={execution}; "
        f"log tail:\\n{solver_log[-5000:]}"
    )
    flow = json.loads(convergence_path.read_text())
    assert execution["converged"] is (run.returncode == 0)
    assert flow["converged"] is (run.returncode == 0)
    face = out / "benchmark_observables/face_fluxes.csv"
    cell = out / "benchmark_observables/cell_fields.csv"
    wall = wall_audit.audit(
        face,
        cell,
        wall_patches=("downstream_wall",),
        density=1.2,
        viscosity=44.2 * 0.0127 / 36000,
    )
    assert wall["status"] == "BLOCKED" and wall["max_y_plus_estimate"] > 1.0
    assert wall["physical_qualification"] == "NOT_QUALIFIED"
    # The coarse-grid reattachment is deliberately not compared with the
    # experiment: the cell wall resolution and inlet profile are unacceptable.
    if run.returncode == 0:
        mass = json.loads(
            (out / "benchmark_observables/mass_conservation.json").read_text()
        )
        assert mass["global_mass_imbalance_relative"] < 1e-7
        result = diag.compare(
            face,
            cell,
            patch="downstream_wall",
            density=1.2,
            viscosity=44.2 * 0.0127 / 36000,
            bulk_speed=44.2,
            step_height=0.0127,
        )
        assert result["status"] == "DIAGNOSTIC_ONLY"
