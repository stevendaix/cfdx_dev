"""The generated BFS geometry and SST field BCs must survive Python -> HDF5."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import generate_th03_bfs_case as bfs
from scripts import run_thermohydraulic_cases as campaign


def test_bfs_mesh_has_complete_sst_scalar_boundary_contract(tmp_path: Path) -> None:
    path = tmp_path / "bfs.cfdx.h5"
    result = bfs.write_case(path, n_up=4, n_down=12, ny=10, step_cells=2)
    assert result["cells"] == (4 + 12) * 10 - 4 * 2
    assert result["physical_qualification"] == "NOT_QUALIFIED"
    with h5py.File(path) as h5:
        assert int(h5.attrs["n_cells"]) == result["cells"]
        setup = json.loads(h5.attrs["case_setup_json"])
        patches = {p["patch_name"]: p for p in setup["boundary_conditions"]}
        for name, patch in patches.items():
            scalars = {s["field"]: s["type"] for s in patch["sst_scalar_conditions"]}
            assert set(scalars) == {"k", "omega"}, name
            if name == "inlet":
                assert scalars == {"k": "fixed_value", "omega": "fixed_value"}
            elif name in ("outlet", "front", "back"):
                assert scalars == {"k": "zero_gradient", "omega": "zero_gradient"}
            else:
                assert scalars == {"k": "wall_k_zero", "omega": "wall_omega_sst"}


def test_benchmark_input_does_not_silently_accept_debug_bfs(tmp_path: Path) -> None:
    path = tmp_path / "bfs.cfdx.h5"
    result = bfs.write_case(path)
    with pytest.raises(
        campaign.CampaignError, match="benchmark_case/scenario provenance"
    ):
        campaign.validate_native_case_physics(
            path, "TH-03", "Re_h_37500", "RANS_SST", result["cells"]
        )
    # Falsifying labels in JSON cannot establish physical validation either;
    # only the preflight input-contract is tested here, never a PASS.


def test_mesh_generation_rejects_nonmatching_step_height() -> None:
    with pytest.raises(ValueError, match="step resolution"):
        bfs.make_mesh(n_up=4, n_down=12, ny=10, step_cells=5)


def test_tanh_wall_mesh_refines_both_sides_of_step_without_shifting_step() -> None:
    import numpy as np

    uniform = bfs.make_mesh(4, 16, 40, 16, 0.0127)
    refined = bfs.make_mesh(4, 16, 40, 16, 0.0127, wall_stretch=3.6)
    y0 = np.unique(uniform["points"][:, 1])
    y1 = np.unique(refined["points"][:, 1])
    assert len(y0) == len(y1) == 41
    assert any(abs(y - 0.0127) < 1e-12 for y in y1)
    assert y1[1] < y0[1] / 4
    assert y1[17] - y1[16] < (y0[17] - y0[16]) / 4
    assert np.min(np.diff(y1)) > 0
    assert refined["n_cells"] == uniform["n_cells"]


def test_wall_mesh_rejects_inadmissible_grading() -> None:
    import pytest

    with pytest.raises(ValueError, match="wall stretching"):
        bfs.make_mesh(4, 16, 10, 2, wall_stretch=4)
    with pytest.raises(ValueError, match="invalid"):
        bfs.make_mesh(4, 16, 40, 16, wall_stretch=100)
