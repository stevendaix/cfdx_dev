"""The full CFDX NACA/N9 campaign must not silently shorten/relabel runs."""
from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "coupled_matrix_campaign", ROOT / "scripts/run_cfdx_coupled_campaign.py"
)
assert spec and spec.loader
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)


def args(tmp_path, mode="naca"):
    return SimpleNamespace(
        build=tmp_path / "build", meshes=tmp_path / "meshes", out=tmp_path / "out",
        only=mode, models=["pcd", "block_local"], n9_models=["pcd", "lsc"],
        levels=[64, 128, 256], iterations=1500, wall_timeout=None,
    )


def test_three_grids_six_commands_same_final_scheme(tmp_path):
    items = list(campaign.build_commands(args(tmp_path)))
    assert len(items) == 6
    for tag, cmd, env in items:
        assert cmd[-1] == "1500"
        assert env["CFDX_NACA_ALGORITHM"] == "COUPLED"
        assert env["CFDX_NACA_LINEAR_TOLERANCE"] == "1e-12"
        assert env["CFDX_NACA_CONVECTION_CELL_LIMITED"] == "1"
        assert env["CFDX_NACA_FARFIELD_PRESSURE_NEUMANN"] == "1"
        assert env["CFDX_NACA_WALL_TRACTION_FACE"] == "1"
        assert env["CFDX_NACA_WALL_GRADIENT_WLS"] == "1"
        if tag.endswith("pcd"):
            assert "CFDX_PCD_FP_CAPTURE_PREFIX" in env
        else:
            assert "CFDX_PCD_FP_CAPTURE_PREFIX" not in env


def test_n9_disjoint_runs(tmp_path):
    items = list(campaign.build_commands(args(tmp_path, "n9")))
    assert [item[0] for item in items] == ["n9_pcd", "n9_lsc"]
    assert all(item[1][1] == "--schur-model" for item in items)


def test_missing_inputs_not_promoted(tmp_path):
    result = campaign.execute(
        "naca_n64_pcd", [str(tmp_path / "not_built"), "--single",
                           str(tmp_path / "missing.h5"), "0.6", "0.25", "1500"],
        {}, tmp_path / "logs", wall_timeout=None,
    )
    assert result["execution"] == "NOT_RUN"
    assert result["physical_certificate"] == "NOT_QUALIFIED"
