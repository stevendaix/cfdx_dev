"""Ensure all six COUPLED Schur shards require actual independent evidence."""
from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path
from types import SimpleNamespace

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("coupled_campaign_strict", ROOT / "scripts/run_cfdx_coupled_campaign.py")
assert spec and spec.loader
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)


def options(tmp_path):
    return SimpleNamespace(build=tmp_path, meshes=tmp_path, out=tmp_path,
                           only="n9", n9_models=list(campaign.N9_MODELS))


def test_all_six_models_strict_by_default_and_no_inherited_naca(tmp_path, monkeypatch):
    monkeypatch.setenv("CFDX_NACA_CONVECTION_CELL_LIMITED", "0")
    monkeypatch.setenv("CFDX_PCD_FP_CAPTURE_PREFIX", "stale")
    runs = list(campaign.build_commands(options(tmp_path)))
    assert len(runs) == 6
    assert {tag for tag, _, _ in runs} == {"n9_" + model for model in campaign.SCHUR_MODELS}
    assert all(cmd[3] == "--require-independent-face" for _, cmd, _ in runs)
    assert all("CFDX_NACA_CONVECTION_CELL_LIMITED" not in env for _, _, env in runs)
    assert all(env["CFDX_PCD_FP_CAPTURE_PREFIX"] != "stale" for _, _, env in runs)


def test_exit_zero_with_missing_physical_stages_is_refused(tmp_path):
    command = [sys.executable, "-c", "print('STARTING COUPLED pcd')"]
    row = campaign.execute("n9_pcd", command, os.environ.copy(), tmp_path, wall_timeout=None)
    assert row["exit_code"] == 0
    assert row["execution"] == "EXIT_ZERO_BUT_N9_PHYSICAL_EVIDENCE_INVALID"
    assert row["physical_certificate"] == "NOT_QUALIFIED"
    assert not row["n9_independent_shard_gate"]["physical_shard_gate_complete"]


def test_four_synthetic_markers_only_validate_parser_not_physics(tmp_path):
    # The command intentionally only PRINTS records, not a real solve.
    # This proves only the parser's positive contract, never CFD accuracy.
    stages = ["Couette", "Poiseuille", "Ghia_Re100", "skew_Couette"]
    records = ["STARTING COUPLED block_local"] + [
        "COMPLETED COUPLED block_local " + stage + " " +
        "iterations=12 momentum_rhs_relative=1e-11 "
        "continuity_normalized=1e-12 mass_local_linf=1e-12 "
        "face_momentum_x_relative=1e-12 face_momentum_y_relative=1e-12"
        for stage in stages
    ]
    result = campaign.verify_n9_shard("\n".join(records), "block_local", 0)
    assert result["physical_shard_gate_complete"]
    assert result["completed_stages"] == stages
    duplicated = records.copy()
    duplicated[1] += " momentum_rhs_relative=1e-11"
    result = campaign.verify_n9_shard("\n".join(duplicated), "block_local", 0)
    assert not result["physical_shard_gate_complete"]
    assert any("duplicated" in error for error in result["issues"])
