"""Fails closed on absent/mixed/incomplete Fp capture replays."""
from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location(
    "pcd_fp_replay_campaign", ROOT / "scripts/pcd_fp_replay_campaign.py"
)
assert SPEC and SPEC.loader
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def row(method, capture, residual=1e-9):
    return dict(method=method, capture=str(capture), rows=3, nnz=7,
                setup_ok=True, first_status=0, first_true_relative_rhs=residual,
                best_true_relative_rhs=residual, second_attempted=False,
                second_improved=False, second_status=0)


def test_missing_captures_rejected(tmp_path):
    with pytest.raises(ValueError, match="missing"):
        mod.compare([], tmp_path)


def test_captured_methods_same_input_and_strict_pass(tmp_path):
    capture = tmp_path / "model.cfp"
    capture.write_bytes(b"CFDXFP01 source data")
    exe = tmp_path / "replay"
    exe.write_text("fixed fixture")
    records = [row(name, capture) for name in mod.METHODS]

    def runner(_args, **_kwargs):
        return SimpleNamespace(returncode=0, stdout="\n".join(map(json.dumps, records)), stderr="")

    report = mod.compare([capture], exe, runner=runner)
    assert len(report["records"]) == len(mod.METHODS)
    assert all(rec["two_pass_meets_true_tolerance"] for rec in report["records"])
    assert all(rec["capture_sha256"] == mod.digest(capture) for rec in report["records"])


def test_incomplete_methods_fail_closed(tmp_path):
    f = tmp_path / "input.cfp"
    f.write_bytes(b"CFDXFP01")
    exe = tmp_path / "replay"
    exe.write_bytes(b"fixture")

    def runner(_args, **_kwargs):
        return SimpleNamespace(returncode=0, stdout=json.dumps(row("ilu0", f)), stderr="")

    with pytest.raises(ValueError, match="incomplete"):
        mod.compare([f], exe, runner=runner)


def test_nonconverged_never_promoted(tmp_path):
    f = tmp_path / "input.cfp"
    f.write_bytes(b"CFDXFP01")
    exe = tmp_path / "replay"
    exe.write_bytes(b"fixture")
    records = [row(method, f, residual=1e-3) for method in mod.METHODS]

    def runner(_args, **_kwargs):
        return SimpleNamespace(returncode=0, stdout="\n".join(map(json.dumps, records)), stderr="")

    report = mod.compare([f], exe, runner=runner)
    assert not any(r["two_pass_meets_true_tolerance"] for r in report["records"])
