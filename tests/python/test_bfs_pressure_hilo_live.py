"""Rejection and pairing contracts for live production high/low diagnostics."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts/qualify_bfs_pressure_hilo_live.py"
)
spec = importlib.util.spec_from_file_location("qualify_bfs_pressure_hilo_live", SCRIPT)
assert spec and spec.loader
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


def _line(prefix, *, face):
    if face:
        return (
            f"{prefix} iter=1 corrector=0 before_relative=2e-10 "
            "after_relative=1e-13 original_csr_after_relative=2e-10 "
            "full_rhie_chow_qualified=0"
        )
    return (
        f"{prefix} iter=1 corrector=0 original_double_relative=4e-10 "
        "final_split_relative=1e-18 final_fp64_relative=2e-10 "
        "physical_fv_qualified=0"
    )


def test_live_report_fails_closed_on_different_operator():
    log = (
        _line("PRESSURE_HILO_FACES_LIVE", face=True)
        + "\n"
        + _line("PRESSURE_HILO_LIVE", face=False)
    )
    report = module.parse_log(log)
    assert report["status"] == "DIAGNOSTIC_ONLY_NOT_QUALIFIED"
    assert report["comparisons"][0]["csr_split_pass"]
    assert report["comparisons"][0]["face_implicit_pass"]
    assert not report["comparisons"][0]["both_face_and_csr_after_face_refinement_pass"]
    assert not report["physical_FV_Rhie_Chow_qualified"]


def test_missing_or_mismatched_state_is_rejected():
    good = _line("PRESSURE_HILO_LIVE", face=False)
    with pytest.raises(ValueError, match="missing or unpaired"):
        module.parse_log(good)
    bad = _line("PRESSURE_HILO_FACES_LIVE", face=True).replace("iter=1", "iter=2")
    with pytest.raises(ValueError, match="different FV states"):
        module.parse_log(good + "\n" + bad)


def test_nonfinite_and_fake_physical_qualification_rejected():
    face = _line("PRESSURE_HILO_FACES_LIVE", face=True)
    csr = _line("PRESSURE_HILO_LIVE", face=False)
    with pytest.raises(ValueError, match="nonfinite"):
        module.parse_log(
            face.replace("after_relative=1e-13", "after_relative=nan") + "\n" + csr
        )
    with pytest.raises(ValueError, match="unexpected physical"):
        module.parse_log(
            face
            + "\n"
            + csr.replace("physical_fv_qualified=0", "physical_fv_qualified=1")
        )
