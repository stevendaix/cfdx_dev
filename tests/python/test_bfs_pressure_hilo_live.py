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


def _bridge_line(*, iterate=1, relative="2e-11", qualified=0):
    return (
        f"PRESSURE_HILO_RHIE_CHOW_BRIDGE iter={iterate} corrector=0 "
        "rhs_base_flux_relative=1e-13 "
        "implicit_linearized_continuity_relative=1e-13 "
        "implicit_operator_relative=1e-13 "
        f"full_rc_compensated_continuity_relative={relative} "
        "full_rc_fp64_continuity_relative=9e-10 "
        "full_rc_face_rounding_l2_over_rhs=2e-10 "
        "rhs_base_contract_pass=1 "
        f"full_rc_qualified={qualified} pressure_state_unchanged=1"
    )


def test_new_full_flux_bridge_is_diagnostic_and_keeps_separate_gates():
    log = "\n".join(
        (
            _line("PRESSURE_HILO_FACES_LIVE", face=True),
            _line("PRESSURE_HILO_LIVE", face=False),
            _bridge_line(),
        )
    )
    report = module.parse_log(log)
    bridge = report["comparisons"][0]["rhie_chow_bridge"]
    assert report["n_rhie_chow_bridge_diagnostics"] == 1
    assert bridge["available"]
    assert bridge["rhs_base_contract_pass"]
    assert bridge["full_rc_compensated_continuity_relative"] > report["threshold"]
    assert report["status"] == "DIAGNOSTIC_ONLY_NOT_QUALIFIED"
    assert not report["complete_flux_gate"]


def test_bridge_rejects_invalid_or_unpaired_results():
    good = "\n".join(
        (
            _line("PRESSURE_HILO_FACES_LIVE", face=True),
            _line("PRESSURE_HILO_LIVE", face=False),
        )
    )
    with pytest.raises(ValueError, match="nonfinite"):
        module.parse_log(good + "\n" + _bridge_line(relative="nan"))
    with pytest.raises(ValueError, match="different FV state"):
        module.parse_log(good + "\n" + _bridge_line(iterate=2))
    with pytest.raises(ValueError, match="cannot qualify"):
        module.parse_log(good + "\n" + _bridge_line(qualified=1))


def test_bridge_deferred_spatial_location_is_preserved_and_validated():
    old = "\n".join(
        (
            _line("PRESSURE_HILO_FACES_LIVE", face=True),
            _line("PRESSURE_HILO_LIVE", face=False),
        )
    )
    spatial = (
        " deferred_divergence_relative=0.72 deferred_face_l2_over_rhs=1.1"
        " worst_deferred_face=12 worst_deferred_cell=5"
        " worst_full_divergence_cell=6"
        " worst_full_divergence_x=-0.03"
        " worst_full_divergence_y=0.011"
        " worst_full_divergence_z=0.004"
    )
    report = module.parse_log(old + "\n" + _bridge_line() + spatial)
    bridge = report["comparisons"][0]["rhie_chow_bridge"]
    assert bridge["deferred_divergence_relative"] == pytest.approx(0.72)
    assert bridge["worst_deferred_face"] == 12
    assert bridge["worst_full_divergence_x"] == pytest.approx(-0.03)
    with pytest.raises(ValueError, match="nonfinite/missing spatial"):
        module.parse_log(
            old
            + "\n"
            + _bridge_line()
            + spatial.replace(
                "worst_full_divergence_y=0.011", "worst_full_divergence_y=nan"
            )
        )
