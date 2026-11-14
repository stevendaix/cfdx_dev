from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
MOD = ROOT / "scripts" / "diagnose_bfs_sst_v22.py"
SPEC = importlib.util.spec_from_file_location("diagnose_bfs_sst_v22", MOD)
assert SPEC and SPEC.loader
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def _line(label, iteration=1):
    if label == "SST_TRANSPORT_DIAGNOSTIC":
        return (f"{label} iteration={iteration} k_returned_nonlinear_raw=0.01 "
                "omega_returned_nonlinear_raw=1.0 k_returned_row_scaled=2.0 "
                "omega_returned_row_scaled=3.0 k_local_relative=0.2 "
                "omega_local_relative=0.4")
    return (f"{label} iteration={iteration} finite=1 outer_update_ok=0 "
            "returned_equations_ok=0")


def test_v22_diagnostic_exposes_failure_not_certification():
    log = "\n".join([_line("SST_TRANSPORT_DIAGNOSTIC"),
                     "SST_JFNK_V22 iteration=1 accepted=1 selected=picard",
                     _line("SST_ACCEPTANCE_GATES"),
                     "Converged NO Iterations 1"])
    report = module.parse_run(log)
    assert not report["qualified_physics"]
    assert not report["outer_converged_as_reported"]
    assert not report["runs"][0]["all_reported_gates_closed"]
    assert report["runs"][0]["samples"][0]["omega_returned_nonlinear_raw"] == 1.0
    assert report["runs"][0]["accepted_steps_in_jfnk_mode"] == 1
    assert report["runs"][0]["newton_steps_selected"] == 0
    assert report["runs"][0]["picard_steps_selected"] == 1


def test_v22_diagnostic_requires_real_samples():
    with pytest.raises(ValueError, match="no complete SST"):
        module.parse_run("SST_JFNK_V22 iteration=1 accepted=1")


def test_v22_diagnostic_rejects_nan_and_negative():
    with pytest.raises(ValueError, match="nonfinite"):
        module.parse_run(_line("SST_TRANSPORT_DIAGNOSTIC").replace(
            "omega_returned_row_scaled=3.0", "omega_returned_row_scaled=NaN"))
    with pytest.raises(ValueError, match="negative"):
        module.parse_run(_line("SST_TRANSPORT_DIAGNOSTIC").replace(
            "omega_returned_row_scaled=3.0", "omega_returned_row_scaled=-2"))


def test_v22_trace_multi_attempts_retains_each():
    log = "\n".join([_line("SST_TRANSPORT_DIAGNOSTIC", 1),
                     _line("SST_TRANSPORT_DIAGNOSTIC", 10),
                     _line("SST_TRANSPORT_DIAGNOSTIC", 1)])
    report = module.parse_run(log)
    assert len(report["runs"]) == 2
    assert len(report["runs"][0]["samples"]) == 2


def test_v22_trace_distinguishes_picard_from_actual_newton_step():
    log = "\n".join([_line("SST_TRANSPORT_DIAGNOSTIC", 1),
                     "SST_JFNK_V22 iteration=1 accepted=1 proposed=0 selected=picard",
                     "SST_JFNK_V22 iteration=10 accepted=1 proposed=1 selected=jfnk",
                     _line("SST_TRANSPORT_DIAGNOSTIC", 10)])
    run = module.parse_run(log)["runs"][0]
    assert run["picard_steps_selected"] == 1
    assert run["newton_steps_selected"] == 1
    assert run["accepted_steps_unknown_family"] == 0
