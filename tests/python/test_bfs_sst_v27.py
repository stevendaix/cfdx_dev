import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/"scripts"))
from diagnose_bfs_sst_v27 import first_attempt


def trace(*, accepted="1", selected="sst_subspace_trust_v27", linf="0.4",
          l2="0.3", radius="0.5", ratio="0.3", slope="-0.2", qualified="0"):
    return ("SST_JFNK_V22 iteration=10 v26_active=1 v27_active=1 "
            "v27_bases=4 v27_evaluations=30 v27_valid=16 "
            f"v27_accepted={accepted} v27_linf={linf} v27_l2={l2} "
            f"v27_radius={radius} v27_ratio={ratio} v27_slope={slope} "
            "v27_choice=projected_l2_cauchy v26_l2_before=0.6 "
            f"old_merit=0.9 new_merit=0.4 selected={selected} "
            f"qualified={qualified}\n"
            "CFDX nonlinear iteration 1 rejected; retry=1\n"
            "SST_JFNK_V22 iteration=10 v27_active=1 qualified=1")


def test_v27_reports_first_attempt_only_and_not_bfs_convergence():
    v=first_attempt(trace())
    assert v["v27_accepted"] and v["v27_selected"]
    assert v["v27_bases"] == 4 and v["v27_evaluations"] == 30
    assert not v["sst_converged"] and not v["bfs_converged"]


def test_v27_rejects_l2_only_improvement_and_bad_direction_ratio():
    with pytest.raises(ValueError): first_attempt(trace(linf="1.0"))
    with pytest.raises(ValueError): first_attempt(trace(l2="0.7"))
    with pytest.raises(ValueError): first_attempt(trace(ratio="-0.1"))
    with pytest.raises(ValueError): first_attempt(trace(slope="0.1"))


def test_v27_rejects_false_selection_and_invalid_trust_radius():
    with pytest.raises(ValueError): first_attempt(trace(accepted="0"))
    with pytest.raises(ValueError): first_attempt(trace(radius="3.0"))
    with pytest.raises(ValueError): first_attempt(trace(qualified="1"))


def test_v27_rejects_missing_first_attempt_and_nan_metrics():
    with pytest.raises(ValueError): first_attempt("SST_JFNK_V22 iteration=9 v27_active=1 qualified=0")
    with pytest.raises(ValueError): first_attempt(trace(ratio="nan"))
    with pytest.raises(ValueError): first_attempt(trace().replace("v27_bases=4", "v27_bases=14"))
