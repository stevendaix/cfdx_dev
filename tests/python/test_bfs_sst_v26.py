import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from diagnose_bfs_sst_v26 import first_attempt


def trace(*, selected="jfnk", l2="0.3", linf="0.7", smooth="1", mismatch="0.08"):
    return ("SST_JFNK_V22 iteration=10 v26_active=1 v26_jvp_available=1 "
            "v26_jvp_evaluations=7 v26_jvp_h1=0.001 v26_jvp_h2=0.0001 "
            f"v26_jvp_h3=0.00001 v26_jvp_max_disagreement={mismatch} "
            f"v26_l2_before=0.6 v26_l2_best_trial=0.25 v26_l2_accepted={smooth} "
            f"v26_l2_selected={l2} v26_linf_selected={linf} v26_alpha=0.5 "
            f"old_merit=0.9 new_merit=0.6 selected={selected} qualified=0\n"
            "CFDX nonlinear iteration 1 rejected; retry=1\n"
            "SST_JFNK_V22 iteration=10 v26_active=1 v26_jvp_available=1 qualified=1")


def test_v26_reports_only_first_attempt_without_certifying_bfs():
    report = first_attempt(trace())
    assert report["maximum_directional_difference"] == pytest.approx(.08)
    assert report["l2_proposal_physically_admissible"]
    assert not report["sst_converged"]
    assert not report["bfs_converged"]
    assert not report["nasa_validated"]


def test_v26_refuses_l2_only_candidate():
    with pytest.raises(ValueError):
        first_attempt(trace(linf="0.95"))
    with pytest.raises(ValueError):
        first_attempt(trace(l2="0.65"))


def test_v26_rejects_invalid_or_unqualified_directional_audit():
    with pytest.raises(ValueError):
        first_attempt(trace(mismatch="nan"))
    with pytest.raises(ValueError):
        first_attempt(trace().replace("v26_jvp_h2=0.0001", "v26_jvp_h2=0.01"))
    with pytest.raises(ValueError):
        first_attempt(trace().replace("v26_jvp_evaluations=7", "v26_jvp_evaluations=3"))


def test_v26_rejects_false_smooth_selection_and_missing_evidence():
    with pytest.raises(ValueError):
        first_attempt(trace(smooth="0", selected="jfnk_smooth_v26"))
    with pytest.raises(ValueError):
        first_attempt("SST_JFNK_V22 iteration=10 old_merit=1 qualified=0")
    with pytest.raises(ValueError):
        first_attempt(trace().replace("qualified=0", "qualified=1"))
