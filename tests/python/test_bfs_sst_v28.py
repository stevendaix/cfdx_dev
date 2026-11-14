import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
from diagnose_bfs_sst_v28 import first_attempt

def trace(active="1", selected="jfnk", committed="0", linf="0.7",
          l2="0.4", radius_after="0.5", qualified="0", ratio="0.8",
          commit_source="none"):
    return ("SST_JFNK_V22 iteration=10 v27_active=1 v28_active=1 "
            f"v28_accepted={active} v28_committed={committed} "
            f"v28_commit_source={commit_source} v27_accepted=1 v27_ratio=0.8 "
            "v28_candidates=3 v28_evaluations=10 "
            f"v28_linf={linf} v28_l2={l2} v28_rho={ratio} "
            "v28_jacobi_weight=0.5 v28_step=0.3 "
            f"v28_radius_before=0.5 v28_radius_after={radius_after} "
            f"v26_l2_before=0.6 old_merit=0.9 new_merit=0.6 alpha=1 "
            f"selected={selected} qualified={qualified}\n"
            "CFDX nonlinear iteration 1 rejected; retry=1\n"
            "SST_JFNK_V22 iteration=10 v28_active=1 qualified=1\n")

def test_accepts_uncommitted_better_direction_without_claiming_convergence():
    result=first_attempt(trace())
    assert result["v28_accepted"] and not result["v28_committed"]
    assert not result["bfs_converged"] and not result["nasa_validated"]

def test_allows_real_committed_step_and_radius_update():
    result=first_attempt(trace(selected="sst_adaptive_trust_v28",committed="1",
                               radius_after="1.0",commit_source="v28"))
    assert result["v28_committed"] and result["radius_after"]==1.0

def test_v27_commit_can_expand_radius_when_v27_wins_final_comparison():
    result=first_attempt(trace(selected="sst_subspace_trust_v27",
                               committed="1",commit_source="v27",
                               radius_after="1.0"))
    assert result["v28_committed"]

def test_v27_can_commit_when_blend_has_no_admissible_candidate():
    r=first_attempt(trace(active="0",selected="sst_subspace_trust_v27",
                          committed="1",commit_source="v27",
                          radius_after="1.0"))
    assert r["v28_committed"] and not r["v28_accepted"]

def test_rejects_fake_v27_commit_source():
    with pytest.raises(ValueError):
        first_attempt(trace(selected="jfnk",committed="1",
                            commit_source="v27"))

def test_rejects_hypothetical_commit_and_unsupported_branch():
    with pytest.raises(ValueError):
        first_attempt(trace(committed="1"))
    with pytest.raises(ValueError):
        first_attempt(trace(selected="unknown"))
    with pytest.raises(ValueError):
        first_attempt(trace(active="0",selected="sst_adaptive_trust_v28"))

def test_rejects_l2_only_and_linf_only_candidates():
    with pytest.raises(ValueError):
        first_attempt(trace(l2="0.8"))
    with pytest.raises(ValueError):
        first_attempt(trace(linf="1.0"))
    with pytest.raises(ValueError):
        first_attempt(trace(ratio="nan"))

def test_rejects_uncommitted_radius_change_and_forged_qualification():
    with pytest.raises(ValueError):
        first_attempt(trace(radius_after="1.0"))
    with pytest.raises(ValueError):
        first_attempt(trace(qualified="1"))
    with pytest.raises(ValueError):
        first_attempt("SST_JFNK_V22 iteration=9 v28_active=1 qualified=0")

from compare_bfs_sst_v28 import trajectory

def test_v28_trajectory_refuses_false_outer_convergence_and_merit_growth():
    t=trace()
    with pytest.raises(ValueError):
        trajectory(t.replace("new_merit=0.6", "new_merit=1.5"),True)
    with pytest.raises(ValueError):
        trajectory(t.replace("v28_active=1", "v28_active=0"),True)
    with pytest.raises(ValueError):
        trajectory(t.split("CFDX nonlinear iteration 1 rejected",1)[0],True)

def test_v28_trajectory_marks_rejected_outer_state_even_when_trust_candidate_exists():
    result=trajectory(trace(),True)
    assert result["first_attempt_rejected"]
    assert not result["outer_accepted"] and not result["bfs_qualified"]
