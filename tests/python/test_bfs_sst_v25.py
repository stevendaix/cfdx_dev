import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
from diagnose_bfs_sst_v25 import first_attempt


def fake(*, physical_step=0.9259, selected="jfnk", certified="0", v25="1"):
    return "\n".join([
        "SST_JFNK_V22 iteration=10 spatial_probe=1 spatial_linear_true_relative=0.166 "
        f"v25_active={v25} v25_direction=1 v25_accepted=1 v25_alpha=0.001 "
        f"v25_best_trial_merit={physical_step} v25_wall_corrected=1 "
        f"old_merit=0.926 new_merit=0.889 selected={selected}",
        "SST_TRANSPORT_DIAGNOSTIC iteration=10 k_local_relative=0.889 omega_local_relative=0.2",
        f"SST_ACCEPTANCE_GATES iteration=10 tolerance=1e-8 qualified_bfs={certified}",
        "CFDX nonlinear iteration 1 rejected; retry=1",
        "SST_JFNK_V22 iteration=10 v25_accepted=0",  # ignored retries
    ])


def test_v25_certificates_remain_false_and_initial_attempt_only():
    r=first_attempt(fake())
    assert r["v25_physical_log_step_found"]
    assert r["spatial_krylov_relative"]==pytest.approx(.166)
    assert r["selected_direction"]=="jfnk"
    assert r["wall_blocks_corrected"]==1
    assert r["sst_full_convergence_certified"] is False
    assert r["bfs_global_convergence_certified"] is False
    assert r["nasa_physical_validation_certified"] is False


def test_v25_refuses_physical_false_promotion():
    with pytest.raises(ValueError): first_attempt(fake(certified="1"))
    with pytest.raises(ValueError): first_attempt(fake(v25="0"))
    with pytest.raises(ValueError): first_attempt(fake(physical_step=1.1))


def test_v25_refuses_missing_or_corrupted_evidence():
    with pytest.raises(ValueError): first_attempt("SST_JFNK_V22 iteration=10")
    with pytest.raises(ValueError): first_attempt(fake().replace("0.166","nan"))
    with pytest.raises(ValueError): first_attempt(fake().replace("0.889","inf"))


def test_v25_does_not_select_unavailable_candidate():
    with pytest.raises(ValueError):
        first_attempt(fake(selected="jfnk_spatial_log_v25").replace(
            "v25_accepted=1", "v25_accepted=0"))
