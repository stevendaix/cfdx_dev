from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from diagnose_bfs_sst_v21 import analyze


def test_projected_wall_accepted_but_not_physically_certified():
    text = ("SST_OWNER_SYNC_V21 iteration=100 owners=69 initial_gap=396.399 "
            "trial_gap=0 final_gap=0 final_relative=0 worst_cell=357 "
            "line_accepted=1 qualified=0\n"
            "SST_WALL_OWNER_CONSTRAINT_V18 iteration=100 owners=69 "
            "max_absolute_mismatch=0 worst_cell=0 qualified=0\n"
            "SST_TRANSPORT_DIAGNOSTIC iteration=100 "
            "k_returned_nonlinear_raw=0.0087 omega_returned_nonlinear_raw=24.4 "
            "k_returned_row_scaled=0.55 omega_returned_row_scaled=1451 "
            "k_local_relative=0.05 omega_local_relative=0.0097 "
            "k_local_worst_cell=357 omega_local_worst_cell=363\n")
    result = analyze(text)
    assert result["owner_projection_stable"]
    assert result["owner_sync"][0]["worst_cell"] == 357
    assert result["transport"][0]["omega_local_worst_cell"] == 363
    assert result["transport"][0]["omega_returned_row_scaled"] == 1451
    assert not result["three_state_contraction_proved"]
    assert not result["driver_seegmiller_qualified"]


def test_nonfinite_or_missing_v21_trace_fails_closed():
    with pytest.raises(ValueError):
        analyze("SST_OWNER_SYNC_V21 iteration=3 owners=69")
    with pytest.raises(ValueError):
        analyze("SST_OWNER_SYNC_V21 iteration=3 owners=69 initial_gap=1 "
                "trial_gap=1 final_gap=nan final_relative=0 worst_cell=357 "
                "line_accepted=1")


def test_unstable_or_rejected_or_absent_wall_projections_never_certify():
    assert not analyze("")["owner_projection_stable"]
    for accept, relative in ((0, 0), (1, 1e-8)):
        sample = ("SST_OWNER_SYNC_V21 iteration=3 owners=69 initial_gap=1 "
                  "trial_gap=1 final_gap=1e-8 final_relative="
                  f"{relative} worst_cell=357 line_accepted={accept}")
        assert not analyze(sample)["owner_projection_stable"]
