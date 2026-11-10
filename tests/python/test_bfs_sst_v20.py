from pathlib import Path
import sys
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from diagnose_bfs_sst_v20 import analyze


def test_budget_and_committed_parsed_without_physical_qualification():
    line = ("SST_ASSEMBLY_V20 iteration=100 field=omega cell=77 convection=1 "
            "diffusion=-2 production=3 destruction=-4 cross=5 balance=3 "
            "closure_error=0 source_identity=0 constrained=0 qualified=0")
    state = ("GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=0.5 p_continuity=0.01 "
             "k=1e-12 omega=1e-12 merit=0.5 converged=0 qualified=0")
    out = analyze(line + "\n" + state)
    assert out["row_budgets"][0]["cell"] == 77
    assert len(out["committed_states"]) == 1
    assert not out["driver_seegmiller_qualified"]
    assert not out["three_state_contraction_proved"]


def test_incomplete_or_nonfinite_log_is_rejected():
    with pytest.raises(ValueError):
        analyze("SST_ASSEMBLY_V20 iteration=1 field=k")
    with pytest.raises(ValueError):
        analyze("SST_ASSEMBLY_V20 iteration=1 field=k cell=77 convection=nan "
                "diffusion=1 production=1 destruction=1 cross=1 balance=1 "
                "closure_error=0 source_identity=0")
