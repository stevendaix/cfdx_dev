from __future__ import annotations
import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2] / "scripts"))
from diagnose_bfs_sst_v24 import first_attempt,compare  # noqa: E402


def fake(spatial:bool, selected:str="jfnk", tolerance:str="1e-8", linear:str="0.16") -> str:
    return "\n".join([
        "SST_JFNK_V22 iteration=10 proposed=1 linear_true_relative=0.24 "
        f"block_probe=1 block_proposed=0 block_linear_true_relative=0.74 "
        f"spatial_probe={int(spatial)} spatial_proposed=0 "
        f"spatial_linear_true_relative={linear} spatial_offcell={5256 if spatial else 0} "
        f"selected={selected}",
        "SST_TRANSPORT_DIAGNOSTIC iteration=10 k_returned_nonlinear_raw=1e-6 "
        "omega_returned_nonlinear_raw=1.7 k_local_relative=0.8899 omega_local_relative=0.2",
        f"SST_ACCEPTANCE_GATES iteration=10 tolerance={tolerance} "
        "outer_update_ok=0 returned_equations_ok=0",
        "SST_JACOBIAN_AUDIT_V24 iteration=10 cell=357 max_relative_difference=0.1154",
        "CFDX nonlinear iteration 1 rejected",
    ])


def test_v24_report_remains_nonqualified():
    r=compare(fake(False),fake(True))
    assert r["v24_spatial_linear_better_than_v23_block"]
    assert r["v24_merit_not_worse_at_iteration_10"]
    assert r["v24"]["jacobian_errors_at_iteration_10"]["357"] == pytest.approx(0.1154)
    assert not r["sst_full_convergence_certified"]
    assert not r["bfs_global_convergence_certified"]
    assert not r["nasa_physical_qualification_certified"]


def test_v24_never_promotes_unconverged_spatial_direction():
    with pytest.raises(ValueError):
        first_attempt(fake(True,selected="jfnk_spatial"))
    with pytest.raises(ValueError):
        compare(fake(False),fake(False))


def test_v24_fails_closed_when_metrics_or_tolerance_differ():
    with pytest.raises(ValueError):
        compare(fake(False),fake(True,tolerance="1e-7"))
    with pytest.raises(ValueError):
        first_attempt(fake(True).replace("omega_local_relative=0.2","omega_local_relative=nan"))
    with pytest.raises(ValueError):
        first_attempt("SST_JFNK_V22 iteration=1 selected=picard")


def test_v24_ignores_later_retries():
    text=fake(True).replace("CFDX nonlinear iteration 1 rejected","CFDX nonlinear iteration 1 rejected\nSST_JFNK_V22 iteration=10 selected=jfnk_spatial")
    assert first_attempt(text)["selected"]=="jfnk"
