"""Negative and positive tests of state counting; no physics promotion."""
import pytest
from scripts.audit_bfs_commits_v17 import audit


def _v13(i: int, merit: float) -> str:
    return (f"GLOBAL_NONLINEAR_COMMITTED_V13 iteration={i} U={merit} "
            f"p_continuity=0.001 k=0.00001 omega=0.00001 nut_closure=0 "
            f"merit={merit} converged=0 qualified=0")


def _v14(i: int, merit: float) -> str:
    return (f"GLOBAL_CONTRACTION_COMMITTED_V14 iteration={i} committed={i} "
            f"contractions={i - 1} merit={merit} converged=0 qualified=0")


def test_one_committed_is_not_three():
    record = audit(_v13(1, 1.783))
    assert record["committed_states"] == 1
    assert not record["three_states_observed"]
    assert not record["strict_same_operator_contraction"]
    assert not record["reference_comparison_permitted"]


def test_three_v13_merits_do_not_prove_same_operator():
    data = "\n".join([_v13(i, v) for i, v in enumerate((1., .8, .6), 1)])
    record = audit(data)
    assert record["committed_states"] == 3
    assert not record["strict_same_operator_contraction"]


def test_v14_ledger_can_prove_contraction_but_not_physics():
    lines = [item for i, value in enumerate((1., .8, .6), 1)
             for item in (_v13(i, value), _v14(i, value))]
    record = audit("\n".join(lines))
    assert record["strict_same_operator_contraction"]
    assert not record["fully_converged"]
    assert not record["physical_qualified"]


def test_missing_or_invalid_evidence_refused():
    assert audit("BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_BEGIN")["committed_states"] == 0
    with pytest.raises(ValueError, match="duplicate"):
        audit("\n".join((_v13(1, 1), _v13(1, .5))))
    with pytest.raises(ValueError, match="nonfinite"):
        audit(_v13(1, float("nan")))
    with pytest.raises(ValueError, match="unsupported qualified"):
        audit(_v13(1, 1).replace("qualified=0", "qualified=1"))
