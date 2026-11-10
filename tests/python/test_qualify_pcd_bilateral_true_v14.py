from __future__ import annotations

import pytest

from scripts.qualify_pcd_bilateral_true_v14 import PcdProofError, certify


def action(
    *,
    projected=4e-9,
    aux=5e-9,
    floor=0.06,
    original=0.06,
    accepted=1,
    operator="deflated",
):
    return (
        "COUPLED_NESTED_PCD_FP bilateral_deflation=1 "
        f"auxiliary_operator={operator} preconditioner_anchored=1 "
        f"auxiliary_true_relative_rhs={aux} projected_original_true_relative={projected} "
        f"left_incompatible_floor={floor} original_true_relative_rhs={original} "
        f"accepted={accepted} status=0 iterations=27"
    )


def test_qualified_projected_not_physics():
    out = certify("\n".join([action()] * 12))
    assert out["actions"] == 12 and out["projected_preconditioner_qualified"]
    assert not out["outer_coupled_converged"] and not out["driver_seegmiller_qualified"]


@pytest.mark.parametrize(
    "line",
    [
        action(operator="anchored"),
        action(projected=1e-2),
        action(aux=1e-1),
        action(original=0),
        action(accepted=0),
        action(floor=float("nan")),
        "COUPLED_NESTED_PCD_FP status=0 accepted=1",
    ],
)
def test_reject_invalid_action(line):
    with pytest.raises(PcdProofError):
        certify("\n".join([line] * 12))


def test_reject_insufficient_evidence():
    with pytest.raises(PcdProofError):
        certify("\n".join([action()] * 2))
