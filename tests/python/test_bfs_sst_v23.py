from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from diagnose_bfs_sst_v23 import compare, first_retry_samples  # noqa: E402


def fake_log(merit: float, *, block: bool, selected: str = "picard") -> str:
    return "\n".join([
        "SST_JFNK_V22 iteration=10 proposed=1 linear_true_relative=0.1 "
        f"block_probe={int(block)} block_proposed={int(block)} "
        f"block_linear_true_relative=0.22 coupled_blocks={712 if block else 0} "
        f"selected={selected}",
        "SST_TRANSPORT_DIAGNOSTIC iteration=10 "
        "k_returned_nonlinear_raw=0.002 omega_returned_nonlinear_raw=0.001 "
        f"k_local_relative={merit} omega_local_relative=0.2",
        "SST_ACCEPTANCE_GATES iteration=10 tolerance=1e-8 "
        "outer_update_ok=0 returned_equations_ok=0",
    ])


def test_compare_does_not_promote_unit_or_10_iteration_result() -> None:
    result = compare(fake_log(0.89, block=False), fake_log(0.89, block=True))
    assert result["v23_nonworse_merit_at_iteration_10"]
    assert not result["physical_convergence_certified"]
    assert not result["nasa_certified"]
    assert result["v23"]["coupled_blocks"] == 712


def test_compare_detects_block_regression() -> None:
    result = compare(fake_log(0.89, block=False), fake_log(0.96, block=True))
    assert not result["v23_nonworse_merit_at_iteration_10"]


def test_uncertified_block_direction_is_never_accepted() -> None:
    text = fake_log(0.2, block=True, selected="jfnk_block").replace(
        "block_linear_true_relative=0.22", "block_linear_true_relative=0.90"
    )
    with pytest.raises(ValueError):
        first_retry_samples(text)


def test_same_tolerance_required() -> None:
    with pytest.raises(ValueError):
        compare(fake_log(0.89, block=False),
                fake_log(0.89, block=True).replace("1e-8", "1e-7"))


def test_no_unexercised_block_promotion() -> None:
    with pytest.raises(ValueError):
        compare(fake_log(0.89, block=False), fake_log(0.89, block=False))
