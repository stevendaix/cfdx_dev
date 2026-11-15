import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from qualify_bfs_global_contraction_v14 import QualificationError, certify


def record(i, merit, contractions=None):
    c = i - 1 if contractions is None else contractions
    return (
        f"GLOBAL_CONTRACTION_COMMITTED_V14 iteration={i} committed={i} "
        f"contractions={c} merit={merit} converged=0 qualified=0"
    )


def test_three_states_contract_without_claiming_physics():
    result = certify("\n".join([record(1, 0.8), record(2, 0.6), record(3, 0.4)]))
    assert result["contractions"] == 2
    assert result["bfs_driver_seegmiller_qualified"] is False


@pytest.mark.parametrize(
    "log",
    [
        record(1, 0.8),
        "\n".join([record(1, 0.8), record(2, 0.8), record(3, 0.4)]),
        "\n".join([record(1, 0.8), record(2, 0.6), record(3, 0.4, 0)]),
        "\n".join([record(1, 0.8), record(2, float("nan")), record(3, 0.4)]),
        "\n".join([record(1, 0.8), record(2, 0.6), record(2, 0.4)]),
    ],
)
def test_fail_closed(log):
    with pytest.raises(QualificationError):
        certify(log)
