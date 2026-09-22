import pytest

from cfdx.probe import nearest_probe


def test_nearest_probe() -> None:
    assert nearest_probe(((0, 0), (2, 0), (0, 2)), (1, 4, 8), (1.8, 0.1)) == 4.0


def test_probe_validates_alignment() -> None:
    with pytest.raises(ValueError):
        nearest_probe(((0,),), (), (0,))
