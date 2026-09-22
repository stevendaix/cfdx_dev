import pytest

from cfdx.derived import derived_magnitude, derived_sum


def test_vector_magnitude() -> None:
    assert derived_magnitude({"x": (3, 0), "y": (4, 0)}) == (5.0, 0.0)


def test_derived_sum() -> None:
    assert derived_sum((1, 2), (3, 4)) == (4.0, 6.0)


def test_mismatched_fields_are_rejected() -> None:
    with pytest.raises(ValueError, match="lengths"):
        derived_sum((1,), (1, 2))
