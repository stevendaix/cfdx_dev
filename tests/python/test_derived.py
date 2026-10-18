import pytest
from cfdx.derived import derived_magnitude, derived_sum, derived_vector_magnitude


def test_vector_magnitude() -> None:
    assert derived_magnitude({"x": (3, 0), "y": (4, 0)}) == (5.0, 0.0)


def test_vector_field_magnitude_per_cell() -> None:
    """The production solver stores U as one cell field with three components.

    ``U -> [[u_x, u_y, u_z], ...]`` so the magnitude is computed per cell from
    the vector array, matching how the dataset actually stores the velocity.
    """
    assert derived_vector_magnitude([[3.0, 0.0, 0.0], [0.0, 4.0, 0.0]]) == (3.0, 4.0)


def test_vector_field_magnitude_rejects_empty() -> None:
    with pytest.raises(ValueError, match="at least one cell"):
        derived_vector_magnitude([])


def test_vector_field_magnitude_rejects_mismatched_components() -> None:
    with pytest.raises(ValueError, match="lengths"):
        derived_vector_magnitude([[3.0, 0.0], [4.0, 0.0, 0.0]])


def test_derived_sum() -> None:
    assert derived_sum((1, 2), (3, 4)) == (4.0, 6.0)


def test_mismatched_fields_are_rejected() -> None:
    with pytest.raises(ValueError, match="lengths"):
        derived_sum((1,), (1, 2))
