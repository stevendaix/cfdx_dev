import pytest

from cfdx.surface import surface_integral, vector_surface_integral


def test_surface_integral() -> None:
    assert surface_integral((2, 3), (0.5, 2.0)) == 7.0


def test_vector_surface_integral() -> None:
    assert vector_surface_integral(((1, 2), (3, 4)), (2, 1)) == (5.0, 8.0)


def test_surface_shape_is_validated() -> None:
    with pytest.raises(ValueError):
        surface_integral((1,), (1, 2))
