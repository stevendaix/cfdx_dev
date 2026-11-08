"""Radial mesh nesting must be geometric, not merely a JSON claim."""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts" / "generate_2d_validation_meshes.py"
)
SPEC = importlib.util.spec_from_file_location("cfdx_nested_naca_mesh", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
GEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(GEN)


def rings(n_surface: int, n_radial: int) -> np.ndarray:
    points, _, _ = GEN.naca_o_grid(
        n_surface, n_radial, radius=20.0, base_radial_layers=24
    )
    m = 2 * n_surface + 1
    return np.asarray(points)[: (n_radial + 1) * m, :2].reshape(n_radial + 1, m, 2)


def test_common_radial_nodes_and_first_layer_spacing_are_exactly_nested():
    # Keep fixed tangential positions: the full N64/N128/N256 family also
    # changes tangential geometry, so only the radial subfamily is nested.
    coarse, medium, fine = rings(32, 24), rings(32, 48), rings(32, 96)
    np.testing.assert_allclose(medium[::2], coarse, atol=2e-14, rtol=0.0)
    np.testing.assert_allclose(fine[::2], medium, atol=2e-14, rtol=0.0)
    coarse_spacing = np.linalg.norm(coarse[1] - coarse[0], axis=1)
    medium_spacing = np.linalg.norm(medium[1] - medium[0], axis=1)
    fine_spacing = np.linalg.norm(fine[1] - fine[0], axis=1)
    np.testing.assert_allclose(coarse_spacing / medium_spacing, 2.0, atol=1e-11)
    np.testing.assert_allclose(medium_spacing / fine_spacing, 2.0, atol=1e-11)


def test_noninteger_base_refinement_fails_before_mesh_creation():
    with pytest.raises(ValueError, match="integer subdivision"):
        GEN.naca_o_grid(32, 48, base_radial_layers=25)


def test_trailing_edge_sampling_keeps_exact_profile_and_nested_x_nodes():
    coarse = np.asarray(GEN.naca0012_surface(64, trailing_edge_regularization=0.04))
    fine = np.asarray(GEN.naca0012_surface(128, trailing_edge_regularization=0.04))
    assert coarse.shape == (129, 2) and fine.shape == (257, 2)
    np.testing.assert_allclose(coarse[:65], fine[:129:2], atol=3e-15, rtol=0)
    # Actual NACA thickness polynomial at the redistributed x coordinates.
    x = fine[:, 0]
    t = (
        5
        * 0.12
        * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x * x
            + 0.2843 * x**3
            - 0.1015 * x**4
        )
    )
    np.testing.assert_allclose(np.abs(fine[:, 1]), t, atol=5e-15, rtol=0)
    np.testing.assert_allclose(fine[0], [1.0, 0.00126], atol=2e-15, rtol=0)
    np.testing.assert_allclose(fine[-1], [1.0, -0.00126], atol=2e-15, rtol=0)


def test_te_regularization_reduces_extreme_aspect_without_inverting_cells():
    baseline = GEN.naca_o_grid_quality(128, 48, base_radial_layers=24)
    modified = GEN.naca_o_grid_quality(
        128, 48, base_radial_layers=24, trailing_edge_regularization=0.04
    )
    assert baseline["maximum_cell_edge_aspect_ratio"] > 200
    assert modified["maximum_cell_edge_aspect_ratio"] < 80
    assert modified["inverted_or_degenerate_cells"] == 0
    assert modified["trailing_edge_regularization"] == 0.04
    assert modified["maximum_nonorthogonality_degrees"] < 70
    assert baseline["worst_aspect_radial_layer"] == 0
    assert baseline["worst_aspect_cell_centre_xy"][0] > 1.0
    assert baseline["cells_with_edge_aspect_above_100"] == 2
    assert modified["cells_with_edge_aspect_above_100"] == 0
    assert modified["p99_cell_edge_aspect_ratio"] < 20


def test_te_regularization_invalid_parameter_fails():
    for value in (-0.1, 0.09, float("nan"), float("inf")):
        with pytest.raises(ValueError, match="trailing-edge regularization"):
            GEN.naca0012_surface(64, trailing_edge_regularization=value)
