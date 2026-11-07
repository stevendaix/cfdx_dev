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
