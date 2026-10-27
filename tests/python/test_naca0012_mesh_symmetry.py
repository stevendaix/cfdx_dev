"""Regression tests for the generated NACA0012 O-grid outer boundary."""

from pathlib import Path
import sys

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from generate_2d_validation_meshes import naca_o_grid


@pytest.mark.parametrize(
    ("n_surface", "n_radial"),
    [(64, 24), (128, 48), (256, 96)],
)
def test_naca0012_outer_boundary_is_reflection_symmetric(n_surface, n_radial):
    points, _, _ = naca_o_grid(n_surface, n_radial, radius=20.0)
    points = np.asarray(points, dtype=np.float64)

    # extrude_2d stores all z=0 points first. The last radial ring has m
    # points, paired by i <-> m - 1 - i because the airfoil polygon has an
    # odd number of vertices (2*n_surface + 1).
    m = 2 * n_surface + 1
    outer = points[n_radial * m : (n_radial + 1) * m, :2]
    mirrored = outer[::-1]

    np.testing.assert_allclose(outer[:, 0], mirrored[:, 0], rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(outer[:, 1], -mirrored[:, 1], rtol=0.0, atol=1.0e-12)
    np.testing.assert_allclose(
        np.linalg.norm(outer, axis=1), 20.0, rtol=0.0, atol=1.0e-12
    )
