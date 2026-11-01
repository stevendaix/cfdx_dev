"""Regression tests for the wall-normal NACA0012 O-grid construction.

These tests qualify geometric invariants only; they do not establish CFD
convergence, force accuracy, or spatial discretization uncertainty.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "generate_2d_validation_meshes.py"
SPEC = importlib.util.spec_from_file_location(
    "cfdx_generate_2d_validation_meshes_wallnormal", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MESHGEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MESHGEN)


@pytest.mark.parametrize(("n_surface", "n_radial"), [(64, 24), (128, 48), (256, 96)])
def test_naca0012_wall_normal_grid_preserves_symmetry_and_valid_cells(
    n_surface, n_radial
):
    points, cells, _ = MESHGEN.naca_o_grid(n_surface, n_radial, radius=20.0)
    points = np.asarray(points, dtype=np.float64)
    m = 2 * n_surface + 1
    planar = points[: (n_radial + 1) * m, :2].reshape(n_radial + 1, m, 2)

    np.testing.assert_allclose(planar[:, :, 0], planar[:, ::-1, 0], atol=5e-12)
    np.testing.assert_allclose(planar[:, :, 1], -planar[:, ::-1, 1], atol=5e-12)
    np.testing.assert_allclose(np.linalg.norm(planar[-1], axis=1), 20.0, atol=5e-12)

    report = MESHGEN.naca_o_grid_quality(n_surface, n_radial, radius=20.0)
    assert len(cells) == m * n_radial
    assert report["inverted_or_degenerate_cells"] == 0
    assert report["minimum_absolute_cell_area"] > 0.0
    assert np.isfinite(report["maximum_nonorthogonality_degrees"])


@pytest.mark.parametrize(("n_surface", "n_radial"), [(64, 24), (128, 48), (256, 96)])
def test_naca0012_first_layer_follows_outward_surface_normal(n_surface, n_radial):
    points, _, _ = MESHGEN.naca_o_grid(n_surface, n_radial, radius=20.0)
    points = np.asarray(points, dtype=np.float64)
    surface = np.asarray(MESHGEN.naca0012_surface(n_surface), dtype=np.float64)
    m = len(surface)
    planar = points[: (n_radial + 1) * m, :2].reshape(n_radial + 1, m, 2)

    alignment = []
    for i in range(m):
        tangent = surface[(i + 1) % m] - surface[(i - 1) % m]
        tangent /= np.linalg.norm(tangent)
        outward_normal = np.array((tangent[1], -tangent[0]))
        first_layer = planar[1, i] - planar[0, i]
        first_layer /= np.linalg.norm(first_layer)
        alignment.append(float(np.dot(first_layer, outward_normal)))

    assert min(alignment) > 0.90
    assert float(np.mean(alignment)) > 0.97


@pytest.mark.parametrize(("n_surface", "n_radial"), [(64, 24), (128, 48), (256, 96)])
def test_naca0012_wall_normal_grid_controls_peak_nonorthogonality(
    n_surface, n_radial
):
    report = MESHGEN.naca_o_grid_quality(n_surface, n_radial, radius=20.0)

    # The former radial-blend family reached 73.24 degrees at N256.
    # Keep a broad bound to detect a regression without overfitting the mesh.
    assert report["maximum_nonorthogonality_degrees"] < 45.0
    assert report["p95_nonorthogonality_degrees"] < 20.0
    assert report["maximum_radial_growth_ratio"] < 2.6
    assert report["inverted_or_degenerate_cells"] == 0
