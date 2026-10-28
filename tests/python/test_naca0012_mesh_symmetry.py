"""Regression tests for the generated NACA0012 O-grid geometry and diagnostics."""

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "scripts" / "generate_2d_validation_meshes.py"
SPEC = importlib.util.spec_from_file_location(
    "cfdx_generate_2d_validation_meshes", SCRIPT
)
assert SPEC is not None and SPEC.loader is not None
MESHGEN = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MESHGEN)


@pytest.mark.parametrize(
    ("n_surface", "n_radial"),
    [(64, 24), (128, 48), (256, 96)],
)
def test_naca0012_outer_boundary_is_reflection_symmetric(n_surface, n_radial):
    points, _, _ = MESHGEN.naca_o_grid(n_surface, n_radial, radius=20.0)
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


@pytest.mark.parametrize(
    ("n_surface", "n_radial"),
    [(64, 24), (128, 48), (256, 96)],
)
def test_naca0012_o_grid_quality_metrics_are_finite(n_surface, n_radial):
    quality = MESHGEN.naca_o_grid_quality(n_surface, n_radial, radius=20.0)

    assert quality["cell_count"] == (2 * n_surface + 1) * n_radial
    assert quality["inverted_or_degenerate_cells"] == 0
    assert quality["minimum_absolute_cell_area"] > 0.0
    assert quality["minimum_first_layer_spacing"] > 0.0
    assert (
        quality["maximum_first_layer_spacing"] >= quality["minimum_first_layer_spacing"]
    )
    assert quality["minimum_radial_growth_ratio"] > 0.0
    assert (
        quality["maximum_radial_growth_ratio"] >= quality["minimum_radial_growth_ratio"]
    )
    assert np.isfinite(quality["maximum_cell_edge_aspect_ratio"])
    assert quality["maximum_cell_edge_aspect_ratio"] >= 1.0
    assert 0.0 <= quality["p95_nonorthogonality_degrees"]
    assert (
        quality["p95_nonorthogonality_degrees"]
        <= quality["maximum_nonorthogonality_degrees"]
        <= 90.0
    )


@pytest.mark.parametrize(
    ("n_surface", "n_radial", "radius"),
    [(1, 24, 20.0), (64, 0, 20.0), (64, 24, 0.0), (64, 24, float("inf"))],
)
def test_naca0012_o_grid_rejects_invalid_parameters(n_surface, n_radial, radius):
    with pytest.raises(ValueError):
        MESHGEN.naca_o_grid(n_surface, n_radial, radius)


def test_naca0012_first_layer_spacing_refines_monotonically():
    levels = [(64, 24), (128, 48), (256, 96)]
    reports = [
        MESHGEN.naca_o_grid_quality(n_surface, n_radial, radius=20.0)
        for n_surface, n_radial in levels
    ]

    spacings = [report["mean_first_layer_spacing"] for report in reports]
    assert spacings[0] > spacings[1] > spacings[2]
    assert all(report["inverted_or_degenerate_cells"] == 0 for report in reports)


def test_naca0012_quality_report_is_valid_json(tmp_path):
    subprocess.run(
        [
            sys.executable,
            str(SCRIPT),
            "naca0012",
            "--quick",
            "--output-dir",
            str(tmp_path),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(
        (tmp_path / "naca0012_mesh_quality.json").read_text(encoding="utf-8")
    )
    assert len(report) == 1
    assert report[0]["inverted_or_degenerate_cells"] == 0
    assert "maximum_nonorthogonality_degrees" in report[0]
