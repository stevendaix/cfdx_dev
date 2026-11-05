"""Manufactured-solution and fail-closed tests for wall-gradient diagnostics."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts import compare_th03_wall_gradients as diag


def make_fixtures(tmp_path: Path, *, skew: bool = False) -> tuple[Path, Path]:
    faces = tmp_path / "faces.csv"
    cells = tmp_path / "cells.csv"
    with faces.open("w", newline="") as f:
        w = csv.DictWriter(
            f,
            fieldnames=[
                "patch",
                "x_m",
                "y_m",
                "z_m",
                "area_x_m2",
                "area_y_m2",
                "area_z_m2",
                "area_m2",
                "wall_normal_distance_m",
                "wall_shear_x_Pa",
                "wall_shear_method",
            ],
        )
        w.writeheader()
        for x, a in ((1, -4), (2, -2), (3, 2), (4, 4)):
            u1 = a * 0.1 + 5 * 0.1**2
            w.writerow(
                {
                    "patch": "downstream_wall",
                    "x_m": x,
                    "y_m": 0,
                    "z_m": 0.5,
                    "area_x_m2": 1 if skew else 0,
                    "area_y_m2": -1,
                    "area_z_m2": 0,
                    "area_m2": 1,
                    "wall_normal_distance_m": 0.1,
                    "wall_shear_x_Pa": 0.02 * u1 / 0.1,
                    "wall_shear_method": "first_order_one_sided_diagnostic",
                }
            )
    with cells.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["x_m", "y_m", "z_m", "u_x_m_s"])
        w.writeheader()
        for x, a in ((1, -4), (2, -2), (3, 2), (4, 4)):
            for d in (0.1, 0.2):
                w.writerow({"x_m": x, "y_m": d, "z_m": 0.5, "u_x_m_s": a * d + 5 * d * d})
    return faces, cells


def test_manufactured_quadratic_wall_derivative_and_root(tmp_path: Path) -> None:
    faces, cells = make_fixtures(tmp_path)
    result = diag.compare(
        faces,
        cells,
        patch="downstream_wall",
        density=2,
        viscosity=0.01,
        bulk_speed=3,
        step_height=1,
    )
    assert result["quadratic_reattachment_x_over_h"] == pytest.approx(2.5, abs=1e-12)
    assert result["first_order_reattachment_x_over_h"] < 2.5
    assert result["reconstruction_spread_x_over_h"] > 0
    assert result["physical_qualification"] == "NOT_QUALIFIED"


def test_quadratic_wall_gradient_rejects_bad_geometry_and_stencil(
    tmp_path: Path,
) -> None:
    faces, cells = make_fixtures(tmp_path, skew=True)
    with pytest.raises(diag.WallFitError, match="planar"):
        diag.compare(
            faces,
            cells,
            patch="downstream_wall",
            density=2,
            viscosity=0.01,
            bulk_speed=3,
            step_height=1,
        )
    with pytest.raises(diag.WallFitError, match="stencil"):
        diag.quadratic_dudn(0.1, 1, 0.1, 1)
