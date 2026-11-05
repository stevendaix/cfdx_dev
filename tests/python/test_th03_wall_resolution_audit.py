"""TH-03 mesh qualification has an independent wall-resolution gate."""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from scripts.audit_th03_wall_resolution import WallResolutionError, audit


def fixture(
    tmp_path: Path, *, speed: float = 1.0, declared: bool = True
) -> tuple[Path, Path]:
    cells = tmp_path / "cells.csv"
    faces = tmp_path / "faces.csv"
    with cells.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["cell_id", "nu_t_m2_s"])
        writer.writeheader()
        writer.writerow({"cell_id": 0, "nu_t_m2_s": 9.0})
    columns = [
        "owner",
        "patch",
        "wall_normal_distance_m",
        "wall_shear_x_Pa",
        "wall_shear_y_Pa",
        "wall_shear_z_Pa",
        "wall_shear_method",
        "wall_y_plus_estimate",
        "owner_nu_t_over_nu",
    ]
    with faces.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerow(
            {
                "owner": 0,
                "patch": "downstream_wall",
                "wall_normal_distance_m": 0.1,
                "wall_shear_x_Pa": speed**2,
                "wall_shear_y_Pa": 0,
                "wall_shear_z_Pa": 0,
                "wall_shear_method": "first_order_one_sided_diagnostic",
                "wall_y_plus_estimate": (0.1 * speed if declared else 999),
                "owner_nu_t_over_nu": 9,
            }
        )
    return faces, cells


def test_fine_wall_preflight_is_not_physical_qualification(tmp_path: Path) -> None:
    f, c = fixture(tmp_path)
    r = audit(f, c, wall_patches=("downstream_wall",), density=1, viscosity=1)
    assert r["status"] == "PREFLIGHT_PASS"
    assert r["physical_qualification"] == "NOT_QUALIFIED"
    assert r["max_y_plus_estimate"] == pytest.approx(0.1)


def test_underresolved_wall_is_blocked(tmp_path: Path) -> None:
    f, c = fixture(tmp_path, speed=50)
    r = audit(f, c, wall_patches=("downstream_wall",), density=1, viscosity=1)
    assert r["status"] == "BLOCKED" and r["max_y_plus_estimate"] == pytest.approx(5)


def test_corrupted_face_yplus_is_rejected(tmp_path: Path) -> None:
    f, c = fixture(tmp_path, declared=False)
    with pytest.raises(WallResolutionError, match="differs"):
        audit(f, c, wall_patches=("downstream_wall",), density=1, viscosity=1)


def test_missing_patch_fails_closed(tmp_path: Path) -> None:
    f, c = fixture(tmp_path)
    with pytest.raises(WallResolutionError, match="no matching wall"):
        audit(f, c, wall_patches=("step_wall",), density=1, viscosity=1)


def test_reject_negative_or_nan_turbulent_viscosity(tmp_path: Path) -> None:
    f, c = fixture(tmp_path)
    for invalid in ("-1", "nan", "inf"):
        c.write_text("cell_id,nu_t_m2_s\n0," + invalid + "\n", encoding="utf-8")
        with pytest.raises(WallResolutionError, match="viscosity"):
            audit(f, c, wall_patches=("downstream_wall",), density=1, viscosity=1)


def test_every_requested_wall_patch_must_be_present(tmp_path: Path) -> None:
    f, c = fixture(tmp_path)
    with pytest.raises(WallResolutionError, match="missing requested wall patches"):
        audit(f, c, wall_patches=("downstream_wall", "top"), density=1, viscosity=1)
