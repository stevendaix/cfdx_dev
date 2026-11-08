"""Fail-closed NACA farfield face-flux diagnostics, including signed cancellation."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import h5py
import numpy as np
import pytest

MODULE = Path(__file__).resolve().parents[2] / "scripts" / "audit_naca_farfield_flux.py"
SPEC = importlib.util.spec_from_file_location("audit_naca_farfield_flux", MODULE)
assert SPEC and SPEC.loader
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)

FIELDS = AUDIT.REQUIRED


def fixture_case(tmp_path: Path, *, matched: bool = False) -> tuple[Path, Path]:
    mesh = tmp_path / "mesh.h5"
    with h5py.File(mesh, "w") as h:
        h.attrs["boundary_patches"] = "airfoil:0:1:0;farfield:1:2:7"
        h.attrs["n_cells"] = "1"
        h.attrs["n_faces"] = "3"
        h.attrs["mesh_hash"] = "fixture-geometry"
        h.attrs["topology_hash"] = "fixture-topology"
        h.create_dataset("owner", data=np.array([0, 0, 0], dtype=np.int64))
        h.create_dataset("neighbour", data=np.array([-1, -1, -1], dtype=np.int64))
        h.create_dataset(
            "patch_face_offsets", data=np.array([0, 1, 3], dtype=np.uint64)
        )
        h.create_dataset("patch_face_ids", data=np.array([0, 1, 2], dtype=np.uint64))
    filename = tmp_path / "flux.csv"
    data = [
        ["0", "0", "-1", "0", "0", "0", "-1", "0", "0", "0"],
        ["1", "0", "-1", "20", "1", "0", "1", "0", "0", "1" if matched else "0.9"],
        [
            "2",
            "0",
            "-1",
            "20",
            "-1",
            "0",
            "-0.5",
            "0",
            "0",
            "-0.5" if matched else "-0.4",
        ],
    ]
    with filename.open("w", newline="") as out:
        w = csv.writer(out)
        w.writerow(FIELDS)
        w.writerows(data)
    return mesh, filename


def test_boundary_face_mismatch_does_not_cancel(tmp_path: Path) -> None:
    mesh, flux = fixture_case(tmp_path)
    result, rows = AUDIT.analyze(mesh, flux)
    assert result["selected_face_count"] == 2
    assert abs(result["mismatch_sum_signed"]) < 1e-14
    assert result["mismatch_l1"] == pytest.approx(0.2)
    assert result["mismatch_linf"] == pytest.approx(0.1)
    assert result["compatible_to_absolute_1e-10"] is False
    assert {row["face"] for row in rows} == {1, 2}


def test_neumann_like_prescribed_velocity_boundary_is_exact(tmp_path: Path) -> None:
    mesh, flux = fixture_case(tmp_path, matched=True)
    result, _ = AUDIT.analyze(mesh, flux)
    assert result["compatible_to_absolute_1e-10"] is True
    assert result["mismatch_linf"] == 0


@pytest.mark.parametrize(
    "altered",
    ["duplicate", "missing", "extra", "owner", "neighbour", "nonfinite", "columns"],
)
def test_invalid_face_flux_fails_closed(tmp_path: Path, altered: str) -> None:
    mesh, flux = fixture_case(tmp_path)
    with flux.open(newline="") as stream:
        rows = list(csv.reader(stream))
    if altered == "duplicate":
        rows[2][0] = "0"
    elif altered == "missing":
        rows.pop()
    elif altered == "extra":
        rows.append(rows[-1])
    elif altered == "owner":
        rows[1][1] = "2"
    elif altered == "neighbour":
        rows[1][2] = "0"
    elif altered == "nonfinite":
        rows[2][9] = "nan"
    elif altered == "columns":
        rows[0][5] = "absent"
    with flux.open("w", newline="") as stream:
        csv.writer(stream).writerows(rows)
    with pytest.raises(ValueError):
        AUDIT.analyze(mesh, flux)


def test_missing_patch_or_corrupt_patch_indices_rejected(tmp_path: Path) -> None:
    mesh, flux = fixture_case(tmp_path)
    with pytest.raises(ValueError):
        AUDIT.analyze(mesh, flux, patch="unknown")
    with h5py.File(mesh, "r+") as h:
        h["patch_face_ids"][2] = 1
    with pytest.raises(ValueError):
        AUDIT.analyze(mesh, flux)


def test_invalid_physical_inputs_rejected(tmp_path: Path) -> None:
    mesh, flux = fixture_case(tmp_path)
    with pytest.raises(ValueError):
        AUDIT.analyze(mesh, flux, density=-1)
    with pytest.raises(ValueError):
        AUDIT.analyze(mesh, flux, velocity=(float("inf"), 0.0, 0.0))
