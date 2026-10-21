"""End-to-end qualification of OpenFOAM ASCII result import through the CLI."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np


def _header(object_name: str, class_name: str) -> str:
    return (
        "FoamFile\n{\n"
        "    version 2.0;\n"
        "    format ascii;\n"
        f"    class {class_name};\n"
        f"    object {object_name};\n"
        "}\n"
    )


def _write_hex_case(case: Path) -> None:
    """Write a single valid 3-D hexahedral OpenFOAM cell with two patches."""
    poly = case / "constant" / "polyMesh"
    poly.mkdir(parents=True)
    points = [
        (0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
        (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
    ]
    faces = [
        (0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4),
        (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7),
    ]
    (poly / "points").write_text(
        _header("points", "vectorField")
        + f"{len(points)}\n(\n"
        + "".join(f"({' '.join(map(str, point))})\n" for point in points)
        + ")\n",
        encoding="utf-8",
    )
    (poly / "faces").write_text(
        _header("faces", "faceList")
        + f"{len(faces)}\n(\n"
        + "".join(f"{len(face)}({' '.join(map(str, face))})\n" for face in faces)
        + ")\n",
        encoding="utf-8",
    )
    (poly / "owner").write_text(
        _header("owner", "labelList") + "6\n(\n0\n0\n0\n0\n0\n0\n)\n",
        encoding="utf-8",
    )
    (poly / "neighbour").write_text(
        _header("neighbour", "labelList") + "0\n(\n)\n", encoding="utf-8"
    )
    (poly / "boundary").write_text(
        _header("boundary", "polyBoundaryMesh")
        + "2\n(\n"
        "    inlet\n    {\n        type patch;\n"
        "        nFaces 3;\n        startFace 0;\n    }\n"
        "    outlet\n    {\n        type wall;\n"
        "        nFaces 3;\n        startFace 3;\n    }\n"
        ")\n",
        encoding="utf-8",
    )


def _write_velocity(
    case: Path, time: str, values: list[tuple[float, float, float]]
) -> None:
    folder = case / time
    folder.mkdir()
    value_list = "\n".join(f"    ({x} {y} {z})" for x, y, z in values)
    (folder / "U").write_text(
        _header("U", "volVectorField")
        + "dimensions [0 1 -1 0 0 0 0];\n"
        + f"internalField nonuniform List<vector> {len(values)}\n(\n"
        + value_list
        + "\n);\n"
        + "boundaryField\n{\n"
        + "    inlet { type fixedValue; value uniform (0 0 0); }\n"
        + "    outlet { type zeroGradient; }\n"
        + "}\n",
        encoding="utf-8",
    )


def test_openfoam_results_cli_imports_latest_cell_velocity(tmp_path: Path) -> None:
    case = tmp_path / "foam_case"
    _write_hex_case(case)
    _write_velocity(case, "2", [(1.0, 2.0, 3.0)])
    expected = np.asarray([[4.25, -0.5, 8.0]], dtype=np.float64)
    _write_velocity(case, "10", [tuple(expected[0])])

    output = tmp_path / "converted"
    env = os.environ.copy()
    package_root = str(Path(__file__).resolve().parents[2] / "src" / "cfdx" / "python")
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (package_root, env.get("PYTHONPATH", "")) if part
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "cfdx.io.cli",
            "convert",
            str(case),
            "--results",
            "--output",
            str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
        env=env,
    )
    assert completed.returncode == 0, (
        f"OpenFOAM conversion CLI failed (exit {completed.returncode}).\n"
        f"stdout:\n{completed.stdout}\nstderr:\n{completed.stderr}"
    )
    assert "OpenFOAM time 10" in completed.stdout

    case_h5 = output / "case.cfdx.h5"
    assert case_h5.is_file()
    with h5py.File(case_h5, "r") as h5:
        assert int(h5.attrs["n_points"]) == 8
        assert int(h5.attrs["n_faces"]) == 6
        assert int(h5.attrs["n_cells"]) == 1
        np.testing.assert_array_equal(h5["owner"][:], np.zeros(6, dtype=np.uint64))
        assert h5["neighbour"].shape == (6,)
        assert set(h5["fields/vector"].keys()) == {"U"}
        np.testing.assert_allclose(h5["fields/vector/U"][:], expected, rtol=0, atol=0)
        assert "boundary_patches" in h5.attrs
        assert "inlet" in h5.attrs["boundary_patches"]
        assert "outlet" in h5.attrs["boundary_patches"]


def test_openfoam_results_cli_rejects_wrong_velocity_cell_count(tmp_path: Path) -> None:
    case = tmp_path / "foam_case"
    _write_hex_case(case)
    _write_velocity(case, "10", [(1.0, 2.0, 3.0), (4.0, 5.0, 6.0)])

    env = os.environ.copy()
    package_root = str(Path(__file__).resolve().parents[2] / "src" / "cfdx" / "python")
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (package_root, env.get("PYTHONPATH", "")) if part
    )
    completed = subprocess.run(
        [
            sys.executable,
            "-m",
            "cfdx.io.cli",
            "convert",
            str(case),
            "--results",
            "--output",
            str(tmp_path / "converted"),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
        env=env,
    )
    assert completed.returncode != 0
    assert "2 values" in completed.stdout or "2 values" in completed.stderr
