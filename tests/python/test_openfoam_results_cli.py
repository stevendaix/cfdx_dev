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
        (0, 0, 0),
        (1, 0, 0),
        (1, 1, 0),
        (0, 1, 0),
        (0, 0, 1),
        (1, 0, 1),
        (1, 1, 1),
        (0, 1, 1),
    ]
    faces = [
        (0, 1, 2, 3),
        (4, 5, 6, 7),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
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
        _header("boundary", "polyBoundaryMesh") + "2\n(\n"
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


def _convert_openfoam_case(
    case: Path, output: Path, env: dict[str, str]
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
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


def _solver_environment() -> dict[str, str]:
    env = os.environ.copy()
    package_root = str(Path(__file__).resolve().parents[2] / "src" / "cfdx" / "python")
    env["PYTHONPATH"] = os.pathsep.join(
        part for part in (package_root, env.get("PYTHONPATH", "")) if part
    )
    return env


def _run_production_solver(
    case_h5: Path, output: Path
) -> subprocess.CompletedProcess[str]:
    solver = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver:
        import pytest

        pytest.skip("CFDX_PRODUCTION_SOLVER is not configured for this test run")
    return subprocess.run(
        [
            solver,
            "--mesh",
            str(case_h5),
            "--output-dir",
            str(output),
            "--iterations",
            "1",
            "--initial-velocity-field",
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
        env=_solver_environment(),
    )


def _converted_case_with_velocity(
    tmp_path: Path, values: list[tuple[float, float, float]]
) -> Path:
    case = tmp_path / "foam_case"
    _write_hex_case(case)
    _write_velocity(case, "10", values)
    converted = tmp_path / "converted"
    completed = _convert_openfoam_case(case, converted, _solver_environment())
    assert completed.returncode == 0, completed.stdout + completed.stderr

    # The OpenFOAM reader intentionally imports mesh/results without guessing
    # CFDX numerical choices. Supply a valid, explicit CFDX selection for this
    # production-solver boundary test rather than weakening the solver contract.
    import json

    setup = {
        "schema_version": 1,
        "source": {
            "solver": "openfoam",
            "version": "test-fixture",
            "case_name": "openfoam_velocity_import",
            "format": "openfoam",
        },
        "physics_model": "incompressible_laminar",
        "turbulence_model": "laminar",
        "energy_model": "isothermal",
        "solver_mode": "steady",
        "units": "SI",
        "boundary_conditions": [
            {"patch_name": "inlet", "type": "inlet", "value_type": "fixed",
             "velocity_vector": [0.0, 0.0, 0.0]},
            {"patch_name": "outlet", "type": "outlet", "value_type": "zero_gradient"}
        ],
        "numerics": {
            "momentum_scheme": "first_order",
            "pressure_scheme": "standard",
            "momentum_interpolation": "linear",
            "transient_scheme": "steady",
            "gradient_operator": "green_gauss_cell",
            "coupled_solver": "SIMPLE",
            "max_iterations": 1,
            "selection": {
                "entries": [
                    {
                        "family": "gradient",
                        "configuration_key": "numerics.gradient.gauss",
                    },
                    {
                        "family": "convection",
                        "configuration_key": "numerics.convection.upwind",
                    },
                    {
                        "family": "pressure_velocity",
                        "configuration_key": "pressure_velocity.simple",
                    },
                    {
                        "family": "linear_solver",
                        "configuration_key": "linear.fgmres",
                    },
                ],
                "required_families": [
                    "gradient",
                    "convection",
                    "pressure_velocity",
                ],
            },
        },
    }
    case_h5 = converted / "case.cfdx.h5"
    with h5py.File(case_h5, "r+") as h5:
        h5.attrs["case_setup_json"] = json.dumps(setup, separators=(",", ":"))
    return case_h5


def test_production_solver_emits_exact_imported_velocity_evidence(
    tmp_path: Path,
) -> None:
    expected = [[4.25, -0.5, 8.0]]
    case_h5 = _converted_case_with_velocity(
        tmp_path, [(expected[0][0], expected[0][1], expected[0][2])]
    )
    output = tmp_path / "solver_output"
    completed = _run_production_solver(case_h5, output)

    assert completed.returncode == 0, (
        f"production solver failed (exit {completed.returncode})\\n"
        f"stdout:\\n{completed.stdout}\\nstderr:\\n{completed.stderr}"
    )
    import json

    imported = json.loads(
        (output / "imported_velocity_field.json").read_text(encoding="utf-8")
    )
    assert imported["field"] == "U"
    assert imported["location"] == "cell"
    assert imported["cell_count"] == 1
    assert imported["values"] == expected

    initialized_path = output / "initialized_velocity_field.json"
    assert initialized_path.is_file(), "solver did not emit post-initialization U"
    initialized = json.loads(initialized_path.read_text(encoding="utf-8"))
    assert initialized["field"] == "U"
    assert initialized["location"] == "cell"
    assert initialized["stage"] == "after_initialization_before_first_residual"
    assert initialized["cell_count"] == 1
    assert initialized["values"] == expected
    assert "Initialization strategy provided" in completed.stdout


def test_production_solver_rejects_missing_imported_velocity(tmp_path: Path) -> None:
    import h5py

    case_h5 = _converted_case_with_velocity(tmp_path, [(1.0, 2.0, 3.0)])
    with h5py.File(case_h5, "r+") as h5:
        del h5["fields/vector/U"]

    completed = _run_production_solver(case_h5, tmp_path / "solver_output")
    assert completed.returncode != 0
    assert "field 'U' is absent" in completed.stderr


def test_production_solver_rejects_nonfinite_imported_velocity(tmp_path: Path) -> None:
    import h5py
    import numpy as np

    case_h5 = _converted_case_with_velocity(tmp_path, [(1.0, 2.0, 3.0)])
    with h5py.File(case_h5, "r+") as h5:
        h5["fields/vector/U"][0, 1] = np.nan

    completed = _run_production_solver(case_h5, tmp_path / "solver_output")
    assert completed.returncode != 0
    assert "imported U contains non-finite values" in completed.stderr


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
    assert "Conversion BLOCKED due to blocking incompatibilities" in completed.stdout
