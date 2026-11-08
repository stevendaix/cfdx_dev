"""Production CLI regression for persisted SST cases and DAT restart equivalence.

The fixture starts from a real converted CFDX case so the C++ case reader,
boundary metadata reader, and production executable are exercised together.
This is software-integration evidence, not a physical SST qualification.
"""

from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest
from cfdx.dat_io import read_dat_restart


def _prepare_persisted_sst_case(case_path: Path) -> None:
    """Declare explicit SST initial data and patch-level scalar BCs in a case."""
    with h5py.File(case_path, "r+") as h5:
        raw = h5.attrs["case_setup_json"]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        setup = json.loads(str(raw))
        setup["physics_model"] = "incompressible_laminar"
        setup["turbulence_model"] = "k_omega_sst"
        setup["energy_model"] = "none"
        setup["radiation_model"] = "none"
        setup["multiphase_model"] = "none"
        setup["transient"] = False
        setup["solver_mode"] = "steady"
        setup["materials"] = [
            {
                "name": "air",
                "density": 1.2,
                "dynamic_viscosity": 0.012,
                "turbulence_model": "k_omega_sst",
            }
        ]
        setup["initial_condition"] = {
            "velocity": 1.0,
            "pressure": 0.0,
            "velocity_vector": [1.0, 0.0, 0.0],
            "scalar_fields": {"k": 0.01, "omega": 2.0},
        }
        # The source SU2 fixture is an implicit-time case. Remove only its
        # temporal selection because this executable is steady; retain the
        # converter's other resolved N1 selections and their validation gate.
        for key in ("numerics", "numerical_config"):
            numerics = setup.get(key)
            if not isinstance(numerics, dict):
                continue
            selection = numerics.get("selection")
            if not isinstance(selection, dict):
                continue
            entries = selection.get("entries", [])
            selection["entries"] = [
                entry for entry in entries if entry.get("family") != "temporal"
            ]
            required = selection.get("required_families")
            if isinstance(required, list):
                selection["required_families"] = [
                    family for family in required if family != "temporal"
                ]

            # The SU2 adapter records JST as unmapped; the production driver
            # intentionally rejects that conversion gap. This integration
            # fixture explicitly chooses a supported convection method rather
            # than weakening the driver's numerical-selection contract.
            if key == "numerics":
                entries = selection.setdefault("entries", [])
                if not any(entry.get("family") == "convection" for entry in entries):
                    entries.append(
                        {
                            "family": "convection",
                            "configuration_key": "numerics.convection.upwind",
                        }
                    )
                required = selection.setdefault("required_families", [])
                if "convection" not in required:
                    required.append("convection")
        # The steady production executable rejects temporal declarations from
        # both the legacy scheme and the resolved numerical-method report.
        numerics = setup.get("numerics")
        if isinstance(numerics, dict):
            numerics["transient_scheme"] = "steady"
        numerical_config = setup.get("numerical_config")
        if isinstance(numerical_config, dict):
            numerical_config["transient_scheme"] = "steady"
            entries = numerical_config.get("entries")
            if isinstance(entries, list):
                numerical_config["entries"] = [
                    entry for entry in entries if entry.get("family") != "temporal"
                ]
        numerical_report = setup.get("numerical_report")
        if isinstance(numerical_report, dict):
            resolved = numerical_report.get("resolved")
            if isinstance(resolved, list):
                numerical_report["resolved"] = [
                    entry for entry in resolved if entry.get("family") != "temporal"
                ]
        for boundary in setup.get("boundary_conditions", []):
            patch_type = str(boundary.get("type", "")).lower()
            value_type = str(boundary.get("value_type", "")).lower()
            if str(boundary.get("patch_name", "")).lower() == "outlet":
                boundary["type"] = "pressure_outlet"
                boundary["value_type"] = "outlet_pressure"
                boundary["pressure"] = 0.0
            if patch_type == "inlet":
                # The production SST path derives inlet k/omega from these
                # explicit, dimensionally consistent inputs. Defaults of zero
                # are intentionally rejected by the case contract.
                boundary["velocity_magnitude"] = 1.0
                boundary["velocity_vector"] = [1.0, 0.0, 0.0]
                boundary["turbulence_intensity"] = 0.05
                boundary["turbulence_length_scale"] = 0.1
            if patch_type == "wall" or value_type in {"wall_no_slip", "wall"}:
                boundary["sst_scalar_conditions"] = [
                    {"field": "k", "type": "wall_k_zero", "value": 0.0},
                    {"field": "omega", "type": "wall_omega_sst", "value": 0.0},
                ]
            else:
                boundary["sst_scalar_conditions"] = [
                    {"field": "k", "type": "fixed_value", "value": 0.01},
                    {"field": "omega", "type": "fixed_value", "value": 2.0},
                ]
        # Recreate the attribute because the serialized JSON may grow beyond
        # the converter's original HDF5 string-attribute width.
        del h5.attrs["case_setup_json"]
        h5.attrs["case_setup_json"] = json.dumps(setup)


def _make_volumetric_sst_restart_fixture(
    converted_case: Path, output_case: Path
) -> None:
    """Build a small *closed, three-dimensional* channel for SST lifecycle tests.

    The two-triangle SU2 converter test above checks metadata persistence, not
    a trustworthy 3D RANS control volume. The production/restart test instead
    uses 12 conforming hexahedral volumes (3 x 2 x 2) with no missing faces,
    two connected internal planes and all six physical patch categories.
    """
    import copy

    from cfdx.io.hdf5_writer import _boundary_patch_metadata, _integrity_hashes

    nx, ny, nz = 3, 2, 2

    def vertex(i: int, j: int, k: int) -> int:
        return (i * (ny + 1) + j) * (nz + 1) + k

    points = np.array(
        [
            [float(i), j / ny, k / nz]
            for i in range(nx + 1)
            for j in range(ny + 1)
            for k in range(nz + 1)
        ],
        dtype=np.float64,
    )
    faces: list[tuple[int, ...]] = []
    owner: list[int] = []
    neighbour: list[int] = []
    face_ids: dict[tuple[int, ...], int] = {}
    cells: list[list[int]] = []

    def add_face(vertices: tuple[int, ...], cell: int) -> int:
        key = tuple(sorted(vertices))
        if key in face_ids:
            face = face_ids[key]
            assert neighbour[face] == -1, "non-manifold regression fixture"
            neighbour[face] = cell
            return face
        face = len(faces)
        face_ids[key] = face
        faces.append(vertices)
        owner.append(cell)
        neighbour.append(-1)
        return face

    for i in range(nx):
        for j in range(ny):
            for k in range(nz):
                a, b = vertex(i, j, k), vertex(i + 1, j, k)
                c, d = vertex(i + 1, j + 1, k), vertex(i, j + 1, k)
                e, f = vertex(i, j, k + 1), vertex(i + 1, j, k + 1)
                g, h = vertex(i + 1, j + 1, k + 1), vertex(i, j + 1, k + 1)
                cell = len(cells)
                cells.append(
                    [
                        add_face((a, d, c, b), cell),
                        add_face((e, f, g, h), cell),
                        add_face((a, b, f, e), cell),
                        add_face((d, h, g, c), cell),
                        add_face((a, e, h, d), cell),
                        add_face((b, c, g, f), cell),
                    ]
                )

    patches = {
        name: []
        for name in (
            "inlet",
            "outlet",
            "wall_bottom",
            "wall_top",
            "wall_front",
            "wall_back",
        )
    }
    for face, vertices in enumerate(faces):
        if neighbour[face] >= 0:
            continue
        x, y, z = points[list(vertices)].mean(axis=0)
        if abs(x) < 1.0e-12:
            patch = "inlet"
        elif abs(x - nx) < 1.0e-12:
            patch = "outlet"
        elif abs(y) < 1.0e-12:
            patch = "wall_bottom"
        elif abs(y - 1) < 1.0e-12:
            patch = "wall_top"
        elif abs(z) < 1.0e-12:
            patch = "wall_front"
        elif abs(z - 1) < 1.0e-12:
            patch = "wall_back"
        else:
            raise AssertionError("unclassified volumetric SST boundary")
        patches[patch].append(face)
    assert all(patches.values())
    assert sum(map(len, patches.values())) == sum(n < 0 for n in neighbour)
    assert sum(n >= 0 for n in neighbour) > 0

    offsets = np.zeros(len(faces) + 1, dtype=np.uint64)
    cell_offsets = np.zeros(len(cells) + 1, dtype=np.uint64)
    for face, vertices in enumerate(faces):
        offsets[face + 1] = offsets[face] + len(vertices)
    for cell, vertices in enumerate(cells):
        cell_offsets[cell + 1] = cell_offsets[cell] + len(vertices)
    patch_kinds = {"inlet": "inlet", "outlet": "outlet"}
    mesh = {
        "points": points,
        "face_vertices": np.array([v for f in faces for v in f], dtype=np.uint64),
        "face_offsets": offsets,
        "owner": np.array(owner, dtype=np.uint64),
        "neighbour": np.array(neighbour, dtype=np.int64),
        "cell_faces": np.array([v for c in cells for v in c], dtype=np.uint64),
        "cell_offsets": cell_offsets,
        "n_points": len(points),
        "n_faces": len(faces),
        "n_cells": len(cells),
        "boundary_patch_assignment": [
            {
                "name": name,
                "type": patch_kinds.get(name, "wall"),
                "face_ids": patch_faces,
            }
            for name, patch_faces in patches.items()
        ],
    }
    with h5py.File(converted_case, "r") as h:
        raw = h.attrs["case_setup_json"]
        setup = json.loads(raw.decode() if isinstance(raw, bytes) else str(raw))
    original = {b["patch_name"]: b for b in setup["boundary_conditions"]}
    boundaries = []
    for name in patches:
        if name in original:
            boundaries.append(copy.deepcopy(original[name]))
        else:
            template = copy.deepcopy(original["wall_top"])
            template["patch_name"] = name
            boundaries.append(template)
    setup["boundary_conditions"] = boundaries
    # A generated 3D fixture is not a SU2 2D import. Record that distinction.
    source_metadata = setup.setdefault("source_metadata", {})
    source_metadata["fixture_kind"] = "volumetric_3d_SST_restart"
    topology_hash, geometry_hash, mesh_hash = _integrity_hashes(mesh)
    metadata, patch_face_ids, patch_face_offsets = _boundary_patch_metadata(mesh)
    with h5py.File(output_case, "w") as h:
        for key, value in {
            "format": "CFDX-HDF5-mesh-v1",
            "format_version": "1",
            "schema_version": "1",
            "cfdx_version": "0.7",
            "mesh_topology": "cfdx-csr-v1",
            "n_points": str(len(points)),
            "n_faces": str(len(faces)),
            "n_cells": str(len(cells)),
            "topology_hash": topology_hash,
            "geometry_hash": geometry_hash,
            "mesh_hash": mesh_hash,
            "boundary_patches": metadata,
            "source_solver": "CFDX test fixture",
            "source_format": "generated_3d",
            "source_version": "1",
            "source_case_path": str(output_case),
            "source_case_name": "sst_volumetric_restart",
            "case_setup_json": json.dumps(setup),
        }.items():
            h.attrs[key] = value
        for name in (
            "points",
            "face_vertices",
            "face_offsets",
            "owner",
            "neighbour",
            "cell_faces",
            "cell_offsets",
        ):
            h.create_dataset(name, data=mesh[name])
        h.create_dataset("patch_face_ids", data=patch_face_ids)
        h.create_dataset("patch_face_offsets", data=patch_face_offsets)


def _run(
    solver: Path,
    case: Path,
    output: Path,
    iterations: int,
    restart: Path | None = None,
) -> subprocess.CompletedProcess[str]:
    command = [
        str(solver),
        "--mesh",
        str(case),
        "--output-dir",
        str(output),
        "--iterations",
        str(iterations),
        "--experimental-sst",
    ]
    if restart is not None:
        command.extend(["--restart", str(restart)])
    return subprocess.run(
        command, capture_output=True, text=True, timeout=120, check=False
    )


def _assert_close_field(left, right, name: str) -> None:
    assert left.dimension == right.dimension, f"{name}: dimension changed"
    a = np.asarray(left.values, dtype=float)
    b = np.asarray(right.values, dtype=float)
    assert a.shape == b.shape, f"{name}: shape changed ({a.shape} != {b.shape})"
    assert np.isfinite(a).all() and np.isfinite(b).all(), f"{name}: non-finite state"
    # The restart is a serialization boundary, not a new discretization.
    np.testing.assert_allclose(a, b, rtol=2e-10, atol=2e-12, err_msg=name)


def _assert_completed_run(
    run: subprocess.CompletedProcess[str], output: Path, iterations: int
) -> None:
    """Assert the documented CLI execution contract for a nominal run.

    The production CLI separates execution success from convergence
    (``execution_summary.h``): 0 = converged, 1 = completed without
    converging, 2 = execution error. This regression exercises persisted
    metadata and restart determinism, not convergence, so both documented
    completion codes are nominal here; an execution error is never accepted.
    The execution summary is cross-checked so a silent failure cannot pass
    as a completed budget run.
    """
    assert run.returncode in (0, 1), (run.stdout + run.stderr)[-16000:]
    summary_path = output / "execution.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    assert summary["process_exit_code"] == run.returncode
    # The fixture premise is an exact budget run: each trajectory must have
    # performed its full iteration budget so the field comparison below
    # compares equivalent iteration counts.
    assert summary["iterations"] == iterations
    assert summary["artifacts"]["restart_dat"] is True
    assert summary["artifacts"]["convergence_json"] is True
    expected_converged = run.returncode == 0
    assert summary["converged"] is expected_converged


def test_persisted_sst_cli_wall_omega_and_dat_restart_equivalence(
    tmp_path: Path,
) -> None:
    solver_text = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver_text:
        pytest.skip("CFDX_PRODUCTION_SOLVER is provided by CTest")
    solver = Path(solver_text)
    assert solver.is_file()

    # Keep this integration fixture deliberately tiny.  The regression is
    # about persisted SST metadata and restart determinism, not airfoil physics;
    # the NACA mesh made each first SST iteration exceed the CI test budget.
    case_dir = tmp_path / "persisted_case"
    case_dir.mkdir()
    source_case = case_dir / "sst_restart_fixture.su2"
    config_case = case_dir / "sst_restart_fixture.cfg"
    source_case.write_text(
        """NDIME= 2
NELEM= 2
5 0 1 2 0
5 0 2 3 1
NPOIN= 4
0.0 0.0 0
1.0 0.0 1
1.0 1.0 2
0.0 1.0 3
NMARK= 4
MARKER_TAG= inlet
MARKER_ELEMS= 1
3 0 3
MARKER_TAG= outlet
MARKER_ELEMS= 1
3 1 2
MARKER_TAG= wall_bottom
MARKER_ELEMS= 1
3 0 1
MARKER_TAG= wall_top
MARKER_ELEMS= 1
3 2 3
""",
        encoding="utf-8",
    )
    config_case.write_text(
        """SOLVER= INC_RANS
KIND_TURB_MODEL= SST
INC_DENSITY_INIT= 1.2
INC_VELOCITY_INIT= (1.0, 0.0, 0.0)
MU_CONSTANT= 0.012
CONV_NUM_METHOD_FLOW= ROE
MUSCL_FLOW= NO
NUM_METHOD_GRAD= GREEN_GAUSS
TIME_DISCRE_FLOW= EULER_IMPLICIT
ITER= 2
""",
        encoding="utf-8",
    )
    case_path = case_dir / "case.cfdx.h5"

    from cfdx.io.converter import convert

    converted = convert(source_case, output=case_path)
    assert converted.case is not None
    _prepare_persisted_sst_case(case_path)
    # Prove that the contract is stored in the HDF5 file, not only in memory.
    with h5py.File(case_path, "r") as h5:
        raw = h5.attrs["case_setup_json"]
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        setup = json.loads(str(raw))
        assert setup["turbulence_model"] == "k_omega_sst"
        assert any(
            condition.get("type") == "wall_omega_sst"
            for boundary in setup["boundary_conditions"]
            for condition in boundary.get("sst_scalar_conditions", [])
        ), "fixture must persist at least one WALL_OMEGA_SST obligation"
        inlets = [
            boundary
            for boundary in setup["boundary_conditions"]
            if str(boundary.get("type", "")).lower() == "inlet"
        ]
        assert inlets, "fixture must persist an inlet boundary"
        for inlet in inlets:
            assert 0.0 < inlet["turbulence_intensity"] < 1.0
            assert inlet["turbulence_length_scale"] > 0.0

    # The tiny two-triangle converter fixture is retained above as a metadata
    # regression. It is not an admissible volumetric SST transport benchmark.
    # Use a connected 3D finite-volume channel for deterministic restart.
    production_case = case_dir / "sst_volumetric_restart.cfdx.h5"
    _make_volumetric_sst_restart_fixture(case_path, production_case)
    with h5py.File(production_case, "r") as h:
        assert len(h["cell_offsets"]) - 1 == 12
        assert (np.asarray(h["neighbour"]) >= 0).any()

    continuous = _run(solver, production_case, tmp_path / "continuous", 2)
    _assert_completed_run(continuous, tmp_path / "continuous", 2)
    continuous_dat = tmp_path / "continuous" / "restart.dat"
    assert continuous_dat.is_file()
    continuous_state = read_dat_restart(continuous_dat)
    assert {"U", "p", "k", "omega"}.issubset(continuous_state.fields)

    first = _run(solver, production_case, tmp_path / "first_segment", 1)
    _assert_completed_run(first, tmp_path / "first_segment", 1)
    first_dat = tmp_path / "first_segment" / "restart.dat"
    assert first_dat.is_file()
    first_state = read_dat_restart(first_dat)
    assert {"U", "p", "k", "omega"}.issubset(first_state.fields)

    # --iterations is a per-invocation iteration budget, not a global target:
    # the resumed segment performs one more iteration after the first segment.
    resumed = _run(
        solver, production_case, tmp_path / "resumed_segment", 1, restart=first_dat
    )
    _assert_completed_run(resumed, tmp_path / "resumed_segment", 1)
    resumed_dat = tmp_path / "resumed_segment" / "restart.dat"
    assert resumed_dat.is_file()
    resumed_state = read_dat_restart(resumed_dat)
    assert resumed_state.cells == continuous_state.cells
    for name in ("U", "p", "k", "omega"):
        _assert_close_field(
            continuous_state.fields[name], resumed_state.fields[name], name
        )
