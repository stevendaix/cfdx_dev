from __future__ import annotations

import csv
import json
import math
import os
import subprocess
from pathlib import Path

import h5py
import numpy as np
import pytest

# Case-level body-force Poiseuille on the full `.cfdx.h5` production path.
#
# This closes the last N9 physical boundary documented in #663/#601: the
# Poiseuille case of tests/validation/test_n9_physical_matrix.cpp is driven
# by a uniform body force, but the production solver never consumed a
# case-declared body force (or material), so the physical case could not be
# expressed through the production `.cfdx.h5` path. The case written here is
# the same physical channel as the N9 matrix Poiseuille entry:
#
#   * 12x16 hexahedral single-layer channel, [0,1]^3
#   * inlet/outlet/front/back: zero-gradient velocity (front/back are EMPTY
#     patches, so the layer is the same 2-D plane channel as the matrix)
#   * bottom/top: no-slip fixed walls, u = 0
#   * pressure zero-gradient everywhere
#   * body_force = (1, 0, 0) per unit mass, rho = 1.0, mu = 0.1 -> nu = 0.1
#   * analytic solution u(y) = G*y*(1-y)/(2*nu), same L2 <= 5e-3 gate as the
#     N9 physical matrix.
#
# The mesh topology mirrors make_channel_mesh(12, 16, 0.0) from the matrix
# test and is written directly in the CFDX-HDF5 mesh-v1 interchange format
# (the Python meshio bridge cannot express named channel patches yet).

NX, NY = 12, 16
BODY_FORCE = (1.0, 0.0, 0.0)  # per unit mass (m/s^2), x-driven
DENSITY = 1.0  # kg/m^3
DYNAMIC_VISCOSITY = 0.1  # Pa s -> nu = 0.1 m^2/s
L2_GATE = 5e-3  # identical to the N9 matrix Poiseuille analytic gate
MAX_ITERATIONS = 3000

_PATCH_TYPE = {"inlet": 1, "outlet": 2, "bottom": 0, "top": 0, "front": 6, "back": 6}
_PATCH_ORDER = ("inlet", "outlet", "bottom", "top", "front", "back")

_FNV_OFFSET_BASIS = 1469598103934665603
_FNV_PRIME = 1099511628211
_UINT64_MASK = 0xFFFFFFFFFFFFFFFF


def _fnv1a_update(hash_value: int, values: np.ndarray) -> int:
    arr = np.ascontiguousarray(values)
    if arr.size == 0:
        return hash_value
    for byte in arr.tobytes(order="C"):
        hash_value = ((hash_value ^ byte) * _FNV_PRIME) & _UINT64_MASK
    return hash_value


def _hash_hex(hash_value: int) -> str:
    return f"{hash_value:016x}"


def _integrity_hashes(m: dict) -> tuple[str, str, str]:
    """Mirror read_case_cfdx_h5()'s FNV-1a digest order byte-for-byte."""
    topology = _FNV_OFFSET_BASIS
    topology = _fnv1a_update(topology, m["face_vertices"])
    if m["n_faces"] > 0:
        topology = _fnv1a_update(topology, m["face_offsets"])
    topology = _fnv1a_update(topology, m["owner"])
    topology = _fnv1a_update(topology, m["neighbour"])
    topology = _fnv1a_update(topology, m["cell_faces"])
    if m["n_cells"] > 0:
        topology = _fnv1a_update(topology, m["cell_offsets"])

    points = m["points"]
    coords = [np.ascontiguousarray(points[:, i], dtype=np.float64) for i in range(3)]
    geometry = _FNV_OFFSET_BASIS
    for c in coords:
        geometry = _fnv1a_update(geometry, c)
    mesh_hash = topology
    for c in coords:
        mesh_hash = _fnv1a_update(mesh_hash, c)
    return _hash_hex(topology), _hash_hex(geometry), _hash_hex(mesh_hash)


def _build_channel() -> dict:
    """Hexahedral single-layer channel, mirroring make_channel_mesh(NX, NY, 0)."""
    nx, ny = NX, NY
    plane = (nx + 1) * (ny + 1)
    points = np.zeros((2 * plane, 3), dtype=np.float64)

    def pid(i: int, j: int, k: int) -> int:
        return (j * (nx + 1) + i) * 2 + k

    for j in range(ny + 1):
        y = j / ny
        for i in range(nx + 1):
            x = i / nx
            points[pid(i, j, 0)] = (x, y, 0.0)
            points[pid(i, j, 1)] = (x, y, 1.0)

    face_map: dict[tuple[int, ...], int] = {}
    faces: list[list[int]] = []
    owner: list[int] = []
    neighbour: list[int] = []
    cell_faces: list[list[int]] = []

    def add_face(vertices: tuple[int, ...], cell: int) -> int:
        key = tuple(sorted(vertices))
        existing = face_map.get(key)
        if existing is not None:
            neighbour[existing] = cell
            return existing
        fid = len(faces)
        faces.append(list(vertices))
        owner.append(cell)
        neighbour.append(-1)
        face_map[key] = fid
        return fid

    for j in range(ny):
        for i in range(nx):
            cell = j * nx + i
            a, b = pid(i, j, 0), pid(i + 1, j, 0)
            c0, d = pid(i + 1, j + 1, 0), pid(i, j + 1, 0)
            e, f = pid(i, j, 1), pid(i + 1, j, 1)
            g, h = pid(i + 1, j + 1, 1), pid(i, j + 1, 1)
            cell_faces.append(
                [
                    add_face((a, d, c0, b), cell),
                    add_face((e, f, g, h), cell),
                    add_face((a, b, f, e), cell),
                    add_face((d, h, g, c0), cell),
                    add_face((a, e, h, d), cell),
                    add_face((b, c0, g, f), cell),
                ]
            )

    # Patch assignment is geometric, exactly like the matrix mesh builder.
    patch_ids: dict[str, list[int]] = {name: [] for name in _PATCH_ORDER}
    eps = 1e-12
    for fid, face in enumerate(faces):
        if neighbour[fid] >= 0:
            continue
        x, y, z = points[face].mean(axis=0)
        if abs(y) < eps:
            name = "bottom"
        elif abs(y - 1.0) < eps:
            name = "top"
        elif abs(z) < eps:
            name = "front"
        elif abs(z - 1.0) < eps:
            name = "back"
        elif abs(x) < eps:
            name = "inlet"
        elif abs(x - 1.0) < eps:
            name = "outlet"
        else:  # pragma: no cover - geometry is exact
            raise AssertionError("unclassified channel boundary face")
        patch_ids[name].append(fid)

    return {
        "points": points,
        "faces": faces,
        "owner": owner,
        "neighbour": neighbour,
        "cell_faces": cell_faces,
        "patch_ids": patch_ids,
        "n_points": len(points),
        "n_faces": len(faces),
        "n_cells": len(cell_faces),
    }


def _case_setup_json() -> str:
    from cfdx.io.schema import (
        BCType,
        BCValueType,
        BoundarySpec,
        CaseSetup,
        InitialCondition,
        MaterialSpec,
        NumericalSelection,
        NumericalSelectionConfig,
        SourceInfo,
    )

    setup = CaseSetup(
        source=SourceInfo(
            solver="cfdx",
            version="n9-body-force",
            case_name="n9_poiseuille_body_force",
            format="cfdx",
        )
    )
    setup.physics_model = "incompressible_laminar"
    setup.turbulence_model = "laminar"
    setup.energy_model = "isothermal"
    setup.solver_mode = "steady"
    setup.units = "SI"

    setup.materials = [
        MaterialSpec(name="fluid", density=DENSITY, dynamic_viscosity=DYNAMIC_VISCOSITY)
    ]
    # The physics closure under test: a case-level uniform body force.
    setup.body_force = list(BODY_FORCE)

    zero = [0.0, 0.0, 0.0]
    setup.boundary_conditions = [
        BoundarySpec(
            patch_name="inlet",
            type=BCType.INLET,
            value_type=BCValueType.ZERO_GRADIENT,
        ),
        BoundarySpec(
            patch_name="outlet",
            type=BCType.OUTLET,
            value_type=BCValueType.ZERO_GRADIENT,
        ),
        BoundarySpec(
            patch_name="bottom",
            type=BCType.WALL,
            value_type=BCValueType.FIXED,
            velocity_vector=list(zero),
        ),
        BoundarySpec(
            patch_name="top",
            type=BCType.WALL,
            value_type=BCValueType.FIXED,
            velocity_vector=list(zero),
        ),
        BoundarySpec(
            patch_name="front",
            type=BCType.EMPTY,
            value_type=BCValueType.ZERO_GRADIENT,
        ),
        BoundarySpec(
            patch_name="back",
            type=BCType.EMPTY,
            value_type=BCValueType.ZERO_GRADIENT,
        ),
    ]

    setup.initial_condition = InitialCondition(
        velocity=0.0, pressure=0.0, velocity_vector=list(zero)
    )

    setup.numerics.max_iterations = MAX_ITERATIONS
    setup.numerics.selection = NumericalSelectionConfig(
        entries=[
            NumericalSelection(
                family="gradient",
                configuration_key="numerics.gradient.gauss",
            ),
            NumericalSelection(
                family="convection",
                configuration_key="numerics.convection.upwind",
            ),
            NumericalSelection(
                family="pressure_velocity",
                configuration_key="pressure_velocity.simple",
            ),
            NumericalSelection(
                family="linear_solver",
                configuration_key="linear.fgmres",
            ),
        ],
        required_families=["gradient", "convection", "pressure_velocity"],
    )

    return json.dumps(setup.model_dump(mode="json"))


def _write_case(path: Path) -> None:
    mesh = _build_channel()
    faces = mesh["faces"]
    cell_faces = mesh["cell_faces"]

    fv = np.asarray([v for f in faces for v in f], dtype=np.uint64)
    fo = np.zeros(len(faces) + 1, dtype=np.uint64)
    for i, f in enumerate(faces):
        fo[i + 1] = fo[i] + len(f)
    own = np.asarray(mesh["owner"], dtype=np.uint64)
    nei = np.asarray(mesh["neighbour"], dtype=np.int64)
    cf = np.asarray([f for c in cell_faces for f in c], dtype=np.uint64)
    co = np.zeros(len(cell_faces) + 1, dtype=np.uint64)
    for i, c in enumerate(cell_faces):
        co[i + 1] = co[i] + len(c)

    patch_face_ids: list[int] = []
    patch_face_offsets = [0]
    patch_meta: list[str] = []
    for name in _PATCH_ORDER:
        ids = mesh["patch_ids"][name]
        assert ids, f"patch {name} has no faces"
        patch_face_ids.extend(ids)
        patch_face_offsets.append(len(patch_face_ids))
        patch_meta.append(f"{name}:{ids[0]}:{len(ids)}:{_PATCH_TYPE[name]}")

    topology_hash, geometry_hash, mesh_hash = _integrity_hashes(
        {
            "points": mesh["points"],
            "face_vertices": fv,
            "face_offsets": fo,
            "owner": own,
            "neighbour": nei,
            "cell_faces": cf,
            "cell_offsets": co,
            "n_faces": mesh["n_faces"],
            "n_cells": mesh["n_cells"],
        }
    )

    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h5.attrs["format_version"] = "1"
        h5.attrs["schema_version"] = "1"
        h5.attrs["cfdx_version"] = "0.7"
        h5.attrs["mesh_topology"] = "cfdx-csr-v1"
        h5.attrs["n_points"] = str(mesh["n_points"])
        h5.attrs["n_faces"] = str(mesh["n_faces"])
        h5.attrs["n_cells"] = str(mesh["n_cells"])
        h5.attrs["topology_hash"] = topology_hash
        h5.attrs["geometry_hash"] = geometry_hash
        h5.attrs["mesh_hash"] = mesh_hash
        h5.attrs["boundary_patches"] = ";".join(patch_meta)
        h5.attrs["source_solver"] = "cfdx"
        h5.attrs["source_format"] = "cfdx"
        h5.attrs["source_version"] = "n9-body-force"
        h5.attrs["source_case_path"] = ""
        h5.attrs["source_case_name"] = "n9_poiseuille_body_force"
        h5.attrs["case_setup_json"] = _case_setup_json()

        h5.create_dataset("points", data=mesh["points"])
        h5.create_dataset("face_vertices", data=fv)
        h5.create_dataset("face_offsets", data=fo)
        h5.create_dataset("owner", data=own)
        h5.create_dataset("neighbour", data=nei)
        h5.create_dataset("cell_faces", data=cf)
        h5.create_dataset("cell_offsets", data=co)
        h5.create_dataset(
            "patch_face_ids", data=np.asarray(patch_face_ids, dtype=np.uint64)
        )
        h5.create_dataset(
            "patch_face_offsets", data=np.asarray(patch_face_offsets, dtype=np.uint64)
        )


def test_n9_poiseuille_body_force_production_path(tmp_path: Path) -> None:
    solver_path = os.environ.get("CFDX_PRODUCTION_SOLVER")
    if not solver_path:
        pytest.skip("CFDX_PRODUCTION_SOLVER is provided by ctest")
    solver = Path(solver_path)
    assert solver.is_file()

    case_path = tmp_path / "n9_poiseuille_body_force.cfdx.h5"
    _write_case(case_path)

    # The serialized case must carry the declared body force through the
    # canonical schema: this is the boundary the test exists to close.
    with h5py.File(case_path, "r") as h5:
        setup = json.loads(h5.attrs["case_setup_json"])
        assert setup["body_force"] == list(BODY_FORCE)
        entries = {
            (entry["family"], entry["configuration_key"])
            for entry in setup["numerics"]["selection"]["entries"]
        }
        assert ("pressure_velocity", "pressure_velocity.simple") in entries

    # Centreline column probes: the plane-channel solution is x-invariant, so
    # the mid-x column samples the whole field. One u_x probe per y row,
    # located exactly at the cell centres of that column.
    probes: list[str] = []
    n_cells_x = NX
    n_cells_y = NY
    x = (n_cells_x // 2 + 0.5) / n_cells_x
    for j in range(n_cells_y):
        y = (j + 0.5) / n_cells_y
        probes.append(f"u_y{j:02d}:{x!r},{y!r},0.5:u_x")

    output_dir = tmp_path / "run"
    probe_csv = tmp_path / "probes.csv"
    completed = subprocess.run(
        [
            str(solver),
            "--mesh",
            str(case_path),
            "--output-dir",
            str(output_dir),
            "--iterations",
            str(MAX_ITERATIONS),
            *[item for probe in probes for item in ("--probe", probe)],
            "--probe-csv",
            str(probe_csv),
        ],
        capture_output=True,
        text=True,
        timeout=800,
        check=False,
    )
    diagnostics = completed.stdout + completed.stderr
    assert completed.returncode == 0, diagnostics[-12000:]
    assert "Resolved numerical selections:" in diagnostics
    assert "scheme[pressure_velocity]=pressure_velocity.simple" in diagnostics
    assert "Converged YES" in diagnostics

    # Final sample per probe: rows are ordered by probe name then iteration.
    final: dict[str, float] = {}
    with probe_csv.open(newline="") as stream:
        for row in csv.reader(stream):
            if not row or row[0].startswith("#"):
                continue
            name, _iteration, _time, value = row
            final[name] = float(value)  # rows are iteration-ordered per probe

    assert len(final) == n_cells_y, diagnostics[-12000:]

    nu = DYNAMIC_VISCOSITY / DENSITY
    g = BODY_FORCE[0]
    squared = 0.0
    for j in range(n_cells_y):
        y = (j + 0.5) / n_cells_y
        exact = g * y * (1.0 - y) / (2.0 * nu)
        error = final[f"u_y{j:02d}"] - exact
        squared += error * error
    l2 = math.sqrt(squared / n_cells_y)
    assert l2 <= L2_GATE, f"body-force Poiseuille analytic L2 gate failed: {l2:.6e}"
