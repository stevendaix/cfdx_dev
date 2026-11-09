"""Generate an auditable CFDX HDF5 *debug* backward-facing step mesh.

A single-layer hex extrusion with EMPTY span planes and explicit step/no-slip
patches. This is for production SST integration and wall-shear debugging, not
a Driver-Seegmiller validation mesh: the NASA reference requires measured
inlet boundary layer, documented farfield geometry, wall y+, and 3-grid
convergence. No synthetic empirical observations are produced.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src/cfdx/python"))

import h5py
import numpy as np
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
    SSTScalarBoundarySpec,
    SSTScalarBoundaryType,
)

FNV_OFFSET = 1469598103934665603
FNV_PRIME = 1099511628211
PATCH_TYPES = {
    "inlet": 1,
    "outlet": 2,
    "upstream_wall": 0,
    "step_wall": 0,
    "downstream_wall": 0,
    "top": 0,
    "front": 6,
    "back": 6,
}


def fnv(hash_value: int, array: np.ndarray) -> int:
    for byte in np.ascontiguousarray(array).tobytes():
        hash_value = ((hash_value ^ byte) * FNV_PRIME) & ((1 << 64) - 1)
    return hash_value


def make_mesh(
    n_up: int,
    n_down: int,
    ny: int,
    step_cells: int,
    height: float = 0.0127,
    wall_stretch: float = 0.0,
) -> dict:
    if not (
        2 <= n_up <= 500
        and 4 <= n_down <= 1500
        and 4 <= ny <= 500
        and 1 <= step_cells < ny // 2
        and height > 0
        and np.isfinite(wall_stretch)
        and 0 <= wall_stretch <= 5.0
    ):
        raise ValueError("invalid step resolution or wall stretching")
    if wall_stretch and min(step_cells, ny - step_cells) < 3:
        raise ValueError("wall stretching requires at least 3 cells in each y segment")
    xs = np.concatenate(
        (
            np.linspace(-6 * height, 0, n_up + 1),
            np.linspace(0, 40 * height, n_down + 1)[1:],
        )
    )

    # Two-sided tanh grading keeps the exact step grid line, while targeting
    # y=0, y=h (including the upstream wall) and the opposite wall y=8h.
    # This is a *debug* grid improvement; y+ must still be checked with CFDX.
    def cluster(steps: int) -> np.ndarray:
        s = np.linspace(0.0, 1.0, steps + 1)
        if not wall_stretch:
            return s
        return 0.5 * (1 + np.tanh(wall_stretch * (2 * s - 1)) / np.tanh(wall_stretch))

    low = height * cluster(step_cells)
    high = height + 7 * height * cluster(ny - step_cells)
    ys = np.concatenate((low, high[1:]))
    if np.min(np.diff(ys)) <= 1e-12 * height:
        raise ValueError("wall stretching produces degenerate y cells")
    nx = n_up + n_down
    points = np.empty(((nx + 1) * (ny + 1) * 2, 3), dtype=np.float64)

    def pid(i, j, k):
        return (j * (nx + 1) + i) * 2 + k

    for j, y in enumerate(ys):
        for i, x in enumerate(xs):
            points[pid(i, j, 0)] = (x, y, 0.0)
            points[pid(i, j, 1)] = (x, y, height)
    faces: list[list[int]] = []
    owners: list[int] = []
    neighbours: list[int] = []
    cell_faces: list[list[int]] = []
    fmap: dict[tuple[int, ...], int] = {}

    def add_face(verts, cell):
        key = tuple(sorted(verts))
        if key in fmap:
            fid = fmap[key]
            if neighbours[fid] >= 0:
                raise AssertionError("nonmanifold face")
            neighbours[fid] = cell
            return fid
        fid = len(faces)
        fmap[key] = fid
        faces.append(list(verts))
        owners.append(cell)
        neighbours.append(-1)
        return fid

    for j in range(ny):
        for i in range(nx):
            if i < n_up and j < step_cells:
                continue
            c = len(cell_faces)
            a, b = pid(i, j, 0), pid(i + 1, j, 0)
            c0, d = pid(i + 1, j + 1, 0), pid(i, j + 1, 0)
            e, f = pid(i, j, 1), pid(i + 1, j, 1)
            g, h = pid(i + 1, j + 1, 1), pid(i, j + 1, 1)
            cell_faces.append(
                [
                    add_face((a, d, c0, b), c),
                    add_face((e, f, g, h), c),
                    add_face((a, b, f, e), c),
                    add_face((d, h, g, c0), c),
                    add_face((a, e, h, d), c),
                    add_face((b, c0, g, f), c),
                ]
            )
    patches = {p: [] for p in PATCH_TYPES}
    tolerance = 1e-9 * height
    for fid, verts in enumerate(faces):
        if neighbours[fid] >= 0:
            continue
        xyz = points[verts].mean(axis=0)
        x, y, z = xyz
        if abs(z) < tolerance:
            name = "front"
        elif abs(z - height) < tolerance:
            name = "back"
        elif abs(y - 8 * height) < tolerance:
            name = "top"
        elif abs(x + 6 * height) < tolerance:
            name = "inlet"
        elif abs(x - 40 * height) < tolerance:
            name = "outlet"
        elif abs(y - height) < tolerance and x < 0:
            name = "upstream_wall"
        elif abs(x) < tolerance and y < height:
            name = "step_wall"
        elif abs(y) < tolerance and x > 0:
            name = "downstream_wall"
        else:
            raise ValueError(f"unclassified boundary face at {xyz}")
        patches[name].append(fid)
    if any(not f for f in patches.values()):
        raise ValueError("empty critical patch")
    fv = np.asarray([v for f in faces for v in f], dtype=np.uint64)
    fo = np.zeros(len(faces) + 1, dtype=np.uint64)
    for i, f in enumerate(faces):
        fo[i + 1] = fo[i] + len(f)
    cf = np.asarray([v for face in cell_faces for v in face], dtype=np.uint64)
    co = np.zeros(len(cell_faces) + 1, dtype=np.uint64)
    for i, f in enumerate(cell_faces):
        co[i + 1] = co[i] + len(f)
    own = np.asarray(owners, dtype=np.uint64)
    nei = np.asarray(neighbours, dtype=np.int64)
    topology = FNV_OFFSET
    for a in (fv, fo, own, nei, cf, co):
        topology = fnv(topology, a)
    geometry = FNV_OFFSET
    for comp in range(3):
        geometry = fnv(geometry, np.asarray(points[:, comp], dtype=np.float64))
    mh = topology
    for comp in range(3):
        mh = fnv(mh, np.asarray(points[:, comp], dtype=np.float64))
    patch_ids = []
    patch_offsets = [0]
    patch_metadata = []
    for name, ids in patches.items():
        patch_ids += ids
        patch_offsets.append(len(patch_ids))
        patch_metadata.append(f"{name}:{ids[0]}:{len(ids)}:{PATCH_TYPES[name]}")
    return {
        "points": points,
        "face_vertices": fv,
        "face_offsets": fo,
        "owner": own,
        "neighbour": nei,
        "cell_faces": cf,
        "cell_offsets": co,
        "patch_face_ids": np.asarray(patch_ids, dtype=np.uint64),
        "patch_face_offsets": np.asarray(patch_offsets, dtype=np.uint64),
        "patch_metadata": ";".join(patch_metadata),
        "topology_hash": f"{topology:016x}",
        "geometry_hash": f"{geometry:016x}",
        "mesh_hash": f"{mh:016x}",
        "n_cells": len(cell_faces),
        "n_faces": len(faces),
    }


def make_setup(height: float, reynolds: float, velocity: float) -> str:
    density = 1.2
    viscosity = density * velocity * height / reynolds
    setup = CaseSetup(
        source=SourceInfo(
            solver="cfdx",
            version="bfs-smoke-v1",
            case_name="TH03_debug_only",
            format="cfdx",
        )
    )
    setup.physics_model = "incompressible_laminar"  # experimental SST opt-in
    setup.turbulence_model = "k_omega_sst"
    setup.energy_model = "isothermal"
    setup.solver_mode = "steady"
    setup.units = "SI"
    setup.materials = [
        MaterialSpec(
            name="air",
            density=density,
            dynamic_viscosity=viscosity,
            turbulence_model="k_omega_sst",
        )
    ]
    k = 1.5 * (velocity * 0.05) ** 2
    omega = (k**0.5) / (0.09**0.25 * height)
    names = list(PATCH_TYPES)
    bc = []
    for name in names:
        if name == "inlet":
            spec = BoundarySpec(
                patch_name=name,
                type=BCType.INLET,
                value_type=BCValueType.FIXED,
                velocity_vector=[velocity, 0, 0],
                velocity_magnitude=velocity,
                turbulence_intensity=0.05,
                turbulence_length_scale=height,
                sst_scalar_conditions=[
                    SSTScalarBoundarySpec(
                        field="k", type=SSTScalarBoundaryType.FIXED_VALUE, value=k
                    ),
                    SSTScalarBoundarySpec(
                        field="omega",
                        type=SSTScalarBoundaryType.FIXED_VALUE,
                        value=omega,
                    ),
                ],
            )
        elif name == "outlet":
            spec = BoundarySpec(
                patch_name=name,
                type=BCType.PRESSURE_OUTLET,
                value_type=BCValueType.OUTLET_PRESSURE,
                pressure=0,
                sst_scalar_conditions=[
                    SSTScalarBoundarySpec(
                        field="k", type=SSTScalarBoundaryType.ZERO_GRADIENT
                    ),
                    SSTScalarBoundarySpec(
                        field="omega", type=SSTScalarBoundaryType.ZERO_GRADIENT
                    ),
                ],
            )
        elif name in ("front", "back"):
            spec = BoundarySpec(
                patch_name=name,
                type=BCType.EMPTY,
                value_type=BCValueType.ZERO_GRADIENT,
                sst_scalar_conditions=[
                    SSTScalarBoundarySpec(
                        field="k", type=SSTScalarBoundaryType.ZERO_GRADIENT
                    ),
                    SSTScalarBoundarySpec(
                        field="omega", type=SSTScalarBoundaryType.ZERO_GRADIENT
                    ),
                ],
            )
        else:
            spec = BoundarySpec(
                patch_name=name,
                type=BCType.WALL,
                value_type=BCValueType.WALL_NO_SLIP,
                velocity_vector=[0, 0, 0],
                sst_scalar_conditions=[
                    SSTScalarBoundarySpec(
                        field="k", type=SSTScalarBoundaryType.WALL_K_ZERO
                    ),
                    SSTScalarBoundarySpec(
                        field="omega", type=SSTScalarBoundaryType.WALL_OMEGA_SST
                    ),
                ],
            )
        bc.append(spec)
    setup.boundary_conditions = bc
    setup.initial_condition = InitialCondition(
        velocity=velocity,
        pressure=0,
        velocity_vector=[velocity, 0, 0],
        scalar_fields={"k": k, "omega": omega},
    )
    setup.numerics.max_iterations = 250
    setup.numerics.raw_settings["sst_transport_max_iterations"] = "600"
    setup.numerics.selection = NumericalSelectionConfig(
        entries=[
            NumericalSelection(
                family="gradient", configuration_key="numerics.gradient.gauss"
            ),
            NumericalSelection(
                family="convection", configuration_key="numerics.convection.upwind"
            ),
            NumericalSelection(
                family="pressure_velocity", configuration_key="pressure_velocity.simple"
            ),
            NumericalSelection(
                family="linear_solver", configuration_key="linear.fgmres"
            ),
        ],
        required_families=["gradient", "convection", "pressure_velocity"],
    )
    return json.dumps(setup.model_dump(mode="json"))


def write_case(
    path: Path,
    *,
    n_up: int = 4,
    n_down: int = 16,
    ny: int = 10,
    step_cells: int = 2,
    height: float = 0.0127,
    reynolds: float = 36000,
    velocity: float = 44.2,
    wall_stretch: float = 0.0,
) -> dict:
    if not (
        np.isfinite(reynolds)
        and reynolds > 0
        and np.isfinite(velocity)
        and velocity > 0
    ):
        raise ValueError("invalid Reynolds or velocity")
    mesh = make_mesh(n_up, n_down, ny, step_cells, height, wall_stretch)
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h5:
        for key, value in {
            "format": "CFDX-HDF5-mesh-v1",
            "format_version": "1",
            "schema_version": "1",
            "cfdx_version": "0.7",
            "mesh_topology": "cfdx-csr-v1",
            "n_points": str(len(mesh["points"])),
            "n_faces": str(mesh["n_faces"]),
            "n_cells": str(mesh["n_cells"]),
            "topology_hash": mesh["topology_hash"],
            "geometry_hash": mesh["geometry_hash"],
            "mesh_hash": mesh["mesh_hash"],
            "boundary_patches": mesh["patch_metadata"],
            "source_solver": "cfdx",
            "source_format": "cfdx",
            "source_version": "bfs-smoke-v1",
            "source_case_path": "",
            "source_case_name": "TH03_debug_only",
            "case_setup_json": make_setup(height, reynolds, velocity),
        }.items():
            h5.attrs[key] = value
        for key, value in mesh.items():
            if isinstance(value, np.ndarray):
                h5.create_dataset(key, data=value)
    return {
        "path": str(path),
        "cells": mesh["n_cells"],
        "faces": mesh["n_faces"],
        "reference_velocity_m_s": velocity,
        "reynolds_step": reynolds,
        "step_height_m": height,
        "wall_stretch": wall_stretch,
        "physical_qualification": "NOT_QUALIFIED",
        "reason": "uniform inlet, simplified channel geometry and no wall y+ or grid-convergence evidence",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--upstream", type=int, default=4)
    parser.add_argument("--downstream", type=int, default=16)
    parser.add_argument("--ny", type=int, default=10)
    parser.add_argument("--step-cells", type=int, default=2)
    parser.add_argument(
        "--wall-stretch",
        type=float,
        default=0.0,
        help="tanh clustering at y=0,h,8h, 0 means uniform",
    )
    args = parser.parse_args()
    print(
        json.dumps(
            write_case(
                args.output,
                n_up=args.upstream,
                n_down=args.downstream,
                ny=args.ny,
                step_cells=args.step_cells,
                wall_stretch=args.wall_stretch,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
