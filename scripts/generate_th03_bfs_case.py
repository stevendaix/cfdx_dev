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
import shutil
import hashlib
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


SCHUR_MODELS = ("pcd", "block_local", "lsc", "bfbt", "simple", "simplec")
COUPLED_PRECONDITIONERS = ("auto", "pcd", "coupled_block_schur", "mgr", "braess_sarazin")


def validate_solver_selection(algorithm: str, schur: str | None,
                              preconditioner: str | None) -> None:
    """Fail before the HDF5 write; never silently substitute a Schur or PC."""
    if algorithm not in ("simple", "simplec", "coupled"):
        raise ValueError("algorithm must be simple, simplec or coupled")
    if algorithm != "coupled":
        if schur is not None or preconditioner is not None:
            raise ValueError("segregated BFS must not request a Schur or coupled preconditioner")
        return
    if schur not in SCHUR_MODELS:
        raise ValueError("coupled BFS requires one of six explicit Schur models")
    if preconditioner not in COUPLED_PRECONDITIONERS:
        raise ValueError("coupled BFS requires an explicit supported preconditioner policy")
    if preconditioner == "braess_sarazin" and schur != "block_local":
        raise ValueError("Braess-Sarazin requires block_local Schur")


def make_setup(height: float, reynolds: float, velocity: float,
               *, algorithm: str = "simple", schur: str | None = None,
               preconditioner: str | None = None, iterations: int = 250) -> str:
    validate_solver_selection(algorithm, schur, preconditioner)
    if iterations < 1 or iterations > 100000:
        raise ValueError("BFS iteration budget must be in [1,100000]")
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
    setup.numerics.max_iterations = iterations
    setup.numerics.raw_settings["sst_transport_max_iterations"] = "600"
    setup.numerics.raw_settings["v8_global_nonlinear_residual"] = "true"
    setup.numerics.raw_settings["v7_post_sst_momentum_gate"] = "true"
    selections = [
        NumericalSelection(
            family="gradient", configuration_key="numerics.gradient.gauss"
        ),
        NumericalSelection(
            family="convection", configuration_key="numerics.convection.second_order_upwind"
        ),
        NumericalSelection(
            family="pressure_velocity", configuration_key=f"pressure_velocity.{algorithm}"
        ),
        NumericalSelection(
            family="linear_solver", configuration_key="linear.fgmres"
        ),
    ]
    if algorithm == "coupled":
        selections.append(NumericalSelection(
            family="schur", configuration_key=f"schur.{schur}"
        ))
        if preconditioner != "auto":
            selections.append(NumericalSelection(
                family="preconditioner",
                configuration_key=f"preconditioner.{preconditioner}"
            ))
    setup.numerics.selection = NumericalSelectionConfig(
        entries=selections,
        required_families=["gradient", "convection", "pressure_velocity"]
            + (["schur"] if algorithm == "coupled" else []),
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
    algorithm: str = "simple",
    schur: str | None = None,
    preconditioner: str | None = None,
    iterations: int = 250,
) -> dict:
    validate_solver_selection(algorithm, schur, preconditioner)
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
            "case_setup_json": make_setup(
                height, reynolds, velocity, algorithm=algorithm, schur=schur,
                preconditioner=preconditioner, iterations=iterations
            ),
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
        "algorithm": algorithm,
        "schur": schur,
        "preconditioner": preconditioner,
        "iterations": iterations,
        "geometry_hash": mesh["geometry_hash"],
        "topology_hash": mesh["topology_hash"],
        "physical_qualification": "NOT_QUALIFIED",
        "reason": "uniform inlet, simplified channel geometry and no wall y+ or grid-convergence evidence",
    }


def reconfigure_bfs_reference_case(
    source_path: Path, output_path: Path, *, algorithm: str,
    schur: str | None, preconditioner: str | None, iterations: int,
) -> dict:
    """Change only numerical selection and certificate opt-ins on an existing BFS case.

    The imported mesh, BCs, SST inlet data and initial condition stay byte-identical
    at the dataset level. A reference case's experimental fidelity is NOT certified
    merely because its mesh/setup file can be used for this campaign.
    """
    validate_solver_selection(algorithm, schur, preconditioner)
    if not (1 <= iterations <= 100000):
        raise ValueError("invalid BFS iteration budget")
    source_path, output_path = Path(source_path), Path(output_path)
    if not source_path.is_file() or source_path.resolve() == output_path.resolve():
        raise ValueError("reference input must exist and differ from configured case")
    with h5py.File(source_path, "r") as source:
        if source.attrs.get("format") != "CFDX-HDF5-mesh-v1":
            raise ValueError("unsupported BFS HDF5 source format")
        if any(not source.attrs.get(key) for key in
               ("topology_hash", "geometry_hash", "case_setup_json")):
            raise ValueError("BFS reference input requires mesh fingerprint and case setup")
        setup = json.loads(source.attrs["case_setup_json"])
        numerics = setup.get("numerics")
        if not isinstance(numerics, dict):
            raise ValueError("reference input has no numerical configuration")
        selection = numerics.get("selection")
        if not isinstance(selection, dict) or not isinstance(selection.get("entries"), list):
            raise ValueError("reference input has no explicit numerical selection")
        entries = selection["entries"]
        retained = [x for x in entries if x.get("family") not in
                    ("pressure_velocity", "schur", "preconditioner")]
        if not any(x.get("family") == "convection" and
                   x.get("configuration_key") in (
                       "numerics.convection.second_order_upwind",
                       "numerics.convection.tvd.vanleer", "numerics.convection.tvd.vanalbada",
                       "numerics.convection.tvd.mc", "numerics.convection.tvd.superbee")
                   for x in retained):
            raise ValueError("reference input must explicitly use a bounded order-two scheme")
        retained.append({"family":"pressure_velocity",
                         "configuration_key":f"pressure_velocity.{algorithm}"})
        if algorithm == "coupled":
            retained.append({"family":"schur", "configuration_key":f"schur.{schur}"})
            if preconditioner != "auto":
                retained.append({"family":"preconditioner",
                                 "configuration_key":f"preconditioner.{preconditioner}"})
        selection["entries"] = retained
        required = [x for x in selection.get("required_families", [])
                    if x not in ("schur", "preconditioner")]
        if algorithm == "coupled":
            required.append("schur")
        selection["required_families"] = list(dict.fromkeys(required))
        numerics["max_iterations"] = iterations
        raw_settings = numerics.setdefault("raw_settings", {})
        raw_settings["v8_global_nonlinear_residual"] = "true"
        raw_settings["v7_post_sst_momentum_gate"] = "true"
        topology, geometry = str(source.attrs["topology_hash"]), str(source.attrs["geometry_hash"])
        n_cells, n_faces = int(source.attrs["n_cells"]), int(source.attrs["n_faces"])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source_path, output_path)
    with h5py.File(output_path,"r+") as configured:
        configured.attrs["case_setup_json"] = json.dumps(setup, allow_nan=False)
        assert str(configured.attrs["topology_hash"]) == topology
        assert str(configured.attrs["geometry_hash"]) == geometry
    return {"path": str(output_path), "cells": n_cells, "faces": n_faces,
            "topology_hash": topology, "geometry_hash": geometry,
            "algorithm": algorithm, "schur": schur,
            "preconditioner": preconditioner,
            "input_kind": "EXTERNAL_BFS_REFERENCE_CANDIDATE",
            "source_case_sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "physical_qualification": "NOT_QUALIFIED",
            "reason": "experimental geometry, inlet profiles, wall y+ and grid study not yet independently qualified"}

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
    parser.add_argument("--algorithm", choices=("simple", "simplec", "coupled"), default="simple")
    parser.add_argument("--schur", choices=SCHUR_MODELS)
    parser.add_argument("--preconditioner", choices=COUPLED_PRECONDITIONERS)
    parser.add_argument("--iterations", type=int, default=250)
    args = parser.parse_args()
    if not str(args.output).endswith(".cfdx.h5"):
        parser.error("native production case must end in .cfdx.h5 to load physical boundary conditions")
    print(
        json.dumps(
            write_case(
                args.output,
                n_up=args.upstream,
                n_down=args.downstream,
                ny=args.ny,
                step_cells=args.step_cells,
                wall_stretch=args.wall_stretch,
                algorithm=args.algorithm, schur=args.schur,
                preconditioner=args.preconditioner, iterations=args.iterations,
            ),
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
