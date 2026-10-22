"""T2 heated-channel Nusselt calibration campaign (N9 fine-gap audit, T2).

Standalone Python driver (no C++ test-tree dependency): it builds the
heated-channel case on the full production `.cfdx.h5` path (mirroring
tests/python/test_n9_poiseuille_body_force_production.py) and runs the
production solver once per refinement mesh, then reports the developed
Nusselt number against the analytic two-plate constant-flux reference
Nu_ref = 140/17.

Locked case decisions (N9/T2 audit, heated-channel campaign #884):
  * channel x in [0, L=10], y in [0, 1], single z layer (EMPTY front/back)
  * uniform body force G = 1.0 per unit mass, nu = 0.1 (rho = 1, mu = 0.1)
  * both plates heat with the same uniform flux q'' = 1.0
  * Pr = 1 (k = 0.1, cp = 1.0), T_in = T_init = 300 K, energy "conjugate"
  * Nu_cell = q''*D_h/(k*(Tw - Tb)) with D_h = 2*H = 2 and the half-cell
    wall correction Tw = T_wall_cell + (q''/k)*(dy/2)
  * Nu_ref = 140/17 (exact, both plates at constant uniform flux)
  * developed region = downstream half, measured at cell centres x_c >= 7.5
  * guard: body-force Poiseuille analytic L2 <= 5e-3 on the probe column,
    checked at EVERY mesh before any Nusselt number is trusted
  * external energy check: mdot*cp*(Tb_out - Tb_in) vs q''*2*L, with mdot
    measured from the u_x probes (not an oracle)

The script never modifies the solver. Evidence (JSON + TXT) is always
written to --output-dir before exiting; the exit code is 1 only when a
mesh fails to run/converge (the err_rel measurement itself is not gated).
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from pathlib import Path

import h5py
import numpy as np

L = 10.0  # channel length (m)
H = 1.0  # channel height (m)
BODY_FORCE = (1.0, 0.0, 0.0)  # per unit mass (m/s^2), x-driven
DENSITY = 1.0  # kg/m^3
DYNAMIC_VISCOSITY = 0.1  # Pa s -> nu = 0.1 m^2/s
THERMAL_CONDUCTIVITY = 0.1  # W/(m K)
SPECIFIC_HEAT = 1.0  # J/(kg K)
HEAT_FLUX = 1.0  # W/m^2, both plates, positive toward the domain
T_IN = 300.0  # K (inlet Dirichlet temperature and initial temperature)
MESHES = [(80, 16), (160, 32), (320, 64)]  # (nx, ny) refinement series
MAX_ITERATIONS = 3000
L2_GATE = 5e-3  # Poiseuille analytic gate, identical to the N9 matrix
DH = 2.0 * H  # hydraulic diameter of the parallel-plate channel (m)
NU_REF = 140.0 / 17.0  # exact two-plate constant-flux Nusselt number
X_DEV = 7.5  # developed-region threshold on cell-centre x (m)
ERR_TARGET = 0.05  # relative Nu error target on the developed region

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


def _build_channel(nx: int, ny: int) -> dict:
    """Hexahedral single-layer channel [0, L] x [0, 1] x [0, 1]."""
    plane = (nx + 1) * (ny + 1)
    points = np.zeros((2 * plane, 3), dtype=np.float64)

    def pid(i: int, j: int, k: int) -> int:
        return (j * (nx + 1) + i) * 2 + k

    for j in range(ny + 1):
        y = j / ny
        for i in range(nx + 1):
            x = i / nx * L
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

    # Patch assignment is geometric, exactly like the N9 matrix mesh builder.
    patch_ids: dict[str, list[int]] = {name: [] for name in _PATCH_ORDER}
    eps = 1e-9
    for fid, face in enumerate(faces):
        if neighbour[fid] >= 0:
            continue
        x, y, z = points[face].mean(axis=0)
        if abs(y) < eps:
            name = "bottom"
        elif abs(y - H) < eps:
            name = "top"
        elif abs(z) < eps:
            name = "front"
        elif abs(z - 1.0) < eps:
            name = "back"
        elif abs(x) < eps:
            name = "inlet"
        elif abs(x - L) < eps:
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


def _ensure_cfdx_package() -> None:
    """Import the repo cfdx package, never scripts/cfdx.py (the CLI).

    Python prepends this script's directory (scripts/) to sys.path when
    running this file, and the repo keeps its cfdx.py CLI entry point
    in scripts/, which shadows the cfdx package: "from cfdx.io.schema
    import ..." then fails with "'cfdx' is not a package" (CI pytest
    never sees this because it runs from the repo root). Drop the
    script directory from sys.path and prepend the in-repo package root
    (src/cfdx/python) so the import resolves to the real package.
    """
    here = Path(__file__).resolve().parent
    sys.path[:] = [p for p in sys.path if p and Path(p).resolve() != here]
    pkg_root = here.parent / "src" / "cfdx" / "python"
    if (pkg_root / "cfdx" / "io" / "schema.py").is_file():
        sys.path.insert(0, str(pkg_root))


def _case_setup_json(nx: int, ny: int) -> str:
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

    case_name = f"t2_heated_channel_{nx}x{ny}"
    setup = CaseSetup(
        source=SourceInfo(
            solver="cfdx",
            version="t2-nu-calibration",
            case_name=case_name,
            format="cfdx",
        )
    )
    setup.physics_model = "incompressible_laminar"
    setup.turbulence_model = "laminar"
    setup.energy_model = "conjugate"
    setup.solver_mode = "steady"
    setup.units = "SI"
    setup.ref_temperature = T_IN

    setup.materials = [
        MaterialSpec(
            name="fluid",
            density=DENSITY,
            dynamic_viscosity=DYNAMIC_VISCOSITY,
            thermal_conductivity=THERMAL_CONDUCTIVITY,
            specific_heat=SPECIFIC_HEAT,
        )
    ]
    setup.body_force = list(BODY_FORCE)

    zero = [0.0, 0.0, 0.0]
    # The thermal closure under test: one Dirichlet inlet temperature plus
    # uniform wall heat fluxes on both plates (raw_params["heat_flux"], the
    # production temperature-BC path). Front/back are EMPTY and skipped.
    setup.boundary_conditions = [
        BoundarySpec(
            patch_name="inlet",
            type=BCType.INLET,
            value_type=BCValueType.ZERO_GRADIENT,
            temperature=T_IN,
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
            raw_params={"heat_flux": repr(HEAT_FLUX)},
        ),
        BoundarySpec(
            patch_name="top",
            type=BCType.WALL,
            value_type=BCValueType.FIXED,
            velocity_vector=list(zero),
            raw_params={"heat_flux": repr(HEAT_FLUX)},
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
        velocity=0.0,
        pressure=0.0,
        velocity_vector=list(zero),
        temperature=T_IN,
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


def _write_case(path: Path, nx: int, ny: int) -> None:
    mesh = _build_channel(nx, ny)
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
        h5.attrs["source_version"] = "t2-nu-calibration"
        h5.attrs["source_case_path"] = ""
        h5.attrs["source_case_name"] = f"t2_heated_channel_{nx}x{ny}"
        h5.attrs["case_setup_json"] = _case_setup_json(nx, ny)

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


def _parse_temperature_vtu(path: Path) -> np.ndarray:
    """Cell-centred temperature field, cell order c = j*nx + i."""
    tree = ET.parse(str(path))
    text = None
    for data_array in tree.getroot().iter("DataArray"):
        if data_array.get("Name") == "temperature":
            text = data_array.text
            break
    if text is None:
        raise RuntimeError("no temperature DataArray in the energy VTU")
    return np.asarray([float(tok) for tok in text.split()], dtype=np.float64)


def _clip_error(text: str, head: int = 600, tail: int = 3000) -> str:
    """Keep head AND tail of a long message; never drop the stop reason."""
    if len(text) <= head + tail:
        return text
    middle = len(text) - head - tail
    return f"{text[:head]}\n...[{middle} chars truncated]...\n{text[-tail:]}"


def _run_mesh(
    solver: Path, nx: int, ny: int, run_dir: Path, evidence_dir: Path
) -> dict:
    """One production solve at (nx, ny); raises on any execution failure."""
    run_dir.mkdir(parents=True, exist_ok=True)
    case_name = f"t2_heated_channel_{nx}x{ny}"
    case_path = run_dir / f"{case_name}.cfdx.h5"
    _write_case(case_path, nx, ny)

    # Mid-x column probes: the plane-channel flow is x-invariant, so the
    # probe column reconstructs the u_x(y) profile for the whole field.
    # One u_x probe per y row, located exactly at the cell centres.
    probes: list[str] = []
    x_probe = (nx // 2 + 0.5) / nx * L
    for j in range(ny):
        y = (j + 0.5) / ny
        probes.append(f"u_y{j:02d}:{x_probe!r},{y!r},0.5:u_x")

    probe_csv = run_dir / "probes.csv"
    completed = subprocess.run(
        [
            str(solver),
            "--mesh",
            str(case_path),
            "--output-dir",
            str(run_dir),
            "--iterations",
            str(MAX_ITERATIONS),
            *[item for probe in probes for item in ("--probe", probe)],
            "--probe-csv",
            str(probe_csv),
        ],
        capture_output=True,
        text=True,
        timeout=3600,
        check=False,
    )
    diagnostics = completed.stdout + completed.stderr
    # stdout carries the run summary (carrier verdict, energy-stage trace,
    # "Converged YES/NO Iterations N", "Energy stage skipped ..."). In run
    # 39778029174 the 200k combined tail was 100% stderr and hid the coarse
    # meshes' stdout summaries entirely. Archive both streams separately,
    # plus every machine-readable json artifact the solver wrote
    # (execution.json / convergence.json carry the exact per-stage status
    # and stop reason).
    (run_dir / f"solver_{nx}x{ny}.log").write_text(diagnostics, encoding="utf-8")
    (evidence_dir / f"t2_solver_tail_{nx}x{ny}.log").write_text(
        diagnostics[-200000:], encoding="utf-8"
    )
    (evidence_dir / f"t2_solver_stdout_{nx}x{ny}.log").write_text(
        completed.stdout[-400000:], encoding="utf-8"
    )
    for artifact in sorted(run_dir.glob("*.json")):
        (evidence_dir / f"t2_{artifact.stem}_{nx}x{ny}.json").write_text(
            artifact.read_text(encoding="utf-8"), encoding="utf-8"
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"solver rc={completed.returncode}; stdout tail:\n"
            f"{completed.stdout[-1500:]}\n...[truncated]...\n"
            f"stderr tail:\n{completed.stderr[-2500:]}"
        )
    for marker in ("Resolved numerical selections:", "Converged YES"):
        if marker not in diagnostics:
            raise RuntimeError(f"missing solver marker {marker!r}")
    if not probe_csv.is_file():
        raise RuntimeError("probe CSV missing")

    final: dict[str, float] = {}
    with probe_csv.open(newline="") as stream:
        for row in csv.reader(stream):
            if not row or row[0].startswith("#"):
                continue
            name, _iteration, _time, value = row
            final[name] = float(value)  # rows are iteration-ordered per probe
    if len(final) != ny:
        raise RuntimeError(f"expected {ny} probe rows, got {len(final)}")
    u = np.asarray([final[f"u_y{j:02d}"] for j in range(ny)], dtype=np.float64)
    if not np.all(np.isfinite(u)):
        raise RuntimeError("non-finite probe velocities")

    # Guard: analytic Poiseuille profile on the probe column, every mesh,
    # BEFORE any Nusselt number is trusted.
    nu_fluid = DYNAMIC_VISCOSITY / DENSITY
    g = BODY_FORCE[0]
    squared = 0.0
    for j in range(ny):
        y = (j + 0.5) / ny
        exact = g * y * (1.0 - y) / (2.0 * nu_fluid)
        error = u[j] - exact
        squared += error * error
    l2 = math.sqrt(squared / ny)
    if l2 > L2_GATE:
        raise RuntimeError(f"Poiseuille L2 gate failed: {l2:.6e} > {L2_GATE}")

    # Temperature field from the energy stage VTU (cell order c = j*nx + i).
    vtu = run_dir / "energy_result_0.vtu"
    if not vtu.is_file():
        raise RuntimeError("energy_result_0.vtu missing (energy stage did not run)")
    temps = _parse_temperature_vtu(vtu)
    if temps.size != nx * ny:
        raise RuntimeError(f"expected {nx * ny} cell temperatures, got {temps.size}")
    if not np.all(np.isfinite(temps)):
        raise RuntimeError("non-finite cell temperatures")
    temp = temps.reshape(ny, nx)

    # Bulk temperature per column from the probe velocity profile (the flow
    # is x-invariant; dy cancels in the discrete bulk ratio).
    u_sum = float(u.sum())
    if u_sum <= 0.0:
        raise RuntimeError("non-positive bulk velocity")
    tb = (u[:, None] * temp).sum(axis=0) / u_sum

    # Half-cell wall correction: Tw = T_wall_cell + (q''/k)*(dy/2).
    dy = H / ny
    wall_offset = (HEAT_FLUX / THERMAL_CONDUCTIVITY) * (dy / 2.0)
    tw_bottom = temp[0, :] + wall_offset
    tw_top = temp[-1, :] + wall_offset
    if np.any(tw_bottom - tb <= 0.0) or np.any(tw_top - tb <= 0.0):
        raise RuntimeError("non-positive wall-to-bulk temperature difference")

    nu_bottom = HEAT_FLUX * DH / (THERMAL_CONDUCTIVITY * (tw_bottom - tb))
    nu_top = HEAT_FLUX * DH / (THERMAL_CONDUCTIVITY * (tw_top - tb))
    nu_local = 0.5 * (nu_bottom + nu_top)

    # Developed region: downstream half at cell centres, x_c >= 7.5.
    x_c = (np.arange(nx) + 0.5) / nx * L
    dev = x_c >= X_DEV
    if not dev.any():
        raise RuntimeError("empty developed region")
    nu_dev = float(nu_local[dev].mean())
    nu_bottom_dev = float(nu_bottom[dev].mean())
    nu_top_dev = float(nu_top[dev].mean())
    err_rel = abs(nu_dev - NU_REF) / NU_REF

    # External energy check of the installation: enthalpy rise vs the total
    # injected wall power (both plates, per unit depth). mdot comes from the
    # probes, not an oracle. Reported as a diagnostic, not gated.
    mdot = DENSITY * float((u * dy).sum())
    q_wall = HEAT_FLUX * 2.0 * L
    enthalpy = mdot * SPECIFIC_HEAT * (float(tb[-1]) - float(tb[0]))
    balance_rel = abs(enthalpy - q_wall) / q_wall

    summary: dict = {}
    summary_path = run_dir / "energy_summary.json"
    if summary_path.is_file():
        try:
            summary = json.loads(summary_path.read_text())
        except (OSError, ValueError):
            summary = {}

    return {
        "l2_poiseuille": l2,
        "probe_x": x_probe,
        "mdot": mdot,
        "tb_in": float(tb[0]),
        "tb_out": float(tb[-1]),
        "q_wall_per_unit_depth": q_wall,
        "enthalpy_rise_per_unit_depth": enthalpy,
        "balance_rel": balance_rel,
        "nu_dev": nu_dev,
        "nu_bottom_dev": nu_bottom_dev,
        "nu_top_dev": nu_top_dev,
        "nu_local_min": float(nu_local[dev].min()),
        "nu_local_max": float(nu_local[dev].max()),
        "n_dev_columns": int(dev.sum()),
        "dx": L / nx,
        "dy": dy,
        "err_rel": err_rel,
        "energy_summary": summary,
    }


def _fit_order(dxs: list[float], errs: list[float]) -> tuple[float, float]:
    """Least-squares fit of err ~ C * dx^p on the (dx, err) series."""
    xs = [math.log(d) for d in dxs]
    ys = [math.log(max(e, 1e-16)) for e in errs]
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    den = sum((x - mx) ** 2 for x in xs)
    if den == 0.0:
        return float("nan"), float("nan")
    p = num / den
    c = math.exp(my - p * mx)
    return p, c


def _write_evidence(out_dir: Path, campaign: dict) -> None:
    run_id = campaign["run_id"]
    json_path = out_dir / f"t2_nu_calibration_{run_id}.json"
    txt_path = out_dir / f"t2_nu_calibration_{run_id}.txt"
    json_path.write_text(json.dumps(campaign, indent=2) + "\n", encoding="utf-8")

    lines: list[str] = []
    lines.append("T2 heated-channel Nu calibration campaign")
    lines.append(f"run_id: {run_id}")
    lines.append(f"solver: {campaign['solver']}")
    lines.append(
        "case: L=10 H=1 G=1 nu=0.1 rho=1 k=0.1 cp=1 q''=1 T_in=300 "
        "energy=conjugate iterations=3000 meshes=80x16,160x32,320x64"
    )
    lines.append(f"Nu_ref = 140/17 = {NU_REF:.10f} ; err_rel target <= {ERR_TARGET}")
    lines.append(f"Poiseuille L2 gate: {L2_GATE} ; developed region: x_c >= {X_DEV}")
    lines.append("")
    for entry in campaign["meshes"]:
        if entry["status"] == "PASS":
            lines.append(
                f"mesh {entry['nx']}x{entry['ny']}: PASS "
                f"l2={entry['l2_poiseuille']:.3e} "
                f"Nu_dev={entry['nu_dev']:.6f} err_rel={entry['err_rel']:.6f} "
                f"(bottom {entry['nu_bottom_dev']:.6f}, top {entry['nu_top_dev']:.6f})"
            )
            lines.append(
                f"  mdot={entry['mdot']:.6f} tb_in={entry['tb_in']:.3f} "
                f"tb_out={entry['tb_out']:.3f} "
                f"balance_rel={entry['balance_rel']:.3e} "
                f"n_dev_columns={entry['n_dev_columns']}"
            )
        else:
            lines.append(f"mesh {entry['nx']}x{entry['ny']}: FAIL")
            lines.append(f"  error: {entry['error']}")
    fit = campaign.get("order_fit")
    if fit is not None:
        lines.append("")
        lines.append(
            f"apparent order fit err ~ C*dx^p: p={fit['p']:.3f} C={fit['C']:.6f}"
        )
    lines.append("")
    lines.append(f"verdict: {campaign['verdict']}")
    lines.append("interpretation (locked protocol):")
    lines.append(
        "- err_rel <= 0.05 on the finest mesh -> wall-treatment defect dominates"
        " the fine-gap budget; no convective-scheme correction required for"
        " the +9.93% audit (T2 closes with this evidence)."
    )
    lines.append(
        "- err_rel > 0.05 and/or apparent order p ~ 1 -> upwind artificial"
        " diffusion dominates; T3 = repeat the campaign with a higher-order"
        " convective scheme."
    )
    txt_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    _ensure_cfdx_package()
    parser = argparse.ArgumentParser(
        description="T2 heated-channel Nu calibration campaign"
    )
    parser.add_argument("--solver", required=True, help="production solver binary")
    parser.add_argument("--output-dir", required=True, help="evidence directory")
    parser.add_argument("--run-id", default="local", help="campaign run id")
    args = parser.parse_args()

    solver = Path(args.solver)
    if not solver.is_file():
        print(f"T2NU_CAMPAIGN FAIL solver_missing={solver}")
        return 1
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    campaign = {
        "run_id": args.run_id,
        "solver": str(solver),
        "case": {
            "L": L,
            "H": H,
            "body_force": list(BODY_FORCE),
            "density": DENSITY,
            "dynamic_viscosity": DYNAMIC_VISCOSITY,
            "thermal_conductivity": THERMAL_CONDUCTIVITY,
            "specific_heat": SPECIFIC_HEAT,
            "heat_flux": HEAT_FLUX,
            "t_in": T_IN,
            "energy_model": "conjugate",
            "max_iterations": MAX_ITERATIONS,
            "meshes": [list(m) for m in MESHES],
            "l2_gate": L2_GATE,
            "dh": DH,
            "x_dev": X_DEV,
        },
        "nu_ref": NU_REF,
        "err_target": ERR_TARGET,
        "meshes": [],
    }

    ok = True
    with tempfile.TemporaryDirectory(prefix="t2_nu_") as tmp:
        work = Path(tmp)
        for nx, ny in MESHES:
            entry: dict = {"nx": nx, "ny": ny}
            try:
                result = _run_mesh(solver, nx, ny, work / f"mesh_{nx}x{ny}", out_dir)
                entry.update(result)
                entry["status"] = "PASS"
                print(
                    f"T2NU_MESH nx={nx} ny={ny} Nu_dev={result['nu_dev']:.6f} "
                    f"err_rel={result['err_rel']:.6f}"
                )
                print(
                    f"T2NU_BALANCE nx={nx} ny={ny} mdot={result['mdot']:.6f} "
                    f"balance_rel={result['balance_rel']:.6e} "
                    f"l2={result['l2_poiseuille']:.3e}"
                )
            except Exception as exc:  # noqa: BLE001 - evidence-first reporting
                message = _clip_error(str(exc))
                entry["status"] = "FAIL"
                entry["error"] = message
                ok = False
                print(f"T2NU_MESH nx={nx} ny={ny} FAILED error={message}")
            campaign["meshes"].append(entry)

    if ok:
        dxs = [entry["dx"] for entry in campaign["meshes"]]
        errs = [entry["err_rel"] for entry in campaign["meshes"]]
        p, c = _fit_order(dxs, errs)
        campaign["order_fit"] = {"dx": dxs, "err": errs, "p": p, "C": c}
        finest = campaign["meshes"][-1]
        within = finest["err_rel"] <= ERR_TARGET
        campaign["verdict"] = (
            "err_rel_within_target" if within else "err_rel_above_target"
        )
        print(f"T2NU_ORDER p={p:.3f} C={c:.6f}")
        print(
            f"T2NU_CAMPAIGN PASS finest_err_rel={finest['err_rel']:.6f} "
            f"target={ERR_TARGET} verdict={campaign['verdict']}"
        )
    else:
        campaign["verdict"] = "execution_failure"
        print("T2NU_CAMPAIGN FAIL execution_failure")

    _write_evidence(out_dir, campaign)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
