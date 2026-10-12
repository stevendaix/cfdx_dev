"""SU2 (.su2 mesh + .cfg config) adapter for CFDX I/O.

Parses the actual SU2 file format:
  - NDIME= <n>           dimensionality
  - NELEM= <n>           element count, then `<type> <n1> <n2> ... <cell_id>` lines
  - NPOIN= <n>           node count, then `x y z idx` lines
  - NMARK= <n>           boundary marker count
  - MARKER_TAG <name>   marker name
  - MARKER_ELEMS <n>    face count for this marker
    then `<type> <n1> ...` lines

.cfg format: `%-commented, `KEY = value` pairs.

Solution files:
  - solution.csv (SU2 CSV format with primitive variables per cell/vertex)
  - history.csv (convergence history)

Corresponds to the C++ header:
  src/cfdx/io/su2/su2_adapter.h
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional, Union

import numpy as np

from cfdx.io.interfaces import (
    SolverAdapter,
    ConversionResult,
)
from cfdx.io.schema import (
    BCType,
    BCValueType,
    CaseSetup,
    MaterialSpec,
    BoundarySpec,
    SourceInfo,
)
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.gap_analysis import GapAnalysis


# SU2 element type codes → (name, n_nodes)
_SU2_ELEM_TYPES: dict[int, tuple[str, int]] = {
    1: ("NODE", 1),
    3: ("LINE", 2),
    5: ("TRIANGLE", 3),
    6: ("QUADRILATERAL", 4),
    9: ("TETRAHEDRON", 4),
    10: ("HEXAHEDRON", 8),
    11: ("PRISM", 6),
    12: ("PYRAMID", 5),
}


# SU2 spells the second-order reconstruction across several keys, and renamed
# MUSCL to SLOPE_LIMITER in v8.  MUSCL_FLOW / SLOPE_LIMITER_FLOW_FLOW switch it
# on for the flow equations; MUSCL_LIMITER_FLOW, MUSCL_LIMITER, SLOPE_LIMITER and
# SLOPE_LIMITER_FLOW name the limiter.  The CFDX convection registry keys
# describe exactly that reconstruction — "MUSCL with psi=<limiter>" — so the
# limiter drives the mapping.
_SU2_MUSCL_LIMITERS: dict[str, str] = {
    "VAN_ALBADA": "tvd_vanalbada",
    "MINMOD": "tvd_minmod",
    "VAN_LEER": "tvd_vanleer",
    "SUPERBEE": "tvd_superbee",
    "MONOTONIZED_CENTRAL": "tvd_mc",
}

# SU2 limiter keys, most specific first.
_SU2_LIMITER_KEYS = (
    "MUSCL_LIMITER_FLOW",
    "SLOPE_LIMITER_FLOW",
    "MUSCL_LIMITER",
    "SLOPE_LIMITER",
)

# SU2 keys that enable the second-order reconstruction for the flow equations.
_SU2_MUSCL_ON_KEYS = ("MUSCL_FLOW", "MUSCL", "SLOPE_LIMITER_FLOW")

# SU2 default limiter when the reconstruction is enabled and no limiter key is
# present.
_SU2_DEFAULT_MUSCL_LIMITER = "VAN_ALBADA"

# SU2 convection families that carry no reconstruction at all.
_SU2_CENTRAL_SCHEMES = {"NO_UPWIND", "NO_CONV"}

# SU2 flux families that are plain first-order upwind when no second-order
# reconstruction is requested, so they map onto the CFDX upwind key.  Families
# with inherent sensor-based limiting (AUSM family, SLAU/SLAU2) and centre-based
# reconstructions (JST, FDS) are deliberately absent: none of them is
# first-order upwind, and CFDX has no equivalent key.
_SU2_FIRST_ORDER_UPWIND_SCHEMES = {"ROE", "HLLC", "HLLEM"}

_TRUE_VALUES = {"yes", "true", "1"}


def _first_param(params: dict[str, str], keys: tuple[str, ...]) -> str:
    for key in keys:
        value = params.get(key, "")
        if value.strip():
            return value
    return ""


def _canonical_convection_scheme(params: dict[str, str]) -> str:
    """Resolve SU2 convection settings onto a canonical CFDX convection family.

    SU2 splits the scheme over CONV_NUM_METHOD_FLOW (flux/scheme family), the
    MUSCL/SLOPE_LIMITER enable flag and the limiter name.  The CFDX registry
    keys describe the reconstruction and its limiter rather than the numerical
    flux, so an active MUSCL reconstruction maps through its limiter and a plain
    Roe-like flux maps to first-order upwind.

    Everything else is returned unresolved: centre-based and sensor-based
    families (JST, FDS, AUSM family, SLAU) and unbounded MUSCL.  CFDX has no
    ``numerics.convection`` key for any of them, and every key it does have
    carries an explicit boundedness claim, so mapping them would misstate the
    discretisation.  The caller turns an unresolved family into a blocking
    conversion gap.
    """
    scheme = params.get("CONV_NUM_METHOD_FLOW", "").strip().upper()
    if scheme in _SU2_CENTRAL_SCHEMES:
        return "central"

    if _first_param(params, _SU2_MUSCL_ON_KEYS).strip().lower() in _TRUE_VALUES:
        limiter = _first_param(params, _SU2_LIMITER_KEYS).strip().upper()
        if not limiter:
            limiter = _SU2_DEFAULT_MUSCL_LIMITER
        return _SU2_MUSCL_LIMITERS.get(limiter, "")

    if scheme:
        if scheme in _SU2_FIRST_ORDER_UPWIND_SCHEMES:
            return "upwind"
        return ""
    return ""


class Su2Adapter(SolverAdapter):
    solver_name = "SU2"
    format_name = "su2"

    def __init__(self, rules: Optional[MappingRules] = None) -> None:
        self.rules = rules or MappingRules.load_default()
        self.setup = CaseSetup(source=SourceInfo(solver="SU2", format="su2"))
        self._points: list[np.ndarray] = []
        self._elements: list[tuple[str, list[int]]] = []
        self._boundaries: list[dict] = []
        self._cfg_params: dict[str, str] = {}
        self._mesh: Optional[dict] = None

    def detect_source(self, case_path: str, info: SourceInfo) -> bool:
        p = Path(case_path)
        candidates = []
        if p.is_file():
            candidates.append(p)
        elif p.is_dir():
            for ext in ("*.su2", "*.cfg"):
                candidates.extend(sorted(p.glob(ext)))
        else:
            for f in (p, Path(str(p) + ".su2"), Path(str(p) + ".cfg")):
                if f.exists():
                    candidates.append(f)

        for cand in candidates:
            try:
                first = cand.read_text(errors="replace")[:2048]
            except OSError:
                continue
            if "NDIME=" in first or first.strip().startswith("%"):
                info.solver = "SU2"
                info.format = "su2" if cand.suffix == ".su2" else "cfg"
                info.case_path = str(cand.parent)
                info.case_name = cand.stem
                info.version = "unknown"
                return True
        return False

    def convert(self, case_path: Union[str, Path], result: ConversionResult) -> bool:
        base = Path(str(case_path))
        su2_file = base
        cfg_file = base

        if base.is_file():
            if base.suffix == ".su2":
                su2_file = base
                cfg_file = base.with_suffix(".cfg")
            elif base.suffix == ".cfg":
                cfg_file = base
                su2_file = base.with_suffix(".su2")
        else:
            su2_file = base / "mesh.su2"
            cfg_file = base / "config.cfg"

        mesh_ok = self.parse_mesh(str(su2_file)) if su2_file.exists() else False
        cfg_ok = self.parse_config(str(cfg_file)) if cfg_file.exists() else False

        if not mesh_ok:
            result.gap_report.unsupported_blocking(
                "mesh", "su2_mesh",
                f"Cannot parse SU2 mesh file: {su2_file}",
                "Ensure the .su2 file uses the standard NDIME/NELEM/NPOIN/NMARK format",
            )
            return False

        if not cfg_ok:
            result.gap_report.unsupported_nonblocking(
                "config", "su2_cfg",
                f"Cannot parse SU2 config file: {cfg_file}",
                "Config parsing is optional; mesh will still be imported",
            )

        # Try to import solution.csv if it exists (non-blocking)
        solution_file = su2_file.with_name(su2_file.stem + "_solution.csv")
        if not solution_file.exists():
            # Try common alternative names
            for alt_name in ["solution.csv", "flow.csv", "surface_flow.csv"]:
                alt = su2_file.with_name(alt_name)
                if alt.exists():
                    solution_file = alt
                    break

        if solution_file.exists():
            self.import_results(str(solution_file), result)
        else:
            result.gap_report.unsupported_nonblocking(
                "results", "solution_csv",
                f"No SU2 solution.csv found alongside mesh (looked for {solution_file.name})",
                "Run SU2 simulation to generate solution.csv, or provide it manually",
            )

        result.source = self.setup.source
        result.setup = self.setup.model_copy(deep=True)
        result.mesh = self._mesh
        self._populate_gap_report(result.gap_report)
        return not result.gap_report.has_blocking()

    def parse_mesh(self, su2_file: str) -> bool:
        text = Path(su2_file).read_text(errors="replace")
        lines = text.splitlines()
        idx = 0
        n_dim = 3
        points: list[list[float]] = []
        elements: list[tuple[str, list[int]]] = []
        boundaries: list[dict] = []

        while idx < len(lines):
            line = lines[idx].strip()
            upper = line.upper()

            if upper.startswith("NDIME="):
                n_dim = int(line.split("=")[1].strip())
                self.setup.mesh_info.dimension = n_dim
                idx += 1
            elif upper.startswith("NELEM="):
                n_elem = int(line.split("=")[1].strip())
                idx += 1
                for _ in range(n_elem):
                    if idx >= len(lines):
                        break
                    tokens = lines[idx].split()
                    idx += 1
                    if not tokens:
                        continue
                    try:
                        elem_type = int(tokens[0])
                    except ValueError:
                        continue
                    name = _SU2_ELEM_TYPES.get(elem_type, ("UNKNOWN", 0))[0]
                    n_nodes = _SU2_ELEM_TYPES.get(elem_type, ("UNKNOWN", 0))[1]
                    node_ids = [int(t) for t in tokens[1 : 1 + n_nodes]]
                    elements.append((name, node_ids))
                self._elements = elements
            elif upper.startswith("NPOIN="):
                n_points = int(line.split("=")[1].strip())
                idx += 1
                for _ in range(n_points):
                    if idx >= len(lines):
                        break
                    tokens = lines[idx].split()
                    idx += 1
                    if len(tokens) < n_dim + 1:
                        continue
                    coords = [float(t) for t in tokens[:n_dim]]
                    if n_dim == 2:
                        coords.append(0.0)
                    points.append(coords)
                self._points = [np.array(p) for p in points]
                self.setup.mesh_info.n_vertices = len(points)
            elif upper.startswith("NMARK="):
                n_mark = int(line.split("=")[1].strip())
                idx += 1
                for _ in range(n_mark):
                    if idx >= len(lines):
                        break
                    tag_line = lines[idx].strip()
                    idx += 1
                    if not tag_line.upper().startswith("MARKER_TAG"):
                        continue
                    marker_name = tag_line.split("=", 1)[1].strip()

                    if idx < len(lines) and lines[idx].strip().upper().startswith("MARKER_ELEMS="):
                        n_elems = int(lines[idx].split("=")[1].strip())
                        idx += 1
                    else:
                        n_elems = 0

                    face_nodes: list[list[int]] = []
                    for _ in range(n_elems):
                        if idx >= len(lines):
                            break
                        ftokens = lines[idx].split()
                        idx += 1
                        if not ftokens:
                            continue
                        try:
                            ftype = int(ftokens[0])
                        except ValueError:
                            continue
                        fn_nodes = _SU2_ELEM_TYPES.get(ftype, ("UNKNOWN", 0))[1]
                        face_nodes.append([int(t) for t in ftokens[1 : 1 + fn_nodes]])

                    boundaries.append({
                        "name": marker_name,
                        "n_elements": n_elems,
                        "node_ids": [n for face in face_nodes for n in face],
                        "face_nodes": face_nodes,
                    })
                self._boundaries = boundaries

                # Build boundary BC specs
                for bnd in boundaries:
                    bc = BoundarySpec(
                        patch_name=bnd["name"],
                        source_zone_name=bnd["name"],
                    )
                    mapped = self.rules.map_boundary_type(bnd["name"].lower())
                    bc.type = mapped
                    bc.value_type = self.rules.map_bc_value_type(bnd["name"].lower())
                    if bnd["name"].lower().find("wall") >= 0:
                        bc.type = BCType.WALL
                        bc.value_type = BCValueType.WALL_NO_SLIP
                    elif bnd["name"].lower().find("far") >= 0 or bnd["name"].lower().find("inlet") >= 0:
                        bc.type = BCType.INLET if "inlet" in bnd["name"].lower() else BCType.OUTLET
                        bc.value_type = BCValueType.FIXED
                    elif bnd["name"].lower().find("sym") >= 0:
                        bc.type = BCType.SYMMETRY
                        bc.value_type = BCValueType.ZERO_GRADIENT
                    self.setup.boundary_conditions.append(bc)
            else:
                idx += 1

        # Build mesh dict (meshio-style)
        cell_data = {}
        cells: dict[str, np.ndarray] = {}
        for elem_type in set(e[0] for e in self._elements):
            conn = [e[1] for e in self._elements if e[0] == elem_type]
            arr = np.array(conn, dtype=np.int32)
            meshio_type = _to_meshio_type(elem_type)
            if meshio_type:
                cells[meshio_type] = arr

        points_arr = np.array(self._points, dtype=np.float64) if self._points else np.empty((0, 3))

        mesh = {"points": points_arr, "cells": cells, "cell_data": cell_data}

        if self._points:
            self._mesh = mesh
            self.setup.mesh_info.n_vertices = len(self._points)
            types_found = sorted(set(e[0] for e in self._elements))
            self.setup.mesh_info.cell_types = types_found
            self.setup.mesh_info.n_cells = len(self._elements)
            self.setup.mesh_info.n_patches = len(self._boundaries)

            # Count boundary faces
            self.setup.mesh_info.n_boundary_faces = sum(b["n_elements"] for b in self._boundaries)

        self.setup.source.case_path = su2_file
        self.setup.source.case_name = Path(su2_file).stem
        return len(self._points) > 0

    def parse_config(self, cfg_file: str) -> bool:
        text = Path(cfg_file).read_text(errors="replace")
        for line in text.splitlines():
            stripped = line.strip()
            if not stripped or stripped.startswith("%"):
                continue
            eq = stripped.find("=")
            if eq == -1:
                continue
            key = stripped[:eq].strip()
            val = stripped[eq + 1:].strip()
            # Remove inline comments
            if ";" in val:
                val = val[: val.find(";")]
            val = val.strip()
            self._cfg_params[key] = val

        self.setup.source.case_path = cfg_file
        self.setup.source.case_name = Path(cfg_file).stem

        # Map known keys
        p = self._cfg_params
        if "MACH_NUMBER" in p:
            try:
                mach = float(p["MACH_NUMBER"])
                self.setup.initial_condition.velocity = mach * 340.0
            except ValueError:
                pass
        if "AOA" in p:
            try:
                alpha = float(p["AOA"]) * np.pi / 180.0
                self.setup.initial_condition.velocity_vector = [
                    float(np.cos(alpha)), float(np.sin(alpha)), 0.0
                ]
            except ValueError:
                pass
        if "FREESTREAM_PRESSURE" in p:
            try:
                self.setup.initial_condition.pressure = float(p["FREESTREAM_PRESSURE"])
            except ValueError:
                pass
        if "FREESTREAM_TEMPERATURE" in p:
            try:
                self.setup.initial_condition.temperature = float(p["FREESTREAM_TEMPERATURE"])
            except ValueError:
                pass
        if "SOLVER" in p:
            solver_val = p["SOLVER"].strip().lower()

            # SU2 uses distinct INC_* solver names for incompressible flows.
            # Check these before the generic "euler"/"navier" substrings:
            # otherwise INC_EULER and INC_NAVIER_STOKES would be classified
            # incorrectly as compressible.
            if solver_val == "inc_euler":
                self.setup.physics_model = "incompressible_euler"
                self.setup.energy_model = "isothermal"
            elif solver_val == "inc_navier_stokes":
                self.setup.physics_model = "incompressible_navier_stokes"
                self.setup.energy_model = "isothermal"
            elif solver_val == "inc_rans":
                self.setup.physics_model = "incompressible_rans"
                self.setup.energy_model = "isothermal"
            elif "euler" in solver_val:
                self.setup.physics_model = "compressible_euler"
                self.setup.energy_model = "compressible"
            elif "navier" in solver_val:
                self.setup.physics_model = "compressible_navier_stokes"
                self.setup.energy_model = "compressible"
            elif "heat" in solver_val:
                self.setup.physics_model = "heat_equation"
                self.setup.energy_model = "buoyant"

        # Incompressible SU2 initial conditions.
        if "INC_DENSITY_INIT" in p:
            try:
                rho = float(p["INC_DENSITY_INIT"])
                self.setup.ref_density = rho
                if self.setup.materials:
                    self.setup.materials[0].density = rho
                else:
                    self.setup.materials.append(MaterialSpec(name="fluid", density=rho))
            except ValueError:
                pass
        if "INC_VELOCITY_INIT" in p:
            try:
                raw = p["INC_VELOCITY_INIT"].strip().strip("()")
                velocity = [float(v.strip()) for v in raw.split(",")]
                self.setup.initial_condition.velocity_vector = velocity
                self.setup.initial_condition.velocity = float(np.linalg.norm(velocity))
                self.setup.ref_velocity = self.setup.initial_condition.velocity
            except (ValueError, TypeError):
                pass
        if "INC_DENSITY_REF" in p:
            try:
                self.setup.ref_density = float(p["INC_DENSITY_REF"])
            except ValueError:
                pass
        if "INC_VELOCITY_REF" in p:
            try:
                self.setup.ref_velocity = float(p["INC_VELOCITY_REF"])
            except ValueError:
                pass
        if "MU_CONSTANT" in p:
            try:
                mu = float(p["MU_CONSTANT"])
                if self.setup.materials:
                    self.setup.materials[0].dynamic_viscosity = mu
                else:
                    self.setup.materials.append(MaterialSpec(name="fluid", dynamic_viscosity=mu))
            except ValueError:
                pass

        # Turbulence model
        if "KIND_TURB_MODEL" in p:
            mapped = self.rules.map_turbulence(p["KIND_TURB_MODEL"])
            self.setup.turbulence_model = mapped if mapped else p["KIND_TURB_MODEL"].lower()
        elif "TURB_MODEL" in p:
            mapped = self.rules.map_turbulence(p["TURB_MODEL"])
            self.setup.turbulence_model = mapped if mapped else p["TURB_MODEL"].lower()
        else:
            self.setup.turbulence_model = "laminar"

        # Numerical schemes
        if "NUM_METHOD_GRAD" in p:
            mapped = self.rules.map_scheme(p["NUM_METHOD_GRAD"])
            self.setup.numerics.gradient_operator = mapped if mapped else "green_gauss_cell"
        if "CONV_NUM_METHOD_FLOW" in p:
            self.setup.numerics.momentum_scheme = _canonical_convection_scheme(p)
            if not self.setup.numerics.momentum_scheme:
                self.setup.numerics.unmapped_settings.append(
                    f"convection: {p['CONV_NUM_METHOD_FLOW'].strip()!r} "
                    f"(no faithful CFDX convection registry equivalent)"
                )
        if "LINEAR_SOLVER" in p:
            self.setup.numerics.linear_solver = p["LINEAR_SOLVER"].strip().lower()
        if "TIME_DISCRE_FLOW" in p:
            td = p["TIME_DISCRE_FLOW"].lower()
            self.setup.numerics.transient_scheme = td
            if "runge" in td:
                self.setup.transient = True
                self.setup.solver_mode = "transient"
        if "ITER" in p:
            try:
                self.setup.numerics.max_iterations = int(float(p["ITER"]))
            except ValueError:
                pass
        if "CFL_NUMBER" in p:
            try:
                self.setup.numerics.residual_target = f"1e-{float(p['CFL_NUMBER'])}"
            except ValueError:
                pass

        # Materials
        if "FREESTREAM_DENSITY" in p:
            try:
                mat = MaterialSpec(name="fluid", density=float(p["FREESTREAM_DENSITY"]))
                if "VISC" in p:
                    mat.dynamic_viscosity = float(p["VISC"])
                self.setup.materials.append(mat)
            except ValueError:
                pass

        self.setup.numerics.raw_settings = dict(p)
        return True

    def parse_solution_csv(self, csv_file: str) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
        """Parse SU2 solution.csv file and return scalar and vector fields.

        Returns:
            tuple of (scalar_fields, vec_fields) where each is a dict mapping
            field name to numpy array of values per cell.
        """
        scalar_fields: dict[str, np.ndarray] = {}
        vec_fields: dict[str, np.ndarray] = {}

        path = Path(csv_file)
        if not path.exists():
            return scalar_fields, vec_fields

        try:
            with open(csv_file, "r") as f:
                reader = csv.reader(f)
                header = next(reader, None)
                if not header:
                    return scalar_fields, vec_fields

                # Normalize header names
                header = [h.strip() for h in header]

                # Identify column indices for known fields
                col_indices: dict[str, int] = {}
                for i, col in enumerate(header):
                    col_indices[col] = i

                # Read all data rows
                data_rows: list[list[float]] = []
                for row in reader:
                    if not row:
                        continue
                    try:
                        data_rows.append([float(v) for v in row])
                    except ValueError:
                        continue

                if not data_rows:
                    return scalar_fields, vec_fields

                data = np.array(data_rows, dtype=np.float64)

                # Map scalar fields (primitive variables)
                scalar_mapping = {
                    "Density": "Density",
                    "Pressure": "Pressure",
                    "Temperature": "Temperature",
                    "Mach": "Mach",
                    "Pressure_Coeff": "Pressure_Coeff",
                    "Laminar_Viscosity": "Laminar_Viscosity",
                    "Eddy_Viscosity": "Eddy_Viscosity",
                    "Turbulent_Kinetic_Energy": "Turbulent_Kinetic_Energy",
                    "Dissipation_Rate": "Dissipation_Rate",
                }

                for su2_name, cfdx_name in scalar_mapping.items():
                    if su2_name in col_indices:
                        col = col_indices[su2_name]
                        scalar_fields[cfdx_name] = data[:, col]

                # Map vector fields (velocity components)
                vel_components = ["Velocity_x", "Velocity_y", "Velocity_z"]
                vel_alt = ["Vel_x", "Vel_y", "Vel_z"]
                vel_su2 = ["VEL_X", "VEL_Y", "VEL_Z"]
                vel_su2_alt = ["Momentum_x", "Momentum_y", "Momentum_z"]

                vel_cols = None
                for comp_set in [vel_components, vel_alt, vel_su2, vel_su2_alt]:
                    if all(c in col_indices for c in comp_set):
                        vel_cols = [col_indices[c] for c in comp_set]
                        break

                if vel_cols is not None:
                    vel_data = data[:, vel_cols]
                    vec_fields["Velocity"] = vel_data

                # Also check for momentum (rho * velocity) - convert to velocity if density available
                if "Velocity" not in vec_fields:
                    mom_cols = None
                    for comp_set in [["Momentum_x", "Momentum_y", "Momentum_z"],
                                      ["MOMENTUM-X", "MOMENTUM-Y", "MOMENTUM-Z"]]:
                        if all(c in col_indices for c in comp_set):
                            mom_cols = [col_indices[c] for c in comp_set]
                            break
                    if mom_cols is not None and "Density" in scalar_fields:
                        mom_data = data[:, mom_cols]
                        rho = scalar_fields["Density"]
                        # Avoid division by zero
                        with np.errstate(divide="ignore", invalid="ignore"):
                            vel_data = mom_data / rho[:, np.newaxis]
                        vel_data = np.nan_to_num(vel_data, nan=0.0, posinf=0.0, neginf=0.0)
                        vec_fields["Velocity"] = vel_data

        except Exception:
            # If parsing fails, return empty dicts (non-blocking)
            pass

        return scalar_fields, vec_fields

    def import_results(self, results_path: str, result: ConversionResult) -> bool:
        """Import SU2 results (solution.csv) into an existing ConversionResult.

        This can be called after convert() to add solution fields.
        """
        base = Path(results_path)
        if base.is_dir():
            solution_file = base / "solution.csv"
        else:
            solution_file = base

        if not solution_file.exists():
            result.gap_report.unsupported_nonblocking(
                "results", "solution_csv",
                f"SU2 solution file not found: {solution_file}",
                "Provide solution.csv or run SU2 to generate results",
            )
            return False

        scalar_fields, vec_fields = self.parse_solution_csv(str(solution_file))

        if not scalar_fields and not vec_fields:
            result.gap_report.unsupported_nonblocking(
                "results", "solution_csv",
                f"No valid fields parsed from SU2 solution.csv: {solution_file}",
                "Check that the CSV has proper headers and numeric data",
            )
            return False

        n_cells_expected = self.setup.mesh_info.n_cells
        mismatches = []
        for name, arr in scalar_fields.items():
            if arr.shape[0] != n_cells_expected:
                mismatches.append(f"Scalar field '{name}' has {arr.shape[0]} values, expected {n_cells_expected} (n_cells)")
            else:
                result.scalar_fields.append((name, arr))

        for name, arr in vec_fields.items():
            if arr.shape[0] != n_cells_expected:
                mismatches.append(f"Vector field '{name}' has {arr.shape[0]} values, expected {n_cells_expected} (n_cells)")
            else:
                result.vec_fields.append((name, arr))

        if mismatches:
            result.gap_report.unsupported_blocking(
                "results", "field_size_mismatch",
                "; ".join(mismatches),
                "Provide a solution file whose field localisation and cardinality match the imported SU2 mesh",
            )
            return False

        result.gap_report.supported(
            "results", "solution_csv",
            f"Imported SU2 solution.csv: {solution_file} ({len(result.scalar_fields)} scalars, {len(result.vec_fields)} vectors)"
        )
        result.field_data_source = str(solution_file)
        return True

    def _populate_gap_report(self, gap: GapAnalysis) -> None:
        gap.supported("mesh", "points", f"Parsed {len(self._points)} vertices from .su2")
        gap.supported("mesh", "elements", f"Parsed {len(self._elements)} elements from .su2")

        for bnd in self._boundaries:
            gap.supported("boundary", bnd["name"], f"Detected boundary marker '{bnd['name']}' with {bnd['n_elements']} faces")

        for bc in self.setup.boundary_conditions:
            gap.supported(
                "boundary_condition", bc.patch_name,
                f"Mapped SU2 marker '{bc.source_zone_name}' to CFDX BC type '{bc.type.value}'",
            )

        if self._cfg_params:
            gap.supported("config", "su2_cfg", f"Parsed {len(self._cfg_params)} config parameters from .cfg")

        unmapped_keys = [k for k, v in self._cfg_params.items()
                         if k not in ("MACH_NUMBER", "AOA", "FREESTREAM_PRESSURE",
                                      "FREESTREAM_TEMPERATURE", "SOLVER", "KIND_TURB_MODEL",
                                      "TURB_MODEL", "NUM_METHOD_GRAD", "CONV_NUM_METHOD_FLOW",
                                      "TIME_DISCRE_FLOW", "MUSCL", "LINEAR_SOLVER", "ITER", "CFL_NUMBER",
                                      "FREESTREAM_DENSITY", "VISC")]
        for k in unmapped_keys:
            gap.unsupported_nonblocking(
                "config", k,
                f"Unmapped SU2 config parameter: {k} = {self._cfg_params[k]}",
                "Manually verify the parameter is applied in CFDX setup",
            )

        if not self._cfg_params:
            gap.unavailable("config", "su2_cfg", "No .cfg config file found — physics and BCs not extracted")

        gap.unsupported_nonblocking(
            "features", "multiphase",
            "SU2 multiphase (VOF) extensions not supported",
            "Use single-phase SU2 configuration",
        )
        gap.unsupported_nonblocking(
            "features", "turbomachinery",
            "SU2 turbomachinery extensions not supported",
            "Use base SU2 configuration without turbo machinery options",
        )
        gap.unsupported_nonblocking(
            "features", "deforming_mesh",
            "SU2 dynamic mesh options not supported",
            "Use static mesh for CFDX conversion",
        )


def _to_meshio_type(su2_type: str) -> str:
    mapping = {
        "TRIANGLE": "triangle",
        "QUADRILATERAL": "quad",
        "TETRAHEDRON": "tetra",
        "HEXAHEDRON": "hexahedron",
        "PRISM": "wedge",
        "PYRAMID": "pyramid",
        "LINE": "line",
        "NODE": "vertex",
    }
    return mapping.get(su2_type, "")
