"""SU2 (.su2 mesh + .cfg config) adapter for CFDX I/O.

Parses the actual SU2 file format:
  - NDIME= <n>           dimensionality
  - NELEM= <n>           element count, then `<type> <n1> <n2> ... <cell_id>` lines
  - NPOIN= <n>           node count, then `x y z idx` lines
  - NMARK= <n>           boundary marker count
  - MARKER_TAG <name>   marker name
  - MARKER_ELEMS <n>    face count for this marker
    then `<type> <n1> ...` lines

.cfg format: `%`-commented, `KEY = value` pairs.

Corresponds to the C++ header:
  src/cfdx/io/su2/su2_adapter.h
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union

import numpy as np
from pydantic import ValidationError

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
    InitialCondition,
    NumericalScheme,
    MeshMetadata,
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
            solver_val = p["SOLVER"].lower()
            if "euler" in solver_val:
                self.setup.physics_model = "compressible_euler"
                self.setup.energy_model = "compressible"
            elif "navier" in solver_val:
                self.setup.physics_model = "compressible_navier_stokes"
                self.setup.energy_model = "compressible"
            elif "heat" in solver_val:
                self.setup.physics_model = "heat_equation"
                self.setup.energy_model = "buoyant"

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
            scheme_val = p["CONV_NUM_METHOD_FLOW"]
            self.setup.numerics.momentum_scheme = scheme_val.lower()
        if "TIME_DISCRE_FLOW" in p:
            td = p["TIME_DISCRE_FLOW"].lower()
            self.setup.numerics.transient_scheme = td
            if "runge" in td:
                self.setup.transient = True
                self.setup.solver_mode = "transient"
        if "MUSCL" in p:
            if p["MUSCL"].lower() in ("yes", "true", "1"):
                self.setup.numerics.momentum_scheme = "MUSCL"
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
                                      "TIME_DISCRE_FLOW", "MUSCL", "ITER", "CFL_NUMBER",
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
