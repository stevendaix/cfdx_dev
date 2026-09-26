"""Fluent legacy ASCII .cas/.dat adapter for CFDX I/O.

Parses the S-expression format used in Fluent legacy ASCII case files:
  - Section 1:  Header (string tags)
  - Section 2:  Nodes (coordinates)
  - Section 3:  Face nodes (connectivity, variable-length per face)
  - Section 4:  Faces (owner/neighbour cell per face)
  - Section 5:  Cells (cell-to-node connectivity)
  - Section 31: Boundary conditions (zone + BC type)
  - Section 39: Materials (name, density, viscosity, etc.)
  - Section 45: Cell zone conditions
  - Section 48: Model settings (turbulence, transient, schemes)

S-expr format: top-level `(id n_values ...data...)` groups; data may
span multiple lines. Nested parens are tracked for correctness.

Also supports HDF5-based .cas.h5 files via pyfluent's CaseFile reader
(no Fluent license required — uses h5py directly for mesh extraction).

Corresponds to the C++ header:
  src/cfdx/io/fluent/fluent_adapter.h
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union

import numpy as np

from cfdx.io.interfaces import (
    SolverAdapter,
    ConversionResult,
    SourceInfo,
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
)
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.gap_analysis import GapAnalysis


def _try_pyfluent() -> bool:
    """Check if pyfluent CaseFile reader is available (no license needed for .cas.h5)."""
    try:
        from ansys.fluent.core.filereader.case_file import CaseFile  # noqa: F401
        return True
    except ImportError:
        return False


def _infer_turbulence_model(rp_vars: dict) -> str:
    """Infer Fluent turbulence model from RP variables."""
    if rp_vars.get("rp-sa?", False):
        return "sa"
    if rp_vars.get("rp-ke?", False):
        if rp_vars.get("ke-realizable?", False):
            return "ke-realizable-viscous"
        if rp_vars.get("rng-ke-on?", False):
            return "rng-ke"
        return "ke-standard-viscous"
    if rp_vars.get("rp-kw?", False):
        if rp_vars.get("kw-sst-on?", False):
            return "k-omega-sst"
        if rp_vars.get("kw-bsl-on?", False):
            return "bsl-ke"
        return "kw-standard-viscous"
    if rp_vars.get("rp-les?", False):
        return "les"
    return "laminar" if rp_vars.get("rp-lam?", False) else "laminar"


def _infer_zone_type(name: str) -> str:
    """Infer Fluent boundary zone type from surface name."""
    lower = name.lower()
    if "inlet" in lower:
        return "velocity-inlet"
    if "outlet" in lower:
        return "pressure-outlet"
    if "wall" in lower:
        return "wall"
    if "sym" in lower:
        return "symmetry"
    if "interior" in lower:
        return "interior"
    if "periodic" in lower or "periodic" in lower:
        return "periodic"
    return "wall"


# Fluent zone type → (zone_type_string, ZoneType)
_FLUENT_ZONE_TYPES: dict[int, str] = {
    0: "interior",
    1: "velocity-inlet",
    2: "pressure-inlet",
    3: "pressure-outlet",
    4: "wall",
    5: "symmetry",
    6: "periodic",
    7: "exterior",
    8: "internal",
    9: "mass-outflow",
    10: "pressure-far-field",
    11: "outlet-vent",
    12: "overset",
    13: "degenerate",
    14: "slip-wall",
}

# Section name → section ID
_FLUENT_SECTIONS: dict[int, str] = {
    1: "header",
    2: "nodes",
    3: "face_nodes",
    4: "faces",
    5: "cells",
    31: "boundary_zones",
    39: "materials",
    45: "cell_zone_conditions",
    48: "model_settings",
    64: "custom_functions",
}


class FluentAdapter(SolverAdapter):
    solver_name = "Fluent"
    format_name = "cas_dat"

    def __init__(self, rules: Optional[MappingRules] = None) -> None:
        self.rules = rules or MappingRules.load_default()
        self.setup = CaseSetup(source=SourceInfo(solver="Fluent", format="cas_dat"))
        self._points: list[list[float]] = []
        self._face_nodes: list[list[int]] = []
        self._face_owner: list[int] = []
        self._face_neighbour: list[int] = []
        self._cell_nodes: list[list[int]] = []
        self._zones: list[dict] = []
        self._materials: list[dict] = []
        self._model_settings: list[str] = []

    def detect_source(self, case_path: str, info: SourceInfo) -> bool:
        p = Path(str(case_path))
        if p.name.endswith(".cas.h5") or p.name.endswith(".cas_h5"):
            if not _try_pyfluent():
                return False
            info.solver = "Fluent"
            info.format = "cas_h5"
            info.case_path = str(p.parent)
            info.case_name = p.stem.replace(".cas", "")
            info.version = "hdf5"
            return True
        if p.is_dir():
            h5_candidates = list(p.glob("*.cas.h5")) + list(p.glob("*.cas_h5"))
            for h5 in h5_candidates:
                if _try_pyfluent():
                    info.solver = "Fluent"
                    info.format = "cas_h5"
                    info.case_path = str(h5.parent)
                    info.case_name = h5.stem.replace(".cas", "")
                    info.version = "hdf5"
                    return True
        cas_file = p
        if cas_file.suffix != ".cas":
            cas_file = p.with_suffix(".cas")
        if not cas_file.exists() and p.suffix == ".cas":
            cas_file = p
        if not cas_file.exists():
            # Check for .cas.h5 fallback
            h5_candidates = list(p.parent.glob("*.cas.h5")) if p.is_dir() else [p.with_suffix(".cas.h5")]
            for h5 in h5_candidates:
                if h5.exists() and _try_pyfluent():
                    info.solver = "Fluent"
                    info.format = "cas_h5"
                    info.case_path = str(h5.parent)
                    info.case_name = h5.stem.replace(".cas", "")
                    info.version = "hdf5"
                    return True
            return False
        try:
            first = cas_file.read_text(errors="replace")[:1024]
        except OSError:
            return False
        if first.startswith("(") or "COMPILED" in first or "cgns" in first.lower():
            info.solver = "Fluent"
            info.format = "cas_dat"
            info.case_path = str(cas_file.parent)
            info.case_name = cas_file.stem
            info.version = "legacy_ascii"
            return True
        return False

    def convert(self, case_path: Union[str, Path], result: ConversionResult) -> bool:
        base = Path(str(case_path))
        cas_file = base
        dat_file = base
        h5_file: Optional[Path] = None

        if base.is_file():
            cas_file = base
            if base.name.endswith(".cas.h5") or base.name.endswith(".cas_h5"):
                h5_file = base
            else:
                dat_file = base.with_suffix(".dat")
        else:
            if base.is_dir():
                cas_candidates = list(base.glob("*.cas"))
                if cas_candidates:
                    cas_file = cas_candidates[0]
                h5_candidates = list(base.glob("*.cas.h5"))
                if h5_candidates:
                    h5_file = h5_candidates[0]
                dat_candidates = list(base.glob("*.dat"))
                if dat_candidates:
                    dat_file = dat_candidates[0]
            else:
                cas_file = base.with_suffix(".cas")
                dat_file = base.with_suffix(".dat")
                h5_candidate = base.with_suffix(".cas.h5")
                if h5_candidate.exists():
                    h5_file = h5_candidate

        if h5_file is not None:
            mesh_ok = self.parse_cas_h5(str(h5_file))
            if not mesh_ok:
                result.gap_report.unsupported_blocking(
                    "mesh", "fluent_cas_h5",
                    f"Cannot parse Fluent .cas.h5 file: {h5_file}",
                    "pyfluent CaseFile is available; see gap report for details",
                )
                return False
        else:
            mesh_ok = self.parse_cas(str(cas_file)) if cas_file.exists() else False
            if not mesh_ok:
                result.gap_report.unsupported_blocking(
                    "mesh", "fluent_cas",
                    f"Cannot parse Fluent .cas file: {cas_file}",
                    "Ensure the file is a legacy ASCII Fluent .cas (not binary .cas.h5)",
                )
                return False

        # Attempt results import — check for .dat.h5 alongside .cas.h5
        if h5_file is not None:
            h5_dat = h5_file.parent / (h5_file.stem.replace(".cas", ".dat.h5") if h5_file.stem.endswith(".cas") else h5_file.parent.name + ".dat.h5")
            if h5_dat.exists():
                self.import_results(str(h5_dat), result)
            else:
                # Also try same basename with .dat.h5
                h5_dat = h5_file.with_suffix(".dat.h5")
                if h5_dat.exists():
                    self.import_results(str(h5_dat), result)
                else:
                    result.gap_report.unsupported_nonblocking(
                        "results", "fluent_dat_h5",
                        "No .dat.h5 solution file found alongside .cas.h5",
                        "Copy the .dat.h5 file alongside the .cas.h5 for solution data import",
                    )
        elif dat_file.exists():
            self.import_results(str(dat_file), result)

        result.source = self.setup.source
        result.setup = self.setup.model_copy(deep=True)

        # Build mesh dict — handle mixed cell types (tetra/hexa/prism/etc.)
        points_arr = np.array(self._points, dtype=np.float64) if self._points else np.empty((0, 3))
        cells: dict[str, np.ndarray] = {}
        if self._cell_nodes:
            type_map = {4: "tetra", 5: "pyramid", 6: "wedge", 8: "hexahedron"}
            by_type: dict[str, list[list[int]]] = {}
            skipped = 0
            for cell in self._cell_nodes:
                n = len(cell)
                key = type_map.get(n)
                if key:
                    by_type.setdefault(key, []).append(cell)
                else:
                    skipped += 1
            for key, cell_list in by_type.items():
                cells[key] = np.array(cell_list, dtype=np.int32)
            if skipped > 0:
                result.gap_report.approximated(
                    "mesh", "polyhedral_cells",
                    f"Skipped {skipped} polyhedral/unknown cells in .cas.h5",
                    "Use Fluent export to VTM or CGNS for polyhedral mesh support",
                )

        mesh = {"points": points_arr, "cells": cells, "cell_data": {}}

        if self._points:
            result.mesh = mesh
            result.setup.mesh_info.n_vertices = len(self._points)
            result.setup.mesh_info.n_cells = len(self._cell_nodes)
            result.setup.mesh_info.n_faces = len(self._face_nodes)
            result.setup.mesh_info.dimension = 3 if any(len(p) >= 3 for p in self._points) else 2

        self._populate_gap_report(result.gap_report)
        return not result.gap_report.has_blocking()

    def parse_cas(self, cas_file: str) -> bool:
        text = Path(cas_file).read_text(errors="replace")
        sections = _parse_s_expressions(text)

        for section_id, data_tokens in sections:
            if section_id == 1:
                self.setup.source_metadata = {"header_tags": list(data_tokens)}
            elif section_id == 2:
                # Nodes: header is (n_nodes n_dims), then coords
                n_nodes = int(data_tokens[0])
                coords = [float(t) for t in data_tokens[2:]]
                dim = 3  # Fluent stores 3D coords (x y z) per node line
                for i in range(n_nodes):
                    base = i * dim
                    pt = coords[base : base + dim]
                    while len(pt) < 3:
                        pt.append(0.0)
                    self._points.append(pt)
            elif section_id == 3:
                # Face nodes: header is (n_faces n_max), then var-length face records
                i = 2  # skip header
                while i < len(data_tokens):
                    n_fn = int(data_tokens[i])
                    i += 1
                    face = [int(data_tokens[i + j]) - 1 for j in range(n_fn)]
                    i += n_fn
                    self._face_nodes.append(face)
            elif section_id == 4:
                # Faces: header is (n_faces n_max), then var-length face data
                i = 2  # skip header
                while i < len(data_tokens):
                    n_cells = int(data_tokens[i])
                    i += 1
                    if i < len(data_tokens):
                        owner = int(data_tokens[i]) - 1
                        i += 1
                    else:
                        owner = 0
                    if n_cells >= 2 and i < len(data_tokens):
                        neighbour = int(data_tokens[i]) - 1
                        i += 1
                    else:
                        neighbour = -1
                    self._face_owner.append(owner)
                    self._face_neighbour.append(neighbour)
            elif section_id == 5:
                # Cells: header is (n_cells n_total), then var-length cell records
                n_cells = int(data_tokens[0])
                i = 2  # skip header
                for _ in range(n_cells):
                    if i >= len(data_tokens):
                        break
                    n_cn = int(data_tokens[i])
                    i += 1
                    cell = [int(data_tokens[i + j]) - 1 for j in range(n_cn)]
                    i += n_cn
                    self._cell_nodes.append(cell)
            elif section_id == 31:
                # Boundary zones: header is (n_zones n_fields_per_zone)
                n_zones = int(data_tokens[0])
                fields_per_zone = int(data_tokens[1])
                i = 2
                for _ in range(n_zones):
                    if i + fields_per_zone > len(data_tokens):
                        break
                    zone_id = int(data_tokens[i])
                    zone_name = data_tokens[i + 1].strip('"')
                    node_type = data_tokens[i + 2]
                    bc_type_str = data_tokens[i + 3].strip('"')
                    n_elems = int(data_tokens[i + 4])
                    zone_type_id = int(data_tokens[i + 5])
                    # Remaining fields are zone-specific surface info
                    self._zones.append({
                        "zone_id": zone_id,
                        "name": zone_name,
                        "node_type": node_type,
                        "element_type": bc_type_str,
                        "n_elements": n_elems,
                        "zone_type": bc_type_str,
                        "zone_type_id": zone_type_id,
                    })
                    i += fields_per_zone
            elif section_id == 39:
                # Materials: header is (n_materials n_fields_per_mat)
                n_mats = int(data_tokens[0])
                i = 2  # skip header
                for _ in range(n_mats):
                    if i + 6 > len(data_tokens):
                        break
                    mat_name = data_tokens[i].strip('"')
                    mat_type = data_tokens[i + 1].strip('"')
                    density = float(data_tokens[i + 2])
                    viscosity = float(data_tokens[i + 3])
                    cond = float(data_tokens[i + 4])
                    cp = float(data_tokens[i + 5])
                    mw = float(data_tokens[i + 6])
                    i += 7
                    self._materials.append({
                        "name": mat_name,
                        "type": mat_type,
                        "density": density,
                        "viscosity": viscosity,
                        "conductivity": cond,
                        "cp": cp,
                        "mw": mw,
                    })
            elif section_id == 45:
                # Cell zone conditions: header is (n_zones n_fields)
                n_zones = int(data_tokens[0])
                fields_per_zone = int(data_tokens[1])
                i = 2
                for _ in range(n_zones):
                    if i + fields_per_zone > len(data_tokens):
                        break
                    zone_id = int(data_tokens[i])
                    zone_name = data_tokens[i + 1].strip('"')
                    node_type = data_tokens[i + 2]
                    elem_type = data_tokens[i + 3].strip('"')
                    n_elems = int(data_tokens[i + 4])
                    self._zones.append({
                        "zone_id": zone_id,
                        "name": zone_name,
                        "node_type": node_type,
                        "element_type": elem_type,
                        "n_elements": n_elems,
                        "zone_type": "cell",
                    })
                    i += fields_per_zone
            elif section_id == 48:
                # Model settings: skip header (2 tokens), rest are key-value pairs
                self._model_settings = list(data_tokens[2:])

        if not self._points:
            return False

        self.setup.source.solver = "Fluent"
        self.setup.source.format = "cas_dat"
        self.setup.source.case_path = cas_file
        self.setup.source.case_name = Path(cas_file).stem
        self.setup.source.version = "legacy_ascii"

        # Build materials
        for mat in self._materials:
            self.setup.materials.append(MaterialSpec(
                name=mat["name"],
                density=mat["density"],
                dynamic_viscosity=mat["viscosity"],
                thermal_conductivity=mat["conductivity"],
                specific_heat=mat["cp"],
                molecular_weight=mat["mw"],
                eos_model="ideal_gas" if mat["type"] == "constant" else "real_gas",
            ))

        # Build boundary conditions from zones
        for zone in self._zones:
            if zone.get("node_type") == "face":
                bc = BoundarySpec(
                    patch_name=zone["name"],
                    source_zone_name=zone["name"],
                    source_type_name=zone.get("zone_type", "unknown"),
                )
                mapped = self.rules.map_boundary_type(zone.get("zone_type", ""))
                bc.type = mapped
                bc.value_type = self.rules.map_bc_value_type(zone.get("zone_type", ""))
                self.setup.boundary_conditions.append(bc)

        # Parse model settings
        _parse_model_settings(self._model_settings, self.setup, self.rules)

        self.setup.mesh_info.n_vertices = len(self._points)
        self.setup.mesh_info.n_cells = len(self._cell_nodes)
        self.setup.mesh_info.n_faces = len(self._face_nodes)

        return True

    def parse_cas_h5(self, cas_h5_file: str) -> bool:
        """Parse Fluent HDF5 case file (.cas.h5) using pyfluent CaseFile.

        Uses pyfluent's offline CaseFile reader — no Fluent license required.
        Reads mesh data (points, connectivity, zones) directly via h5py.
        """
        if not _try_pyfluent():
            self.setup.source.version = "pyfluent_not_available"
            return False

        try:
            from ansys.fluent.core.filereader.case_file import CaseFile
        except ImportError:
            return False

        try:
            reader = CaseFile(case_file_name=cas_h5_file)
        except Exception as e:
            self.setup.source.version = f"pyfluent_error: {e}"
            return False

        mesh = reader.get_mesh()
        surface_ids = list(mesh.get_surface_ids())
        surface_names = mesh.get_surface_names()

        volume_sid = None
        for sid in surface_ids:
            locs = mesh.get_surface_locs(int(sid))
            v = mesh.get_vertices(int(sid))
            n_points = len(v) // 3
            if n_points > 0 and volume_sid is None:
                volume_sid = int(sid)
            elif n_points > (len(mesh.get_vertices(int(volume_sid))) // 3 if volume_sid else 0):
                volume_sid = int(sid)

        if volume_sid is None:
            return False

        v = mesh.get_vertices(volume_sid)
        points = v.reshape(-1, 3)
        points = np.ascontiguousarray(points, dtype=np.float64)
        self._points = points.tolist()

        conn = mesh.get_connectivity(volume_sid)
        i = 0
        while i < len(conn):
            n_cn = int(conn[i])
            cell = [int(conn[i + j]) - 1 for j in range(1, n_cn + 1)]
            i += 1 + n_cn
            self._cell_nodes.append(cell)

        for idx, sid in enumerate(surface_ids):
            sid_int = int(sid)
            if sid_int == volume_sid:
                continue
            name = surface_names[idx]
            locs = mesh.get_surface_locs(sid_int)
            face_conn = mesh.get_connectivity(sid_int)
            n_faces = 0
            j = 0
            while j < len(face_conn):
                n_fn = int(face_conn[j])
                j += 1 + n_fn
                n_faces += 1
            self._zones.append({
                "zone_id": sid_int,
                "name": name,
                "node_type": "face",
                "element_type": name,
                "n_elements": n_faces,
                "zone_type": _infer_zone_type(name),
                "zone_type_id": 4,
            })

        self.setup.source.solver = "Fluent"
        self.setup.source.format = "cas_h5"
        self.setup.source.case_path = cas_h5_file
        self.setup.source.case_name = Path(cas_h5_file).stem.replace(".cas", "")
        self.setup.source.version = "hdf5"

        # Extract turbulence model and physics from RP variables
        rp_vars = reader.rp_vars()
        turb_model = _infer_turbulence_model(rp_vars)
        if turb_model:
            mapped = self.rules.map_turbulence(turb_model)
            self.setup.turbulence_model = mapped if mapped else turb_model

        if rp_vars.get("pseudo-transient-formulation?", False):
            self.setup.transient = True
            self.setup.solver_mode = "transient"
        else:
            self.setup.transient = False
            self.setup.solver_mode = "steady"

        # Extract materials from RP variables
        _parse_rp_materials(rp_vars, self.setup)

        # Extract reference values, gravity, solver settings from RP vars
        _extract_cas_h5_reference_values(rp_vars, self.setup)

        # Extract boundary conditions from surface zones + BC RP vars
        _parse_rp_boundary_conditions(rp_vars, self.setup, self._zones)

        self.setup.mesh_info.n_vertices = len(self._points)
        self.setup.mesh_info.n_cells = len(self._cell_nodes)
        self.setup.mesh_info.n_faces = sum(z.get("n_elements", 0) for z in self._zones)
        self.setup.mesh_info.dimension = 3

        return True

    def import_results(self, dat_file: str, result: ConversionResult) -> bool:
        """Import Fluent .dat.h5 solution data via pyfluent FileSession (no license needed)."""
        if not _try_pyfluent():
            result.gap_report.unsupported_nonblocking(
                "results", "fluent_dat_h5",
                "pyfluent not installed; cannot read .dat.h5 solution data",
                "Install ansys-fluent-core for offline result import",
            )
            return False

        try:
            from ansys.fluent.core.file_session import FileSession
        except ImportError:
            return False

        try:
            session = FileSession()
            case_file = self.setup.source.case_path
            session.read_case(case_file)
            session.read_data(dat_file)
        except Exception as e:
            result.gap_report.unsupported_nonblocking(
                "results", "fluent_dat_h5",
                f"Cannot read Fluent .dat.h5 file: {e}",
                "Ensure both .cas.h5 and .dat.h5 are from the same case",
            )
            return False

        # Extract available scalar fields
        field_info = session.fields.field_info
        scalar_fields = field_info.get_scalar_fields_info()
        vector_fields = field_info.get_vector_fields_info()

        result.available_fields = {
            "scalar": list(scalar_fields.keys()),
            "vector": list(vector_fields.keys()),
        }
        result.field_data_source = dat_file
        result.gap_report.supported(
            "results", "fluent_dat_h5",
            f"Extracted {len(scalar_fields)} scalar and {len(vector_fields)} vector fields via pyfluent FileSession",
        )

        # Map common Fluent variable names to CFDX field names
        field_map = {
            "SV_P": "pressure",
            "SV_T": "temperature",
            "SV_U": "velocity_x",
            "SV_V": "velocity_y",
            "SV_W": "velocity_z",
            "SV_K": "turbulent_kinetic_energy",
            "SV_O": "turbulent_dissipation_rate",
        }

        for fluent_name, cfdx_name in field_map.items():
            if fluent_name in scalar_fields:
                result.gap_report.supported(
                    "results", cfdx_name,
                    f"Field '{fluent_name}' available via FileSession (use field_data.get_scalar_field_data)",
                )

        return True

    def _populate_gap_report(self, gap: GapAnalysis) -> None:
        source_version = self.setup.source.version

        if source_version == "legacy_ascii" or source_version == "cas_dat":
            self._populate_gap_report_legacy(gap)
        elif source_version == "hdf5" or source_version == "cas_h5":
            self._populate_gap_report_hdf5(gap)
        else:
            gap.unsupported_blocking(
                "mesh", "unknown",
                f"Unknown source version: {source_version}",
                "Use a valid Fluent case file or install pyfluent for .cas.h5 support",
            )

    def _populate_gap_report_legacy(self, gap: GapAnalysis) -> None:
        gap.supported("mesh", "nodes", f"Parsed {len(self._points)} vertices from .cas")
        gap.supported("mesh", "face_nodes", f"Parsed {len(self._face_nodes)} face-node connectivities from .cas")
        gap.supported("mesh", "cells", f"Parsed {len(self._cell_nodes)} cells from .cas")

        for zone in self._zones:
            gap.supported("boundary_zone", zone.get("name", "?"), f"Detected zone '{zone.get('name', '?')}' type={zone.get('zone_type', 'unknown')}")

        for mat in self._materials:
            gap.supported("material", mat["name"], f"Parsed material '{mat['name']}' (ρ={mat['density']}, μ={mat['viscosity']})")

        for bc in self.setup.boundary_conditions:
            gap.supported("boundary_condition", bc.patch_name, f"Mapped Fluent zone '{bc.source_zone_name}' to {bc.type.value}")

        gap.unsupported_nonblocking(
            "features", "binary_cas",
            "Binary .cas.h5 format not supported (only legacy ASCII .cas)",
            "Export as legacy ASCII: File → Export → Case...",
        )
        gap.unsupported_nonblocking(
            "features", "dpm",
            "Discrete Phase Model (DPM) particles not supported",
            "Use single-phase flow for CFDX conversion",
        )
        gap.unsupported_nonblocking(
            "features", "udfs",
            "UDFs and custom field functions not supported",
            "Pre-compute and export field values via Fluent surface or volume monitor",
        )
        gap.unsupported_nonblocking(
            "features", "multiphase",
            "Multiphase models (VOF/Mixture/Eulerian) not supported",
            "Use single-phase configuration for CFDX conversion",
        )
        gap.unsupported_nonblocking(
            "features", "sliding_mesh",
            "Sliding mesh and dynamic mesh not supported",
            "Use steady-state or static mesh for CFDX conversion",
        )
        gap.unsupported_nonblocking(
            "features", "results",
            "Fluent .dat results import requires meshio or neutral export; not parsed natively",
            "Export results as .cgns, .vtk, or .plt via Fluent File → Export",
        )

    def _populate_gap_report_hdf5(self, gap: GapAnalysis) -> None:
        gap.supported("mesh", "nodes", f"Parsed {len(self._points)} vertices from .cas.h5 via pyfluent CaseFile")
        gap.supported("mesh", "cells", f"Parsed {len(self._cell_nodes)} cells from .cas.h5")
        for zone in self._zones:
            gap.supported("boundary_zone", zone.get("name", "?"), f"Detected zone '{zone.get('name', '?')}' type={zone.get('zone_type', 'unknown')}")
        for mat in self.setup.materials:
            gap.supported("material", mat.name, f"Parsed material '{mat.name}' (ρ={mat.density}, μ={mat.dynamic_viscosity})")
        for bc in self.setup.boundary_conditions:
            gap.supported("boundary_condition", bc.patch_name, f"Mapped Fluent zone '{bc.source_zone_name}' to {bc.type.value}")
        gap.supported("mesh", "pyfluent_reader", "pyfluent CaseFile reader provided offline access to HDF5 mesh data")
        gap.supported("physics", "turbulence", f"Inferred turbulence model: {self.setup.turbulence_model}")
        gap.supported("physics", "reference_values", "Extracted reference length, density, velocity, temperature, pressure")
        gap.approximated(
            "features", "dpm",
            "Discrete Phase Model (DPM) particles not extracted from .cas.h5 (CaseFile mesh-only reader)",
            "Export DPM data separately via Fluent surface monitors or .cgns/.vtk export",
        )
        gap.approximated(
            "features", "udfs",
            "UDFs and custom field functions not available via CaseFile (mesh-only reader)",
            "Pre-compute and export field values via Fluent surface or volume monitor",
        )
        gap.unsupported_nonblocking(
            "features", "solution_data",
            "CaseFile reader provides mesh only; no .dat.h5 solution data extracted",
            "Use FileSession or export results as .cgns, .vtk, or .plt via Fluent File → Export",
        )
        gap.unsupported_nonblocking(
            "features", "multiphase",
            "Multiphase models (VOF/Mixture/Eulerian) not supported via CaseFile",
            "Use single-phase configuration for CFDX conversion",
        )
        gap.unsupported_nonblocking(
            "features", "sliding_mesh",
            "Sliding mesh and dynamic mesh not supported via CaseFile",
            "Use steady-state or static mesh for CFDX conversion",
        )
        gap.unsupported_nonblocking(
            "features", "results",
            "Fluent .dat.h5 results import requires FileSession (license needed); not parsed natively",
            "Export results as .cgns, .vtk, or .plt via Fluent File → Export",
        )


def _parse_s_expressions(text: str) -> list[tuple[int, list[str]]]:
    """Parse top-level S-expressions from Fluent .cas ASCII text.

    Returns a list of (section_id, [tokens]).
    Each section: `(id n_values tok1 tok2 ... )`
    Tokens are split on whitespace; string literals in double quotes are kept
    as single tokens.
    """
    sections: list[tuple[int, list[str]]] = []
    i = 0
    n = len(text)

    while i < n:
        # Find opening paren
        if text[i] == "(":
            depth = 1
            i += 1
            start = i
            tokens: list[str] = []
            while i < n and depth > 0:
                c = text[i]
                if c == "(":
                    depth += 1
                    i += 1
                elif c == ")":
                    depth -= 1
                    if depth == 0:
                        break
                    i += 1
                elif c == '"':
                    # Read quoted string
                    j = i + 1
                    while j < n and text[j] != '"':
                        j += 1
                    tokens.append(f'"{text[i+1:j]}"')
                    i = j + 1
                elif c in (" ", "\t", "\n", "\r"):
                    i += 1
                else:
                    j = i
                    while j < n and text[j] not in "()\t\n\r\" ":
                        j += 1
                    tokens.append(text[i:j])
                    i = j

            if tokens and depth == 0:
                try:
                    sec_id = int(tokens[0])
                    # Second token is n_values; the actual data starts from token[2]
                    # But some sections have n_values as the second token
                    data_tokens = tokens[1:]
                    sections.append((sec_id, data_tokens))
                except ValueError:
                    pass
        else:
            i += 1

    return sections


def _parse_model_settings(settings: list[str], setup: CaseSetup, rules: MappingRules) -> None:
    """Parse model settings tokens from Fluent section 48."""
    i = 0
    while i < len(settings):
        key = settings[i].lower()
        # Look for key-value pairs
        if i + 1 < len(settings):
            val = settings[i + 1]
            if key == "ke-realizable-viscous" or key == "ke-standard-viscous":
                mapped = rules.map_turbulence(key)
                setup.turbulence_model = mapped if mapped else key
            elif key == "transient-formulation":
                setup.solver_mode = val
            elif key == "time-step":
                try:
                    setup.time_step = float(val)
                except ValueError:
                    pass
            elif key == "number-of-time-steps":
                try:
                    setup.max_time_steps = int(float(val))
                except ValueError:
                    pass
            elif key == "end-time":
                try:
                    setup.end_time = float(val)
                except ValueError:
                    pass
            elif key == "simple":
                setup.numerics.coupled_solver = val
            i += 2
        else:
            key = settings[i].lower()
            if key in ("ke-realizable-viscous", "ke-standard-viscous",
                       "kw-viscous", "mixing-length", "transition-to-turbulence"):
                mapped = rules.map_turbulence(key)
                setup.turbulence_model = mapped if mapped else key
            elif key == "simple":
                setup.numerics.coupled_solver = "SIMPLE"
            i += 1


def _parse_rp_materials(rp_vars: dict, setup: CaseSetup) -> None:
    """Extract materials from Fluent RP variables (materials key)."""
    materials_raw = rp_vars.get("materials")
    if not materials_raw or not isinstance(materials_raw, list):
        return

    for entry in materials_raw:
        if not isinstance(entry, list) or len(entry) < 2:
            continue
        name = str(entry[0])
        mat_type = str(entry[1])

        mat = MaterialSpec(name=name)
        for prop in entry[2:]:
            if isinstance(prop, (list, tuple)) and len(prop) >= 2:
                pname = str(prop[0])
                if isinstance(prop[1], tuple) and len(prop[1]) == 2:
                    model, value = prop[1]
                    if model == "constant" and isinstance(value, (int, float)):
                        if pname == "density":
                            mat.density = float(value)
                        elif pname == "viscosity":
                            mat.dynamic_viscosity = float(value)
                        elif pname == "thermal-conductivity":
                            mat.thermal_conductivity = float(value)
                        elif pname == "specific-heat":
                            mat.specific_heat = float(value)
                        elif pname == "molecular-weight":
                            mat.molecular_weight = float(value)
        mat.eos_model = "ideal_gas" if mat_type == "fluid" else "solid"
        setup.materials.append(mat)


def _parse_rp_boundary_conditions(rp_vars: dict, setup: CaseSetup, zones: list[dict]) -> None:
    """Extract boundary condition info from RP variables and zones."""
    for zone in zones:
        bc = BoundarySpec(
            patch_name=zone.get("name", ""),
            source_zone_name=zone.get("name", ""),
            source_type_name=zone.get("zone_type", ""),
        )
        mapped = MappingRules.load_default().map_boundary_type(zone.get("zone_type", ""))
        bc.type = mapped
        bc.value_type = MappingRules.load_default().map_bc_value_type(zone.get("zone_type", ""))
        setup.boundary_conditions.append(bc)


def _extract_cas_h5_reference_values(rp_vars: dict, setup: CaseSetup) -> None:
    """Extract reference values, gravity, and solver settings from RP vars."""
    setup.ref_length = float(rp_vars.get("reference-length", 1.0))
    setup.ref_density = float(rp_vars.get("reference-density", 1.0))
    setup.ref_velocity = float(rp_vars.get("reference-velocity", 1.0))
    setup.ref_temperature = float(rp_vars.get("reference-temperature", 293.15))
    setup.ref_pressure = float(rp_vars.get("reference-pressure", 101325.0))

    if rp_vars.get("gravity?", False):
        gx = float(rp_vars.get("gravity/x", 0.0))
        gy = float(rp_vars.get("gravity/y", 0.0))
        gz = float(rp_vars.get("gravity/z", 0.0))
        setup.gravity_vector = f"[{gx}, {gy}, {gz}]"
    else:
        setup.gravity_vector = "[0, 0, 0]"

    setup.energy_model = "energy" if rp_vars.get("energy?", False) else "isothermal"

    pressure_relax = rp_vars.get("pressure/relax", None)
    if pressure_relax is not None:
        try:
            setup.numerics.under_relaxation_pressure = float(pressure_relax)
        except (TypeError, ValueError):
            pass

    if rp_vars.get("piso/scheme", None) is not None:
        setup.numerics.coupled_solver = "PISO"
    elif rp_vars.get("simplec/skew-iter", None) is not None:
        setup.numerics.coupled_solver = "SIMPLEC"

    if rp_vars.get("piso/coupling?", False):
        setup.numerics.coupled_solver = "PISO"
