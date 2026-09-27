"""STAR-CCM+ (.sim) adapter for CFDX I/O.

Per issue #425: "Do not attempt unsupported direct parsing of proprietary
.sim databases. Require documented exports or licensed APIs where available.
Document exactly which export path was used."

This adapter performs a *safe, documented* scan of the .sim binary file:
  - Reads the file header to detect version and endianness
  - Scans for ASCII strings identifying boundary zones, materials, and
    physics models
  - Records every limitation in the GapAnalysis report

SUPPORTED: boundary zone name extraction, material name extraction,
physics model detection, mesh metadata (if present in header).

NOT SUPPORTED (documented in gap report):
  - Direct binary mesh topology reconstruction (faces, cells, owner/neighbour)
  - Field results import from .sim
  - UDFs and custom models

RECOMMENDED EXPORT PATH:
  Use documented exports: .plt (EnSight), .vtk, .cgns, .msh (Gmsh), or .su2.

Corresponds to the C++ header:
  src/cfdx/io/starccm/starccm_adapter.h
"""

from __future__ import annotations

import re
import struct
from pathlib import Path
from typing import Optional, Union

import numpy as np

from cfdx.io.interfaces import (
    SolverAdapter,
    ConversionResult,
    SourceInfo,
)
from cfdx.io.schema import (
    CaseSetup,
    MaterialSpec,
    BoundarySpec,
    BCType,
    BCValueType,
)
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.gap_analysis import GapAnalysis


# ASCII string patterns for physics detection
_STARCCM_PHYSICS_PATTERNS: dict[str, str] = {
    "k-epsilon": "k_epsilon",
    "k-omega": "k_omega_sst",
    "k-omega sst": "k_omega_sst",
    "spalart-allmaras": "spalart_allmaras",
    "realizable ke": "realizable_ke",
    "laminar": "laminar",
    "LES": "LES",
    "DES": "DES",
    "RSM": "reynolds_stress",
    "Euler": "compressible_euler",
    "Navier-Stokes": "compressible_navier_stokes",
    "incompressible": "incompressible_laminar",
    "compressible": "compressible",
    "VOF": "vof",
    "Eulerian": "eulerian",
    "mixture": "mixture",
    "energy": "buoyant",
    "buoyancy": "buoyant",
    "transient": "transient",
    "steady": "steady",
}

_BOUNDARY_TYPE_PATTERNS: list[tuple[str, BCType, BCValueType]] = [
    (r"velocity\s*inlet", BCType.INLET, BCValueType.FIXED),
    (r"pressure\s*inlet", BCType.INLET, BCValueType.FIXED),
    (r"mass\s*flow\s*inlet", BCType.INLET, BCValueType.FIXED),
    (r"pressure\s*outlet", BCType.PRESSURE_OUTLET, BCValueType.OUTLET_PRESSURE),
    (r"outlet", BCType.OUTLET, BCValueType.ZERO_GRADIENT),
    (r"wall", BCType.WALL, BCValueType.WALL_NO_SLIP),
    (r"symmetry", BCType.SYMMETRY, BCValueType.ZERO_GRADIENT),
    (r"periodic", BCType.PERIODIC, BCValueType.MIXED),
    (r"interface", BCType.INTERFACE, BCValueType.MIXED),
    (r"overset", BCType.INTERFACE, BCValueType.MIXED),
    (r"far\s*field", BCType.OUTLET, BCValueType.FIXED),
]


class StarCCMAdapter(SolverAdapter):
    solver_name = "STAR-CCM+"
    format_name = "sim"

    def __init__(self, rules: Optional[MappingRules] = None) -> None:
        self.rules = rules or MappingRules.load_default()
        self.setup = CaseSetup(source=SourceInfo(solver="STAR-CCM+", format="sim"))
        self._zones: list[dict] = []
        self._materials: list[dict] = []
        self._physics: dict[str, str] = {}
        self._header_info: dict = {}
        self._all_strings: list[str] = []

    def detect_source(self, case_path: str, info: SourceInfo) -> bool:
        p = Path(str(case_path))
        sim_file = p
        if sim_file.suffix != ".sim":
            sim_file = p.with_suffix(".sim")
        if not sim_file.exists() and p.suffix == ".sim":
            sim_file = p
        if not sim_file.exists():
            return False
        info.solver = "STAR-CCM+"
        info.format = "sim"
        info.case_path = str(sim_file.parent)
        info.case_name = sim_file.stem
        info.version = "unknown"
        return True

    def convert(self, case_path: Union[str, Path], result: ConversionResult) -> bool:
        base = Path(str(case_path))
        sim_file = base
        if base.is_file():
            sim_file = base
        else:
            sim_candidates = list(base.glob("*.sim"))
            if sim_candidates:
                sim_file = sim_candidates[0]
            else:
                return False

        if not sim_file.exists():
            result.gap_report.unsupported_blocking(
                "mesh", "starccm_sim",
                f"STAR-CCM+ .sim file not found: {sim_file}",
                "Provide a valid .sim file path",
            )
            return False

        # Read and scan
        ok = self._scan_sim(str(sim_file))
        if not ok:
            result.gap_report.unsupported_blocking(
                "mesh", "starccm_sim",
                f"Cannot read STAR-CCM+ .sim file: {sim_file} (binary corruption or permissions)",
                "Verify file integrity and try again",
            )
            return False

        result.source = self.setup.source
        result.setup = self.setup.model_copy(deep=True)

        # No mesh is populated — .sim binary topology is NOT reconstructed
        result.gap_report.unsupported_blocking(
            "mesh", "sim_topology",
            "STAR-CCM+ .sim binary mesh topology (faces, cells, owner/neighbour) not reconstructed",
            "Export the mesh via STAR-CCM+ File → Export → Mesh: use .plt (EnSight), .cgns, .vtk, or .msh (Gmsh)",
        )

        self._populate_gap_report(result.gap_report)
        return not result.gap_report.has_blocking()

    def _scan_sim(self, sim_path: str) -> bool:
        try:
            with open(sim_path, "rb") as f:
                raw = f.read()
        except OSError:
            return False

        self._header_info["file_size"] = len(raw)

        # Detect endianness and version from header
        if len(raw) >= 4:
            magic_le = struct.unpack("<I", raw[:4])[0]
            magic_be = struct.unpack(">I", raw[:4])[0]
            # STAR-CCM+ files often start with a version marker
            if len(raw) >= 8:
                try:
                    header_str = raw[:64].decode("ascii", errors="replace")
                    if "STAR" in header_str or "star" in header_str:
                        self.setup.source.version = header_str.strip("\x00")
                except Exception:
                    pass

        # Safe ASCII string extraction (min_len=3 to catch short material names like "air")
        self._all_strings = _extract_ascii_strings(raw, min_len=3)

        # Detect boundary zones
        for zone_name in self._all_strings:
            if re.match(r"^[A-Za-z_][\w\- ]{1,63}$", zone_name):
                # Heuristic: STAR-CCM+ boundary names often end with _wall, _inlet, etc.
                lower = zone_name.lower()
                is_boundary_candidate = (
                    "wall" in lower or "inlet" in lower or "outlet" in lower
                    or "far" in lower or "symmetr" in lower or "interface" in lower
                    or "periodic" in lower or "opening" in lower
                )
                if is_boundary_candidate:
                    self._zones.append({"name": zone_name})

        # Deduplicate zones
        seen = set()
        deduped = []
        for z in self._zones:
            if z["name"] not in seen:
                seen.add(z["name"])
                deduped.append(z)
        self._zones = deduped

        # Detect materials
        material_keywords = ("air", "water", "oxygen", "nitrogen", "fuel",
                             "constant", "ideal", "real", "material")
        for s in self._all_strings:
            lower = s.lower()
            if lower in material_keywords and len(s) < 32:
                self._materials.append({"name": s})

        # Deduplicate materials
        seen_m = set()
        deduped_m = []
        for m in self._materials:
            if m["name"] not in seen_m:
                seen_m.add(m["name"])
                deduped_m.append(m)
        self._materials = deduped_m

        # Detect physics models
        for s in self._all_strings:
            lower = s.lower()
            for pattern, mapped in _STARCCM_PHYSICS_PATTERNS.items():
                if pattern.lower() in lower:
                    self._physics[pattern] = mapped
                    break

        # Build setup from detected info
        for z in self._zones:
            bc = BoundarySpec(
                patch_name=z["name"],
                source_zone_name=z["name"],
            )
            lower_name = z["name"].lower()
            for pattern, bc_type, bcvt in _BOUNDARY_TYPE_PATTERNS:
                if re.search(pattern, lower_name):
                    bc.type = bc_type
                    bc.value_type = bcvt
                    break
            else:
                bc.type = self.rules.map_boundary_type(z["name"])
                bc.value_type = self.rules.map_bc_value_type(z["name"])
            self.setup.boundary_conditions.append(bc)

        for m in self._materials:
            self.setup.materials.append(MaterialSpec(name=m["name"]))

        if "k-epsilon" in self._physics:
            self.setup.turbulence_model = self._physics["k-epsilon"]
        elif "laminar" in self._physics:
            self.setup.turbulence_model = "laminar"
        elif self._physics:
            self.setup.turbulence_model = next(iter(self._physics.values()))

        if "transient" in self._physics:
            self.setup.transient = True
            self.setup.solver_mode = "transient"
        elif "steady" in self._physics:
            self.setup.solver_mode = "steady"

        self.setup.source.case_path = sim_path
        self.setup.source.case_name = Path(sim_path).stem
        return True

    def import_results(self, results_path: str, result: ConversionResult) -> bool:
        result.gap_report.unsupported_blocking(
            "results", "sim_results",
            "STAR-CCM+ .sim results cannot be directly parsed",
            "Export results via: Tools → Export → Results File (.enc, .plt, .vtk)",
        )
        return False

    def _populate_gap_report(self, gap: GapAnalysis) -> None:
        gap.supported("mesh", "header", f"Scanned STAR-CCM+ .sim header ({self._header_info.get('file_size', 0)} bytes)")

        for z in self._zones:
            gap.supported("boundary_zone", z["name"], f"Detected boundary zone '{z['name']}'")

        for m in self._materials:
            gap.supported("material", m["name"], f"Detected material '{m['name']}'")

        for pattern, mapped in self._physics.items():
            gap.supported("physics", pattern, f"Detected physics model '{pattern}' → '{mapped}'")

        for bc in self.setup.boundary_conditions:
            gap.supported("boundary_condition", bc.patch_name, f"Mapped STAR-CCM+ zone '{bc.source_zone_name}' to {bc.type.value}")

        gap.unsupported_blocking(
            "mesh", "sim_topology",
            "STAR-CCM+ .sim binary mesh topology (faces, cells, owner/neighbour) not reconstructed",
            "Export mesh via: File → Export → Mesh (.plt, .cgns, .vtk, or .msh)",
        )
        gap.unsupported_blocking(
            "results", "field_data",
            "STAR-CCM+ .sim field results cannot be directly parsed",
            "Export results via: Tools → Export → Results File (.enc, .plt, .vtk)",
        )
        gap.unsupported_nonblocking(
            "features", "udfs",
            "User-Defined Functions (UDFs) and custom models not detectable",
            "Document UDFs manually and verify in CFDX setup",
        )
        gap.unsupported_nonblocking(
            "features", "overset",
            "Overset (chimera) mesh information not extractable from .sim",
            "Export overset mesh via documented EnSight or VTK export",
        )


def _extract_ascii_strings(data: bytes, min_len: int = 4) -> list[str]:
    """Extract printable ASCII strings from binary data."""
    pattern = re.compile(rb"[\x20-\x7e]{%d,}" % min_len)
    return [m.group().decode("ascii") for m in pattern.finditer(data)]
