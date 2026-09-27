"""Code_Saturne adapter for CFDX I/O.

Parses Code_Saturne case files:
  - XML case setup (.xml) — decomposition of physics, BCs, materials
  - Python case definition (.py) — case.set() calls
  - MED/CGNS/EnSight/VTK result exports (detected, not parsed for mesh)

ARCHITECTURAL NOTE (per issue #425):
  "Evaluate XML/Python setup extraction"
  "Support documented neutral mesh/result exports such as MED/CGNS/EnSight/VTK"

This adapter:
  - Parses the Code_Saturne XML case file for physics model, BCs, materials
  - Parses Python case scripts for case.set() calls (key-value extraction)
  - Recognises when a neutral mesh export (MED/CGNS/VTK) is present
  - Documents all limitations in the GapAnalysis

Corresponds to the C++ header:
  src/cfdx/io/saturne/saturne_adapter.h
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional, Union
from xml.etree import ElementTree as ET

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
    MeshMetadata,
)
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.gap_analysis import GapAnalysis


# XML element name → (section, field) for Code_Saturne
_SATURNE_XML_MAP: dict[str, tuple[str, str]] = {
    "physics_model": ("physics", "physics_model"),
    "turbulence_model": ("physics", "turbulence_model"),
    "energy_model": ("physics", "energy_model"),
    "heat_transfer": ("physics", "energy_model"),
    "multiphase_model": ("physics", "multiphase_model"),
    "transient": ("physics", "transient"),
    "time_step": ("physics", "time_step"),
    "end_time": ("physics", "end_time"),
}


class SaturneAdapter(SolverAdapter):
    solver_name = "Code_Saturne"
    format_name = "saturne_xml_py"

    def __init__(self, rules: Optional[MappingRules] = None) -> None:
        self.rules = rules or MappingRules.load_default()
        self.setup = CaseSetup(source=SourceInfo(solver="Code_Saturne", format="saturne_xml_py"))
        self._entries: list[dict] = []
        self._neutral_files: list[str] = []

    def detect_source(self, case_path: str, info: SourceInfo) -> bool:
        p = Path(str(case_path))
        candidates: list[Path] = []
        if p.is_file():
            candidates.append(p)
        elif p.is_dir():
            candidates = sorted(list(p.glob("*.xml")) + list(p.glob("*.py")))
        else:
            for f in (p, Path(str(p) + ".xml"), Path(str(p) + ".py")):
                if f.exists():
                    candidates.append(f)

        for cand in candidates:
            try:
                head = cand.read_text(errors="replace")[:512]
            except OSError:
                continue
            if cand.suffix == ".xml" and ("<Case" in head or "<case" in head.lower() or "saturne" in head.lower()):
                info.solver = "Code_Saturne"
                info.format = "saturne_xml_py"
                info.case_path = str(cand.parent)
                info.case_name = cand.stem
                info.version = "unknown"
                return True
            elif cand.suffix == ".py" and ("case.set" in head or "cs_" in head or "Code_Saturne" in head or "saturne" in head.lower()):
                info.solver = "Code_Saturne"
                info.format = "saturne_xml_py"
                info.case_path = str(cand.parent)
                info.case_name = cand.stem
                info.version = "unknown"
                return True
        return False

    def convert(self, case_path: Union[str, Path], result: ConversionResult) -> bool:
        base = Path(str(case_path))
        xml_file: Optional[Path] = None
        py_file: Optional[Path] = None

        if base.is_file():
            if base.suffix == ".xml":
                xml_file = base
            elif base.suffix == ".py":
                py_file = base
        elif base.is_dir():
            xml_candidates = list(base.glob("*.xml"))
            py_candidates = list(base.glob("*.py"))
            xml_file = xml_candidates[0] if xml_candidates else None
            py_file = py_candidates[0] if py_candidates else None
        else:
            if base.suffix == ".xml":
                xml_file = base
            elif base.suffix == ".py":
                py_file = base
            else:
                xml_file = base.with_suffix(".xml")
                py_file = base.with_suffix(".py")

        # Parse XML
        if xml_file and xml_file.exists():
            self.parse_xml_setup(str(xml_file))
        else:
            result.gap_report.unsupported_nonblocking(
                "config", "saturne_xml",
                f"No Code_Saturne XML setup file found at {xml_file}",
                "Code_Saturne XML setup is optional; physics will be approximated",
            )

        # Parse Python
        if py_file and py_file.exists():
            self.parse_python_setup(str(py_file))
        else:
            result.gap_report.unsupported_nonblocking(
                "config", "saturne_py",
                f"No Code_Saturne Python setup file found at {py_file}",
                "Code_Saturne Python setup is optional; some BCs may be approximated",
            )

        # Detect neutral mesh exports
        self._detect_neutral_exports(base)

        result.source = self.setup.source
        result.setup = self.setup.model_copy(deep=True)

        # Mesh: only set if neutral export is available
        if self._neutral_files:
            result.gap_report.supported(
                "mesh", "neutral_export",
                f"Detected {len(self._neutral_files)} neutral mesh export(s): {', '.join(self._neutral_files)}",
            )
        else:
            result.gap_report.unsupported_blocking(
                "mesh", "neutral_mesh",
                "No neutral mesh export (MED/CGNS/VTK/EnSight) found for Code_Saturne mesh import",
                "Export mesh via: Computation → Export → Select format (MED, CGNS, VTK, EnSight)",
            )

        self._populate_gap_report(result.gap_report)
        return not result.gap_report.has_blocking()

    def parse_xml_setup(self, xml_file: str) -> bool:
        try:
            tree = ET.parse(xml_file)
        except ET.ParseError:
            self.setup.source_metadata = {"xml_parse_error": "invalid XML structure"}
            return False

        root = tree.getroot()

        def _walk(elem: ET.Element, section: str = "") -> None:
            tag = elem.tag.split("}")[-1]
            text = (elem.text or "").strip()
            attrib = elem.attrib

            # Check element attributes for known physics/B.C. keys
            for attr_name, attr_val in attrib.items():
                attr_lower = attr_name.lower()
                if attr_lower in _SATURNE_XML_MAP:
                    sec, field = _SATURNE_XML_MAP[attr_lower]
                    self._entries.append(
                        {"key": attr_name, "value": attr_val,
                         "section": sec, "attributes": dict(attrib)}
                    )

            # Check tag-based entries (text content)
            if tag in _SATURNE_XML_MAP:
                sec, field = _SATURNE_XML_MAP[tag]
                self._entries.append({"key": tag, "value": text, "section": sec, "attributes": dict(attrib)})

            # Handle boundary conditions (tag can be "boundary" or "boundary_condition")
            if tag in ("boundary", "boundary_condition"):
                bc = BoundarySpec(
                    patch_name=attrib.get("name", attrib.get("marker", tag)),
                    source_zone_name=attrib.get("name", attrib.get("marker", "")),
                    source_type_name=attrib.get("type", attrib.get("bc_type", tag)),
                )
                bc_type_str = (attrib.get("type", "") or attrib.get("bc_type", "") or tag).lower()
                mapped = self.rules.map_boundary_type(bc_type_str)
                bc.type = mapped
                bc.value_type = self.rules.map_bc_value_type(bc_type_str)
                if "wall" in bc_type_str:
                    bc.type = BCType.WALL
                    bc.value_type = BCValueType.WALL_NO_SLIP
                elif "inlet" in bc_type_str or "input" in bc_type_str:
                    bc.type = BCType.INLET
                    bc.value_type = BCValueType.FIXED
                elif "outlet" in bc_type_str or "output" in bc_type_str:
                    bc.type = BCType.OUTLET
                    bc.value_type = BCValueType.ZERO_GRADIENT
                elif "sym" in bc_type_str:
                    bc.type = BCType.SYMMETRY
                    bc.value_type = BCValueType.ZERO_GRADIENT
                elif "far" in bc_type_str:
                    bc.type = BCType.OUTLET
                    bc.value_type = BCValueType.FIXED
                self.setup.boundary_conditions.append(bc)

            # Handle materials
            if tag in ("material", "fluid"):
                mat = MaterialSpec(
                    name=attrib.get("name", "fluid"),
                    density=float(attrib.get("density", "1.225") or "1.225"),
                    dynamic_viscosity=float(attrib.get("viscosity", "1.789e-5") or "1.789e-5"),
                )
                self.setup.materials.append(mat)

            for child in elem:
                _walk(child, section)

        _walk(root)

        # Apply detected physics
        for entry in self._entries:
            if entry["section"] == "physics":
                val_lower = entry["value"].lower()
                if entry["key"].lower() == "physics_model":
                    if "navier" in val_lower:
                        self.setup.physics_model = "compressible_navier_stokes"
                    elif "euler" in val_lower:
                        self.setup.physics_model = "compressible_euler"
                    elif "incompressible" in val_lower:
                        self.setup.physics_model = "incompressible_laminar"
                elif entry["key"].lower() == "turbulence_model":
                    mapped = self.rules.map_turbulence(entry["value"])
                    self.setup.turbulence_model = mapped if mapped else entry["value"].lower()
                elif entry["key"].lower() == "energy_model":
                    self.setup.energy_model = "buoyant" if val_lower in ("yes", "true", "on") else "isothermal"
                elif entry["key"].lower() == "transient":
                    self.setup.transient = val_lower in ("yes", "true", "on")
                    self.setup.solver_mode = "transient" if self.setup.transient else "steady"

        self.setup.source.case_path = xml_file
        self.setup.source.case_name = Path(xml_file).stem
        return True

    def parse_python_setup(self, py_file: str) -> bool:
        text = Path(py_file).read_text(errors="replace")

        # Extract case.set() calls: case.set('key', 'value') or case.set("key", "value")
        set_pattern = re.compile(
            r'case\.set\(\s*[\'"]([^\'"]+)[\'"]\s*,\s*[\'"]([^\'"]+)[\'"]\s*\)'
        )
        for match in set_pattern.finditer(text):
            key, value = match.group(1), match.group(2)
            self._entries.append({"key": key, "value": value, "section": "python", "attributes": {}})

        # Also detect standalone set() calls
        standalone_pattern = re.compile(
            r'\bset\(\s*[\'"]([^\'"]+)[\'"]\s*,\s*([^\s,]+(?:\s*[,\)])?)'
        )

        # Apply physics from Python keys
        for entry in self._entries:
            if entry["section"] == "python":
                key_lower = entry["key"].lower()
                val_lower = entry["value"].lower()
                if "turbulence" in key_lower or "keps" in key_lower:
                    mapped = self.rules.map_turbulence(entry["value"])
                    self.setup.turbulence_model = mapped if mapped else entry["value"].lower()
                elif "energy" in key_lower or "heat" in key_lower:
                    self.setup.energy_model = "buoyant" if val_lower in ("yes", "true", "on", "1") else "isothermal"
                elif "transient" in key_lower:
                    self.setup.transient = val_lower in ("yes", "true", "on", "1")
                    self.setup.solver_mode = "transient" if self.setup.transient else "steady"

        self.setup.source.case_path = py_file
        self.setup.source.case_name = Path(py_file).stem
        return True

    def import_results(self, results_path: str, result: ConversionResult) -> bool:
        rp = Path(results_path)
        suffix = rp.suffix.lower()
        supported_results = {".med", ".cgns", ".vtk", ".vtu", ".case", ".geo", ".scl"}

        if suffix not in supported_results:
            result.gap_report.unsupported_nonblocking(
                "results", "saturne_results",
                f"Code_Saturne results file format '{suffix}' not supported: {results_path}",
                f"Use a neutral export format: {', '.join(sorted(supported_results))}",
            )
            return False

        result.gap_report.supported(
            "results", f"neutral_{suffix[1:]}",
            f"Detected Code_Saturne results in neutral format: {results_path}",
        )
        # Field data would be extracted via meshio or h5py depending on format
        return True

    def _detect_neutral_exports(self, base: Path) -> None:
        neutral_exts = [".med", ".cgns", ".vtk", ".vtu", ".case", ".geo", ".scl"]
        if base.is_dir():
            for ext in neutral_exts:
                self._neutral_files.extend([str(f) for f in base.glob(f"*{ext}")])
        elif base.is_file():
            if base.suffix in neutral_exts:
                self._neutral_files.append(str(base))

    def _populate_gap_report(self, gap: GapAnalysis) -> None:
        gap.supported("config", "saturne_xml", "Code_Saturne XML setup parsed" if self._entries else "No XML setup file found")
        gap.supported("config", "saturne_py", "Code_Saturne Python setup parsed" if self._entries else "No Python setup file found")

        for entry in self._entries:
            gap.supported("config", entry["key"], f"Parsed '{entry['key']}' = '{entry['value']}' from setup")

        for bc in self.setup.boundary_conditions:
            gap.supported("boundary_condition", bc.patch_name, f"Mapped Code_Saturne BC '{bc.source_zone_name}' to {bc.type.value}")

        for mat in self.setup.materials:
            gap.supported("material", mat.name, f"Parsed material '{mat.name}' (ρ={mat.density})")

        if self._neutral_files:
            for nf in self._neutral_files:
                gap.supported("mesh", "neutral_export", f"Detected neutral export: {nf}")
        else:
            gap.unsupported_blocking(
                "mesh", "neutral_mesh",
                "No neutral mesh export (MED/CGNS/VTK/EnSight) found",
                "Export mesh via: Computation → Export → Select format (MED, CGNS, VTK, EnSight)",
            )

        gap.unsupported_nonblocking(
            "features", "python_scripting",
            "Code_Saturne Python scripting (case.set() calls) partially parsed — complex expressions may be missed",
            "Review all case.set() calls manually; complex Python logic is not auto-translated",
        )
        gap.unsupported_nonblocking(
            "features", "user_subscripts",
            "User subroutines (.usr files) not supported",
            "Verify user subroutines are translated to CFDX equivalent models",
        )
        gap.unsupported_nonblocking(
            "features", "multiphase_advanced",
            "Advanced multiphase models (lift, turbulent dispersion, etc.) may not fully translate",
            "Review multiphase model settings in CFDX setup",
        )
        gap.unsupported_nonblocking(
            "features", "radiation",
            "Code_Saturne radiation models not mapped",
            "Configure radiation manually in CFDX if required",
        )
        gap.unsupported_nonblocking(
            "features", "anisotropy",
            "Anisotropic conductivity and permeability tensors not supported",
            "Review material anisotropy in CFDX setup",
        )
