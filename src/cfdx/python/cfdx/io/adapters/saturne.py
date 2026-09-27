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
  - Performs SYNTACTIC EXTRACTION of Python case scripts (see the note below)
  - Recognises when a neutral mesh export (MED/CGNS/VTK) is present
  - Documents all limitations in the GapAnalysis

IMPORTANT — SYNTACTIC EXTRACTION vs VALIDATED SETUP MAPPING:
  `parse_python_setup()` is a *syntactic extractor*, not a validated setup
  mapper. It recognises one narrow syntactic shape:

      case.set("<literal key>", "<literal value>")

  and nothing else. A Code_Saturne case script is ordinary Python, so a real
  script can express the same setting in ways this adapter cannot resolve:

    - intermediate variables   ``ref = "k-epsilon"; case.set("model", ref)``
    - computed/concatenated     ``case.set("n", str(4 * mesh.n_cells))``
    - function calls            ``case.set("period", period(0.1, 2.0))``
    - dictionaries              ``params = {...}; case.set("d", params)``
    - loops                     ``for p in walls: case.set(p, "wall")``
    - imports / module state    ``from user_usr import *; cs_ = ...``
    - conditional logic         ``if REYNOLDS > 1e5: case.set(...)``

  Anything that is not a two-string-literal `case.set()` call is therefore
  *silently skipped*. Absence of an extracted entry means "not recognised",
  NOT "not set in the case". The parser deliberately does not attempt to
  evaluate the script; `self._python_extraction` records what was seen so the
  GapAnalysis can report partial, unvalidated coverage.

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


# Python constructs that a Code_Saturne case script may legitimately use to
# express a setting, but that `parse_python_setup()` CANNOT resolve. Their
# presence downgrades the extraction status from "literal_only" to "partial".
# These are counted for reporting only — no attempt is made to evaluate them.
_PYTHON_UNRESOLVABLE_PATTERNS: dict[str, re.Pattern] = {
    "import": re.compile(r"^\s*(?:from\s+\S+\s+import|import\s+\S+)", re.MULTILINE),
    "function_def": re.compile(r"^\s*def\s+\w+\s*\(", re.MULTILINE),
    "for_loop": re.compile(r"^\s*for\s+\w+.*:", re.MULTILINE),
    "while_loop": re.compile(r"^\s*while\s+.*:", re.MULTILINE),
    "conditional": re.compile(r"^\s*(?:if|elif)\s+.*:", re.MULTILINE),
    "dict_literal": re.compile(r"=\s*\{[^{}]*\}"),
    "function_call_value": re.compile(
        r"case\.set\(\s*[\"'][^\"']+[\"']\s*,\s*[\"'][^\"']+[\"']\s*,\s*[^\"')]"
    ),
    "variable_value": re.compile(
        r"case\.set\(\s*[\"'][^\"']+[\"']\s*,\s*(?![\"'])[A-Za-z_]"
    ),
    "f_string": re.compile(r"\bf[\"']"),
    "call_expression": re.compile(r"=\s*[A-Za-z_]\w*\([^\n]*\)"),
}

# Code_Saturne exports this adapter recognises.
_SATURNE_RESULT_EXTS: frozenset[str] = frozenset(
    {".med", ".cgns", ".vtk", ".vtu", ".case", ".geo", ".scl"}
)
# Subset of the above whose field data can ACTUALLY be read (via meshio).
# The rest are detect-only: presence is recorded, contents are not parsed.
_SATURNE_READABLE_RESULT_EXTS: frozenset[str] = frozenset({".vtk", ".vtu"})


class SaturneAdapter(SolverAdapter):
    solver_name = "Code_Saturne"
    format_name = "saturne_xml_py"

    def __init__(self, rules: Optional[MappingRules] = None) -> None:
        self.rules = rules or MappingRules.load_default()
        self.setup = CaseSetup(source=SourceInfo(solver="Code_Saturne", format="saturne_xml_py"))
        self._entries: list[dict] = []
        self._neutral_files: list[str] = []
        # Status of the Python setup extraction. "status" is always one of:
        #   "not_attempted" — no Python case file was present
        #   "partial"        — literal case.set() calls were extracted, but the
        #                      script contains constructs this adapter cannot
        #                      resolve (see _PYTHON_UNRESOLVABLE_PATTERNS)
        #   "literal_only"   — every case.set() call seen used two string
        #                      literals; still NOT a validated mapping
        self._python_extraction: dict = {"status": "not_attempted", "entries": 0, "unresolved": {}}

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

        # Mesh: only reported when a neutral export is available. Detection is
        # not parsing — no mesh geometry is read during convert().
        if self._neutral_files:
            result.gap_report.approximated(
                "mesh", "neutral_export",
                f"Detected {len(self._neutral_files)} neutral mesh export(s) "
                f"({', '.join(self._neutral_files)}); file presence only, mesh not read",
                "Call import_results() on a readable export (.vtk/.vtu) to load mesh and field data",
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
        """Syntactically extract literal ``case.set()`` calls from a case script.

        This is *syntactic extraction only*, not a validated setup mapping: the
        script is never executed and non-literal arguments are never evaluated.
        See the module docstring for the full statement of the limitation and
        ``_PYTHON_UNRESOLVABLE_PATTERNS`` for what gets skipped.
        """
        text = Path(py_file).read_text(errors="replace")

        # The ONLY shape understood: case.set('<literal key>', '<literal value>')
        set_pattern = re.compile(
            r'case\.set\(\s*[\'"]([^\'"]+)[\'"]\s*,\s*[\'"]([^\'"]+)[\'"]\s*\)'
        )
        matched_spans = []
        for match in set_pattern.finditer(text):
            key, value = match.group(1), match.group(2)
            self._entries.append({"key": key, "value": value, "section": "python", "attributes": {}})
            matched_spans.append(match.span())

        # Count case.set() calls whose arguments this regex could NOT capture
        # (variables, expressions, function calls, dicts, ...). These are
        # dropped, not approximated.
        all_set_calls = list(re.finditer(r'\bcase\.set\s*\(', text))
        captured = sum(1 for s in matched_spans)
        unresolved_calls = max(0, len(all_set_calls) - captured)

        # Record which unresolvable constructs are present in the script, so the
        # GapAnalysis can state concretely why coverage is partial.
        unresolved: dict[str, int] = {}
        for label, pattern in _PYTHON_UNRESOLVABLE_PATTERNS.items():
            n = len(pattern.findall(text))
            if n:
                unresolved[label] = n

        if unresolved or unresolved_calls:
            status = "partial"
        else:
            status = "literal_only"

        self._python_extraction = {
            "status": status,
            "entries": len(matched_spans),
            "unresolved": unresolved,
            "unresolved_case_set_calls": unresolved_calls,
            "file": str(py_file),
        }
        self.setup.source_metadata = {**self.setup.source_metadata, "python_extraction": self._python_extraction}

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
        """Import field data from a Code_Saturne neutral result export.

        Returns True ONLY when field data was actually read into
        ``result.scalar_fields`` / ``result.vec_fields``. Merely noticing that a
        neutral export exists is not a successful import, so the detect-only
        formats record a GapAnalysis finding and return False.
        """
        rp = Path(results_path)
        suffix = rp.suffix.lower()

        if suffix not in _SATURNE_RESULT_EXTS:
            result.gap_report.unsupported_nonblocking(
                "results", "saturne_results",
                f"Code_Saturne results file format '{suffix}' not supported: {results_path}",
                f"Use a neutral export format: {', '.join(sorted(_SATURNE_RESULT_EXTS))}",
            )
            return False

        if not rp.exists():
            result.gap_report.unsupported_nonblocking(
                "results", "saturne_results",
                f"Code_Saturne results file does not exist: {results_path}",
                "Verify the export path and re-export from Code_Saturne",
            )
            return False

        if suffix not in _SATURNE_READABLE_RESULT_EXTS:
            result.gap_report.unsupported_nonblocking(
                "results", f"neutral_{suffix[1:]}",
                f"Detected Code_Saturne results in neutral format ({suffix}) but no reader "
                f"is wired up: presence was recorded, no field data was read",
                f"Re-export as a readable format ({', '.join(sorted(_SATURNE_READABLE_RESULT_EXTS))}), "
                f"or convert the {suffix[1:].upper()} file with an external tool first",
            )
            return False

        try:
            import meshio
        except ImportError:
            result.gap_report.unsupported_nonblocking(
                "results", f"neutral_{suffix[1:]}",
                f"Code_Saturne results in {suffix} require the optional 'meshio' dependency, "
                "which is not installed; no field data was read",
                "Install meshio (pip install meshio) and re-run the import",
            )
            return False

        try:
            mesh = meshio.read(str(rp))
        except (Exception, SystemExit) as e:
            # meshio's internal readers call sys.exit(1) on an unreadable file,
            # so SystemExit must be caught here or it kills the whole CLI.
            reason = "meshio could not parse the file" if isinstance(e, SystemExit) else str(e)
            result.gap_report.unsupported_nonblocking(
                "results", f"neutral_{suffix[1:]}",
                f"Failed to read Code_Saturne results from {results_path}: {reason}",
                "Re-export the results from Code_Saturne, or use another conversion tool",
            )
            return False

        scalars: list[tuple[str, np.ndarray]] = []
        vecs: list[tuple[str, np.ndarray]] = []
        for name, data in (mesh.point_data or {}).items():
            arr = np.asarray(data)
            if arr.ndim == 2 and arr.shape[1] in (2, 3):
                vecs.append((name, arr))
            else:
                scalars.append((name, arr))

        if not scalars and not vecs:
            result.gap_report.unsupported_nonblocking(
                "results", f"neutral_{suffix[1:]}",
                f"Read {results_path} but found no point data fields; no field data imported",
                "Ensure the export includes solution fields (velocity, pressure, ...)",
            )
            return False

        result.scalar_fields.extend(scalars)
        result.vec_fields.extend(vecs)
        result.available_fields = {
            "scalar": [n for n, _ in scalars],
            "vector": [n for n, _ in vecs],
        }
        result.field_data_source = str(rp)
        if result.mesh is None and getattr(mesh, "points", None) is not None:
            result.mesh = {
                "points": np.asarray(mesh.points),
                "cells": {c.type: np.asarray(c.data) for c in mesh.cells},
            }

        result.gap_report.supported(
            "results", f"neutral_{suffix[1:]}",
            f"Read {len(scalars)} scalar and {len(vecs)} vector field(s) from {results_path}",
        )
        return True

    def _detect_neutral_exports(self, base: Path) -> None:
        neutral_exts = sorted(_SATURNE_RESULT_EXTS)
        if base.is_dir():
            for ext in neutral_exts:
                self._neutral_files.extend([str(f) for f in base.glob(f"*{ext}")])
        elif base.is_file():
            if base.suffix in neutral_exts:
                self._neutral_files.append(str(base))

    def _populate_gap_report(self, gap: GapAnalysis) -> None:
        xml_entries = [e for e in self._entries if e["section"] != "python"]

        if xml_entries:
            gap.supported(
                "config", "saturne_xml",
                f"Code_Saturne XML setup parsed ({len(xml_entries)} recognised element(s)/attribute(s))",
            )
        else:
            gap.unsupported_nonblocking(
                "config", "saturne_xml",
                "No Code_Saturne XML setup entries were recognised",
                "Provide the XML case file, or map physics/BCs manually in the CFDX setup",
            )

        self._report_python_extraction(gap)

        for entry in self._entries:
            gap.supported("config", entry["key"], f"Parsed '{entry['key']}' = '{entry['value']}' from setup")

        for bc in self.setup.boundary_conditions:
            gap.supported("boundary_condition", bc.patch_name, f"Mapped Code_Saturne BC '{bc.source_zone_name}' to {bc.type.value}")

        for mat in self.setup.materials:
            gap.supported("material", mat.name, f"Parsed material '{mat.name}' (ρ={mat.density})")

        if self._neutral_files:
            for nf in self._neutral_files:
                gap.approximated(
                    "mesh", "neutral_export",
                    f"Detected neutral export: {nf} (presence only, contents not read)",
                )
            gap.unsupported_nonblocking(
                "results", "neutral_field_data",
                "Neutral exports were detected but no field data was read; "
                f"readable formats are {', '.join(sorted(_SATURNE_READABLE_RESULT_EXTS))}",
                "Call import_results() with a .vtk/.vtu export to populate field data",
            )
        else:
            gap.unsupported_blocking(
                "mesh", "neutral_mesh",
                "No neutral mesh export (MED/CGNS/VTK/EnSight) found",
                "Export mesh via: Computation → Export → Select format (MED, CGNS, VTK, EnSight)",
            )

    def _report_python_extraction(self, gap: GapAnalysis) -> None:
        """Report the Python setup extraction as partial, never as full support."""
        ex = self._python_extraction
        status = ex["status"]

        if status == "not_attempted":
            gap.unsupported_nonblocking(
                "config", "saturne_py",
                "No Code_Saturne Python case script was found or parsed",
                "Provide the case .py file, or set up boundaries/physics manually",
            )
            return

        detail = (
            f"Syntactic extraction only: {ex['entries']} literal case.set() call(s) "
            f"recognised in {ex.get('file', 'the case script')}. This is NOT a validated "
            f"setup mapping — the script is never executed."
        )

        if status == "literal_only":
            gap.approximated(
                "config", "saturne_py",
                detail + " No unresolvable constructs were detected, but the mapping remains unvalidated.",
            )
            return

        reasons = ", ".join(f"{k}×{v}" for k, v in sorted(ex["unresolved"].items()))
        missed = ex.get("unresolved_case_set_calls", 0)
        gap.approximated(
            "config", "saturne_py",
            detail + f" Unresolvable constructs found ({reasons}); "
            f"{missed} case.set() call(s) with non-literal arguments were skipped entirely.",
        )
        gap.unsupported_nonblocking(
            "features", "python_scripting",
            "Code_Saturne Python scripting is only syntactically scanned. Settings expressed as "
            "intermediate variables, expressions, function calls, dictionaries, loops, imports or "
            "computed values are not extracted and cannot be distinguished from unset settings.",
            "Review the case script manually and confirm every boundary and physics option in the CFDX setup",
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
