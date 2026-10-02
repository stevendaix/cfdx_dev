"""Unified CFDX conversion API.

This module is the single entry point for turning an external solver case into
the CFDX canonical representation.  It sits on top of the low-level solver
adapters in :mod:`cfdx.io.adapters` and normalises whatever those adapters
produce into one stable :class:`ConversionResult` shape:

    from cfdx.io.converter import get_converter, convert, ConversionResult

    converter = get_converter("fluent")
    result = converter.convert(source="case.cas.h5", output="case.cfdx.h5")

    result.case             # CaseSetup
    result.mesh             # mesh dict or None
    result.fields           # list of (name, array)
    result.boundaries       # list[BoundarySpec]
    result.materials        # list[MaterialSpec]
    result.solver_settings  # dict
    result.gap_analysis     # GapAnalysis
    result.success          # bool

Design rules:

* Adapters are never required to know about the canonical layout — this module
  adapts them, and wraps their exceptions into blocking gap findings.
* Nothing is ever silently dropped.  Whatever the source format does not carry
  is recorded in :class:`~cfdx.io.gap_analysis.GapAnalysis`.
* Adapters are stateful, therefore a fresh converter instance is created for
  every :func:`get_converter` call.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Optional, Union

import numpy as np

from cfdx.io.gap_analysis import GapAnalysis, GapAnalysisReport
from cfdx.io.hdf5_writer import write_case_cfdx_h5
from cfdx.io.numerical_selection import build_numerical_selection
from cfdx.io.interfaces import ConversionResult as AdapterResult
from cfdx.io.schema import (
    CFDX_SCHEMA_VERSION,
    BoundarySpec,
    CaseSetup,
    MaterialSpec,
    SourceInfo,
)


# --------------------------------------------------------------------------
# Canonical result
# --------------------------------------------------------------------------


@dataclass
class ConversionResult:
    """Canonical CFDX representation of a converted solver case."""

    source: SourceInfo = field(default_factory=SourceInfo)
    case: Optional[CaseSetup] = None
    mesh: Optional[dict] = None
    fields: list[tuple[str, np.ndarray]] = field(default_factory=list)
    boundaries: list[BoundarySpec] = field(default_factory=list)
    materials: list[MaterialSpec] = field(default_factory=list)
    solver_settings: dict[str, Any] = field(default_factory=dict)
    gap_analysis: GapAnalysis = field(default_factory=GapAnalysis)
    output_path: Optional[str] = None

    # -- convenience -------------------------------------------------------

    @property
    def success(self) -> bool:
        """True when no blocking incompatibility was recorded."""
        return not self.gap_analysis.has_blocking()

    @property
    def has_mesh(self) -> bool:
        return self.mesh is not None and self.mesh.get("points") is not None

    @property
    def field_names(self) -> list[str]:
        return [name for name, _ in self.fields]

    @property
    def boundary_names(self) -> list[str]:
        return [bc.patch_name for bc in self.boundaries]

    @property
    def material_names(self) -> list[str]:
        return [mat.name for mat in self.materials]

    def find_field(self, name: str) -> Optional[np.ndarray]:
        for field_name, data in self.fields:
            if field_name == name:
                return data
        return None

    def to_adapter_result(self) -> AdapterResult:
        """Convert back to the low-level adapter result shape.

        Used by the HDF5 writer, which consumes the ``cfdx.io.interfaces``
        flavour of the result.
        """
        legacy = AdapterResult(source=self.source, setup=self.case)
        legacy.mesh = self.mesh
        legacy.gap_report = self.gap_analysis
        legacy.scalar_fields = []
        legacy.vec_fields = []
        for name, data in self.fields:
            arr = np.asarray(data)
            if arr.ndim >= 2:
                legacy.vec_fields.append((name, arr))
            else:
                legacy.scalar_fields.append((name, arr))
        return legacy

    def summary(self) -> dict:
        """JSON-serialisable summary of the conversion."""
        points = None
        n_cells = 0
        cell_types: list[str] = []
        if self.mesh is not None and self.mesh.get("points") is not None:
            pts = np.asarray(self.mesh["points"])
            points = int(pts.shape[0])
            cells = self.mesh.get("cells") or {}
            cell_types = list(cells.keys())
            n_cells = sum(len(v) for v in cells.values()) if cells else 0
            if not n_cells:
                n_cells = int(self.mesh.get("n_cells", 0) or 0)

        gap = self.gap_analysis
        return {
            "solver": self.source.solver,
            "version": self.source.version,
            "format": self.source.format,
            "source": self.source.case_path,
            "output": self.output_path,
            "success": self.success,
            "mesh": {
                "points": points,
                "cells": n_cells,
                "cell_types": cell_types,
            },
            "fields": self.field_names,
            "boundaries": self.boundary_names,
            "materials": self.material_names,
            "gap_analysis": {
                "supported": gap.n_supported(),
                "approximated": gap.n_approximated(),
                "unsupported_nonblocking": gap.n_unsupported_nonblocking(),
                "unsupported_blocking": gap.n_unsupported_blocking(),
                "unavailable": gap.n_unavailable(),
                "total": len(gap.findings),
            },
        }


# --------------------------------------------------------------------------
# Normalisation helpers
# --------------------------------------------------------------------------

_SOLVER_SETTING_KEYS = (
    "physics_model",
    "turbulence_model",
    "energy_model",
    "multiphase_model",
    "radiation_model",
    "transient",
    "time_step",
    "end_time",
    "max_time_steps",
    "solver_mode",
    "units",
    "gravity_vector",
    "ref_length",
    "ref_density",
    "ref_velocity",
    "ref_temperature",
    "ref_pressure",
)


def _copy_gap(src: GapAnalysis, dst: GapAnalysis) -> None:
    for f in src.findings:
        dst.add(f.severity, f.category, f.feature, f.detail, f.suggestion)


def _solver_settings(
    setup: Optional[CaseSetup], legacy: AdapterResult
) -> dict[str, Any]:
    settings: dict[str, Any] = {}
    if setup is not None:
        for key in _SOLVER_SETTING_KEYS:
            settings[key] = getattr(setup, key)
        settings["numerics"] = setup.numerics.model_dump(mode="json")
        settings["initial_condition"] = setup.initial_condition.model_dump(mode="json")
        settings["mesh_info"] = setup.mesh_info.model_dump(mode="json")
    if legacy.available_fields:
        settings["available_fields"] = dict(legacy.available_fields)
    if legacy.field_data_source:
        settings["field_data_source"] = legacy.field_data_source
    return settings


def _record_missing_data(result: ConversionResult) -> None:
    """Record absent data in the gap analysis.

    A category already mentioned by the adapter is left alone so we do not
    duplicate findings; anything the adapter said nothing about is still
    recorded so nothing disappears silently.
    """
    gap = result.gap_analysis
    reported = {f.category.split("_")[0].split(".")[0].lower() for f in gap.findings}

    checks = (
        (
            "mesh",
            "mesh_geometry",
            result.has_mesh,
            "Source contains no readable mesh geometry",
        ),
        (
            "fields",
            "field_data",
            bool(result.fields),
            "Source contains no field/solution data",
        ),
        (
            "material",
            "material_properties",
            bool(result.materials),
            "Source contains no material definitions; CFDX air defaults are implied",
        ),
        (
            "boundary",
            "boundary_conditions",
            bool(result.boundaries),
            "Source contains no boundary condition definitions",
        ),
    )

    for category, feature, present, detail in checks:
        if present or category in reported:
            continue
        gap.unavailable(category, feature, detail)


def _normalize(
    legacy: AdapterResult,
    path: Path,
    converter_name: str,
) -> ConversionResult:
    """Normalise an adapter result into the canonical :class:`ConversionResult`."""
    setup = legacy.setup
    source = legacy.source
    if not source.solver:
        source.solver = converter_name
    if not source.case_path:
        source.case_path = str(path)

    result = ConversionResult(source=source, case=setup, mesh=legacy.mesh)
    if setup is not None:
        result.boundaries = list(setup.boundary_conditions)
        result.materials = list(setup.materials)
        if not source.solver:
            source.solver = setup.source.solver
        if not source.format:
            source.format = setup.source.format
    for name, data in list(legacy.scalar_fields) + list(legacy.vec_fields):
        result.fields.append((name, data))
    result.solver_settings = _solver_settings(setup, legacy)

    if setup is not None:
        setup.numerics.selection = build_numerical_selection(setup.numerics)
        result.solver_settings["numerics"]["selection"] = setup.numerics.selection.model_dump(mode="json")

    _copy_gap(legacy.gap_report, result.gap_analysis)
    _record_missing_data(result)
    return result


# --------------------------------------------------------------------------
# HDF5 output
# --------------------------------------------------------------------------


def _gap_to_json(gap: GapAnalysis) -> str:
    return json.dumps(
        {
            "has_blocking": gap.has_blocking(),
            "summary": {
                "supported": gap.n_supported(),
                "approximated": gap.n_approximated(),
                "unsupported_nonblocking": gap.n_unsupported_nonblocking(),
                "unsupported_blocking": gap.n_unsupported_blocking(),
                "unavailable": gap.n_unavailable(),
                "total": len(gap.findings),
            },
            "findings": [
                {
                    "severity": f.severity.value,
                    "category": f.category,
                    "feature": f.feature,
                    "detail": f.detail,
                    "suggestion": f.suggestion,
                }
                for f in gap.findings
            ],
        }
    )


def _write_csr_h5(result: ConversionResult, output: Path, csr: dict) -> None:
    """Write a result whose mesh is already in CFDX CSR form."""
    import h5py

    with h5py.File(output, "w") as f:
        f.attrs["format_version"] = "1"
        f.attrs["schema_version"] = str(CFDX_SCHEMA_VERSION)
        f.attrs["cfdx_version"] = "0.7"
        f.attrs["source_solver"] = result.source.solver or "unknown"
        f.attrs["source_format"] = result.source.format or "unknown"
        f.attrs["source_version"] = result.source.version or "unknown"
        f.attrs["source_case_path"] = result.source.case_path or ""
        f.attrs["source_case_name"] = result.source.case_name or ""
        if result.case is not None:
            f.attrs["case_setup_json"] = json.dumps(
                result.case.model_dump(mode="json"), default=str
            )
        f.attrs["gap_report_json"] = _gap_to_json(result.gap_analysis)

        for name in (
            "points",
            "face_vertices",
            "face_offsets",
            "owner",
            "neighbour",
            "cell_faces",
            "cell_offsets",
        ):
            if name in csr and csr[name] is not None:
                f.create_dataset(name, data=np.asarray(csr[name]))
        f.attrs["mesh_topology"] = "cfdx-csr-v1"

        fields_grp = f.create_group("fields")
        if result.fields:
            sgrp = fields_grp.create_group("scalar")
            vgrp = fields_grp.create_group("vector")
            for name, data in result.fields:
                arr = np.asarray(data)
                if arr.size == 0:
                    continue
                if arr.ndim >= 2:
                    vgrp.create_dataset(name, data=arr.astype(np.float64))
                else:
                    sgrp.create_dataset(name, data=arr.astype(np.float64))


def write_result(
    result: ConversionResult,
    output: Union[str, Path],
    write_reports: bool = True,
) -> Path:
    """Write a canonical result to ``output`` (HDF5 by default)."""
    output = Path(output)
    if output.parent and not output.parent.exists():
        output.parent.mkdir(parents=True, exist_ok=True)

    csr = (result.mesh or {}).get("cfdx_csr")
    if csr:
        _write_csr_h5(result, output, csr)
    else:
        write_case_cfdx_h5(result.to_adapter_result(), output)

    if write_reports:
        report = GapAnalysisReport(result.gap_analysis)
        stem = output.stem
        (output.parent / f"{stem}_gap_analysis.md").write_text(report.to_markdown())
        (output.parent / f"{stem}_gap_analysis.json").write_text(report.to_json())

    result.output_path = str(output)
    return output


# --------------------------------------------------------------------------
# Converters
# --------------------------------------------------------------------------


class SolverConverter:
    """Base class for all unified CFDX converters."""

    solver_name: str = "unknown"
    format_name: str = "unknown"
    extensions: tuple[str, ...] = ()
    accepts_directory: bool = False

    def convert(
        self,
        source: Union[str, Path],
        output: Optional[Union[str, Path]] = None,
        write_reports: bool = True,
    ) -> ConversionResult:
        """Convert ``source`` into the CFDX canonical representation.

        Args:
            source: Solver case file or case directory.
            output: Optional output file (``*.cfdx.h5`` by convention).
            write_reports: Write the Markdown/JSON gap analysis next to the output.
        """
        result = self.read(source)
        if output is not None:
            write_result(result, output, write_reports=write_reports)
        return result

    def read(self, source: Union[str, Path]) -> ConversionResult:
        path = Path(str(source))
        if not path.exists():
            result = ConversionResult(source=SourceInfo(solver=self.solver_name))
            result.gap_analysis.unsupported_blocking(
                "source", "source_path", f"Source path does not exist: {path}"
            )
            _record_missing_data(result)
            return result
        if path.is_dir() and not self.accepts_directory:
            result = ConversionResult(source=SourceInfo(solver=self.solver_name))
            result.gap_analysis.unsupported_blocking(
                "source",
                "source_path",
                f"{self.solver_name} conversion expects a file, got directory: {path}",
            )
            _record_missing_data(result)
            return result
        return self._read(path)

    def _read(self, path: Path) -> ConversionResult:
        raise NotImplementedError

    def __repr__(self) -> str:  # pragma: no cover - debugging aid
        return f"<{type(self).__name__} solver={self.solver_name!r}>"


class AdapterConverter(SolverConverter):
    """Wraps a :class:`cfdx.io.interfaces.SolverAdapter` implementation."""

    adapter_class: type

    def _make_adapter(self):
        return self.adapter_class()

    def _read(self, path: Path) -> ConversionResult:
        adapter = self._make_adapter()
        legacy = AdapterResult(source=SourceInfo())
        try:
            adapter.detect_source(str(path), legacy.source)
        except Exception as exc:  # detection is best-effort
            legacy.gap_report.unsupported_nonblocking(
                "source",
                "source_detection",
                f"Source detection failed: {type(exc).__name__}: {exc}",
                "Provide the solver name explicitly via get_converter(<solver>)",
            )

        try:
            adapter.convert(str(path), legacy)
        except Exception as exc:
            legacy.gap_report.unsupported_blocking(
                "adapter",
                "exception",
                f"{self.solver_name} adapter raised {type(exc).__name__}: {exc}",
                "Report this as a bug with the source case attached",
            )

        if legacy.source is None:
            legacy.source = SourceInfo()
        return _normalize(legacy, path, self.solver_name)


class FluentConverter(AdapterConverter):
    solver_name = "Fluent"
    format_name = "cas_dat"
    extensions = (".cas", ".dat", ".cas.h5", ".dat.h5", ".cas_h5")

    def _make_adapter(self):
        from cfdx.io.adapters.fluent import FluentAdapter

        return FluentAdapter()


class Su2Converter(AdapterConverter):
    solver_name = "SU2"
    format_name = "su2"
    extensions = (".su2", ".cfg")
    accepts_directory = True

    def _make_adapter(self):
        from cfdx.io.adapters.su2 import Su2Adapter

        return Su2Adapter()


class StarCcmConverter(AdapterConverter):
    solver_name = "STAR-CCM+"
    format_name = "sim"
    extensions = (".sim",)

    def _make_adapter(self):
        from cfdx.io.adapters.starccm import StarCCMAdapter

        return StarCCMAdapter()


class SaturneConverter(AdapterConverter):
    solver_name = "Code_Saturne"
    format_name = "saturne_xml_py"
    extensions = (".xml", ".py", ".med", ".cgns")
    accepts_directory = True

    def _make_adapter(self):
        from cfdx.io.adapters.saturne import SaturneAdapter

        return SaturneAdapter()


class OpenFoamConverter(AdapterConverter):
    solver_name = "OpenFOAM"
    format_name = "openfoam"
    extensions = ()
    accepts_directory = True

    def _make_adapter(self):
        from cfdx.io.adapters.openfoam import OpenFOAMAdapter

        return OpenFOAMAdapter()


class CfdxConverter(SolverConverter):
    """Reads/normalises an existing ``*.cfdx.h5`` file."""

    solver_name = "CFDX"
    format_name = "cfdx_h5"
    extensions = (".cfdx.h5",)

    def _read(self, path: Path) -> ConversionResult:
        import h5py

        result = ConversionResult(source=SourceInfo(solver="CFDX", format="cfdx_h5"))
        try:
            f = h5py.File(path, "r")
        except OSError as exc:
            result.gap_analysis.unsupported_blocking(
                "source",
                "hdf5_read",
                f"Cannot open CFDX HDF5 file {path}: {exc}",
            )
            _record_missing_data(result)
            return result

        with f:
            attrs = f.attrs
            result.source.solver = str(attrs.get("source_solver", "CFDX") or "CFDX")
            result.source.format = str(attrs.get("source_format", "cfdx_h5") or "cfdx_h5")
            result.source.version = str(attrs.get("source_version", "") or "")
            result.source.case_path = str(attrs.get("source_case_path", "") or "")
            result.source.case_name = str(attrs.get("source_case_name", "") or "")
            if not result.source.case_path:
                result.source.case_path = str(path)

            setup_json = attrs.get("case_setup_json")
            if isinstance(setup_json, bytes):
                setup_json = setup_json.decode("utf-8", "replace")
            if setup_json:
                try:
                    result.case = CaseSetup.model_validate_json(str(setup_json))
                except Exception as exc:
                    result.gap_analysis.unsupported_nonblocking(
                        "schema", "case_setup_json",
                        f"Cannot parse stored case setup: {exc}",
                        "Re-convert from the original solver case",
                    )
            if result.case is not None:
                result.boundaries = list(result.case.boundary_conditions)
                result.materials = list(result.case.materials)
                result.solver_settings = _solver_settings(result.case, AdapterResult(source=result.source))

            gap_json = attrs.get("gap_report_json")
            if isinstance(gap_json, bytes):
                gap_json = gap_json.decode("utf-8", "replace")
            if gap_json:
                _load_gap_json(str(gap_json), result.gap_analysis)

            csr = {
                name: f[name][()]
                for name in (
                    "points",
                    "face_vertices",
                    "face_offsets",
                    "owner",
                    "neighbour",
                    "cell_faces",
                    "cell_offsets",
                )
                if name in f
            }
            if "points" in csr:
                points = np.asarray(csr["points"]).reshape(-1, 3)
                result.mesh = {
                    "points": points,
                    "cells": {},
                    "cfdx_csr": csr,
                    "n_cells": int(len(csr["cell_offsets"]) - 1) if "cell_offsets" in csr else 0,
                }

            for kind in ("scalar", "vector"):
                grp = f.get(f"fields/{kind}")
                if grp is None:
                    continue
                for name in grp:
                    result.fields.append((name, np.asarray(grp[name][()])))

        result.gap_analysis.supported(
            "schema", "cfdx_h5", f"Normalised CFDX HDF5 container: {path.name}"
        )
        _record_missing_data(result)
        return result


def _load_gap_json(text: str, gap: GapAnalysis) -> None:
    from cfdx.io.schema import Severity

    try:
        data = json.loads(text)
    except json.JSONDecodeError as exc:
        gap.unsupported_nonblocking(
            "schema", "gap_report_json", f"Cannot parse stored gap report: {exc}"
        )
        return
    for item in data.get("findings", []):
        severity = Severity(item.get("severity", "unavailable"))
        gap.add(
            severity,
            item.get("category", "unknown"),
            item.get("feature", "unknown"),
            item.get("detail", ""),
            item.get("suggestion", ""),
        )


# --------------------------------------------------------------------------
# Registry
# --------------------------------------------------------------------------

_CONVERTER_REGISTRY: dict[str, type[SolverConverter]] = {
    "fluent": FluentConverter,
    "openfoam": OpenFoamConverter,
    "su2": Su2Converter,
    "starccm": StarCcmConverter,
    "saturne": SaturneConverter,
    "cfdx": CfdxConverter,
    "cfdx_h5": CfdxConverter,
    "cfdx_native": CfdxConverter,
}

_CONVERTER_ALIASES = {
    "ansys": "fluent",
    "ansys_fluent": "fluent",
    "open_foam": "openfoam",
    "foam": "openfoam",
    "su_2": "su2",
    "starccm+": "starccm",
    "star_ccm": "starccm",
    "star-ccm": "starccm",
    "starccm_plus": "starccm",
    "code_saturne": "saturne",
    "codesaturne": "saturne",
    "cfdx_hdf5": "cfdx",
}

_EXTENSION_MAP = {
    ".cas": "fluent",
    ".dat": "fluent",
    ".cas.h5": "fluent",
    ".dat.h5": "fluent",
    ".cas_h5": "fluent",
    ".dat_h5": "fluent",
    ".su2": "su2",
    ".cfg": "su2",
    ".sim": "starccm",
    ".xml": "saturne",
    ".py": "saturne",
    ".med": "saturne",
    ".cgns": "saturne",
    ".cfdx.h5": "cfdx",
    ".cfdx_h5": "cfdx",
}

_DIRECTORY_MARKERS = (
    (("constant", "polyMesh"), "openfoam"),
    (("system",), "openfoam"),
    (("mesh.su2",), "su2"),
    (("case.set",), "saturne"),
)


def available_converters() -> list[str]:
    """Names accepted by :func:`get_converter`."""
    return sorted(set(_CONVERTER_REGISTRY.keys()))


def register_converter(name: str, converter_class: type[SolverConverter]) -> None:
    """Register (or replace) a converter implementation."""
    _CONVERTER_REGISTRY[name.strip().lower()] = converter_class


def get_converter(solver_name: Union[str, type]) -> SolverConverter:
    """Return a converter instance for ``solver_name``.

    Args:
        solver_name: ``"fluent"``, ``"openfoam"``, ``"su2"``, ``"starccm"``,
            ``"saturne"`` or ``"cfdx"`` (case-insensitive, aliases allowed).
            A :class:`SolverConverter` subclass or instance is also accepted.

    Raises:
        ValueError: when the solver is unknown.
    """
    if isinstance(solver_name, SolverConverter):
        return solver_name
    if isinstance(solver_name, type) and issubclass(solver_name, SolverConverter):
        return solver_name()

    key = str(solver_name).strip().lower()
    key = _CONVERTER_ALIASES.get(key, key)
    converter_class = _CONVERTER_REGISTRY.get(key)
    if converter_class is None:
        raise ValueError(
            f"Unknown solver converter: {solver_name!r}. "
            f"Available: {', '.join(available_converters())}"
        )
    return converter_class()


def detect_solver(source: Union[str, Path]) -> Optional[str]:
    """Detect the converter key for ``source`` from extension or layout.

    Returns a registry key (``"fluent"``, ``"su2"``, ...) or ``None`` when the
    format cannot be determined.
    """
    path = Path(str(source))
    name = path.name.lower()

    for ext in sorted(_EXTENSION_MAP, key=len, reverse=True):
        if name.endswith(ext):
            return _EXTENSION_MAP[ext]

    if path.is_dir():
        for markers, key in _DIRECTORY_MARKERS:
            if all((path / m).exists() for m in markers):
                return key
        for pattern, key in (
            ("*.su2", "su2"),
            ("*.cas", "fluent"),
            ("*.sim", "starccm"),
            ("*.xml", "saturne"),
        ):
            if next(path.glob(pattern), None) is not None:
                return key
        return "openfoam"

    if path.suffix.lower() == ".h5":
        return "cfdx"

    return _sniff_content(path)


def _sniff_content(path: Path) -> Optional[str]:
    try:
        head = path.read_text(errors="replace")[:4096]
    except OSError:
        return None
    if "NDIME=" in head or head.lstrip().startswith("%"):
        return "su2"
    if head.startswith("(") or "COMPILED" in head:
        return "fluent"
    if "<Case" in head or "<case" in head.lower() or "case.set" in head:
        return "saturne"
    return None


def convert(
    source: Union[str, Path],
    output: Optional[Union[str, Path]] = None,
    solver: Optional[str] = None,
    write_reports: bool = True,
) -> ConversionResult:
    """Convert ``source`` to the CFDX canonical representation.

    The format is auto-detected unless ``solver`` is given.

    Args:
        source: Solver case file or case directory.
        output: Optional output file (``*.cfdx.h5`` by convention).
        solver: Force a specific converter name.
        write_reports: Write Markdown/JSON gap analysis next to the output.

    Raises:
        ValueError: when the format cannot be detected or the solver is unknown.
    """
    if solver is None:
        solver = detect_solver(source)
    if solver is None:
        raise ValueError(
            f"Cannot detect solver type for: {source}. "
            "Pass solver='fluent'|'openfoam'|'su2'|'starccm'|'saturne' explicitly."
        )
    converter = get_converter(solver)
    return converter.convert(source, output=output, write_reports=write_reports)


__all__ = [
    "ConversionResult",
    "SolverConverter",
    "AdapterConverter",
    "FluentConverter",
    "Su2Converter",
    "StarCcmConverter",
    "SaturneConverter",
    "OpenFoamConverter",
    "CfdxConverter",
    "get_converter",
    "convert",
    "detect_solver",
    "register_converter",
    "available_converters",
    "write_result",
]
