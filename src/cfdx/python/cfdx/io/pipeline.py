"""CFDX I/O pipeline orchestrator.

Coordinates conversion of external solver data into the CFDX intermediate
representation. Selects the appropriate adapter based on file extension,
runs conversion, and exports GapAnalysis reports.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Optional, Union

from cfdx.io.adapters.fluent import FluentAdapter
from cfdx.io.adapters.su2 import Su2Adapter
from cfdx.io.adapters.starccm import StarCCMAdapter
from cfdx.io.adapters.saturne import SaturneAdapter
from cfdx.io.interfaces import ConversionResult
from cfdx.io.schema import SourceInfo
from cfdx.io.gap_analysis import GapAnalysisReport


_ADAPTER_REGISTRY = {
    "fluent": (FluentAdapter, (".cas",)),
    "su2": (Su2Adapter, (".su2", ".cfg")),
    "starccm": (StarCCMAdapter, (".sim",)),
    "saturne": (SaturneAdapter, (".xml", ".py")),
}


def detect_adapter(path: Union[str, Path]) -> Optional[str]:
    """Detect solver type from file extension and content."""
    p = Path(str(path))
    ext = p.suffix.lower()

    # Extension-based detection
    ext_map = {
        ".cas": "fluent",
        ".dat": "fluent",
        ".sim": "starccm",
    }
    if ext in ext_map:
        return ext_map[ext]

    if ext == ".su2":
        return "su2"
    if ext == ".xml":
        return "saturne"
    if ext == ".py":
        # Could be a Code_Saturne Python setup
        return "saturne"
    if ext == ".cfg":
        return "su2"

    # Content-based detection
    try:
        head = p.read_text(errors="replace")[:2048]
    except OSError:
        return None

    if head.startswith("(") and "COMPILED" in head:
        return "fluent"
    if "NDIME=" in head:
        return "su2"
    if p.is_file() and p.parent.name.lower() == "starccm":
        return "starccm"
    if "<Case" in head or "<case" in head.lower():
        return "saturne"
    if "case.set" in head:
        return "saturne"

    return None


def create_adapter(solver_name: str):
    """Create an adapter instance by solver name."""
    entry = _ADAPTER_REGISTRY.get(solver_name.lower())
    if entry is None:
        raise ValueError(f"Unknown solver: {solver_name}")
    adapter_cls = entry[0]
    return adapter_cls()


def convert(
    case_path: Union[str, Path],
    output_path: Optional[Union[str, Path]] = None,
    solver: Optional[str] = None,
    write_reports: bool = True,
) -> ConversionResult:
    """Convert an external solver case to CFDX IR.

    Args:
        case_path: Path to the solver case file or directory.
        output_path: Optional directory to write GapAnalysis reports and mesh.
        solver: Force a specific solver (auto-detect if None).
        write_reports: Whether to write Markdown/JSON gap reports.

    Returns:
        ConversionResult with mesh, setup, fields, and gap_report.
    """
    path = Path(str(case_path))

    if solver is None:
        solver = detect_adapter(path)
    if solver is None:
        raise ValueError(f"Cannot detect solver type for: {path}")

    adapter = create_adapter(solver)

    info = SourceInfo()
    adapter.detect_source(str(path), info)

    result = ConversionResult(source=info)
    success = adapter.convert(str(path), result)

    if output_path is not None:
        out_dir = Path(str(output_path))
        out_dir.mkdir(parents=True, exist_ok=True)

        if write_reports:
            report = GapAnalysisReport(result.gap_report)
            (out_dir / "gap_analysis.md").write_text(report.to_markdown())
            (out_dir / "gap_analysis.json").write_text(report.to_json())

        if result.mesh is not None:
            import meshio

            points = result.mesh.get("points", None)
            cells = result.mesh.get("cells", {})
            if points is not None:
                mesh = meshio.Mesh(
                    points=points,
                    cells=[(k, v) for k, v in cells.items()],
                )
                mesh.write(str(out_dir / "mesh.vtk"))

        if result.setup is not None:
            setup_dict = result.setup.model_dump(mode="json")
            with open(out_dir / "case_setup.json", "w") as f:
                json.dump(setup_dict, f, indent=2)

    return result


def get_available_adapters() -> list[str]:
    """Return list of available solver adapter names."""
    return list(_ADAPTER_REGISTRY.keys())
