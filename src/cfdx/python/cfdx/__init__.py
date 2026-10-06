#!/usr/bin/env python3
"""CFDX Python package — solver I/O conversion and case orchestration.

Two surfaces used to live in two different trees that both installed the
package name ``cfdx``:

- ``src/cfdx/python/cfdx`` — solver I/O adapters, canonical numerical
  selection and Gap Analysis reporting (issue #425). This is the tree
  ``pyproject.toml`` packages, so it is what ``pip install cfdx`` shipped.
- ``python/cfdx`` — the orchestration API (session, runner, metrics,
  post-processing, GUI). Not packaged at all, and shadowed by the tree above
  because ``tests/python/conftest.py`` puts it first on ``sys.path``.

They are now one package. ``import cfdx`` exposes both.
"""

__version__ = "0.7.0"

from cfdx.io.interfaces import (
    SolverAdapter,
    ConversionResult,
    SourceInfo,
)
from cfdx.io.gap_analysis import (
    GapAnalysis,
    GapAnalysisReport,
    Severity,
    Finding,
)
from cfdx.io.schema import (
    CaseSetup,
    BoundarySpec,
    MaterialSpec,
    InitialCondition,
    NumericalScheme,
    MeshMetadata,
    BCType,
    BCValueType,
    CFDX_SCHEMA_VERSION,
)
from cfdx.io.mapping_rules import MappingRules

from .case import Case, ExecutionConfig
from .session import CFDXSession, CaseNode, ChangeImpact, Checkpoint, SimulationState
from .tui import TuiRenderer
from .runner import ProcessResult, SolverRunner
from .metrics import SolverMetrics, SolverMetricsParser
from .renderer import Renderer, RenderObject, NullRenderer
from .execution import ExecutionController, ExecutionError
from .case_io import read_case, read_case_with_dat, save_case, save_case_with_dat
from .project import Project
from .setup_model import MeshSelection, Parameter, ParameterType, SetupDiagnostic, typed_parameter\nfrom .setup_schema import SetupField, SetupDomain, SetupSchema, build_setup_schema, setup_fields_for_selection
from .validation import ValidationReport, validate_case
from .mesh_model import MeshCatalog, MeshPatch, read_mesh_catalog
from .physics_setup import PhysicsSpec, FieldSpec, physics_spec, default_parameters
from .boundary_setup import BoundaryFieldSpec, BoundaryTypeSpec, scalar_fields_for_boundary_type, boundary_defaults

# The GUI needs PySide6, and the 3-D view and result-series plumbing need
# pyvista/pyvistaqt. Those are optional extras, and the CI jobs that install
# only the base dependencies must keep importing ``cfdx`` without them, so the
# GUI entry points resolve on first access instead of at import time.
_LAZY_EXPORTS = {
    "CFDXMainWindow": "cfdx.gui",
    "create_application": "cfdx.gui",
    "launch": "cfdx.gui",
}

__all__ = [
    "SolverAdapter",
    "ConversionResult",
    "SourceInfo",
    "GapAnalysis",
    "Severity",
    "Finding",
    "CaseSetup",
    "BoundarySpec",
    "MaterialSpec",
    "InitialCondition",
    "NumericalScheme",
    "MeshMetadata",
    "BCType",
    "BCValueType",
    "MappingRules",
    "GapAnalysisReport",
    "CFDX_SCHEMA_VERSION",
    "Case", "ExecutionConfig", "CFDXSession", "CaseNode", "ChangeImpact",
    "Checkpoint", "SimulationState", "TuiRenderer", "ProcessResult",
    "SolverRunner", "SolverMetrics", "SolverMetricsParser",
    "Renderer", "RenderObject", "NullRenderer",
    "ExecutionController", "ExecutionError", "read_case", "read_case_with_dat",
    "save_case", "save_case_with_dat", "Project", "MeshSelection", "Parameter", "ParameterType",
    "SetupDiagnostic", "typed_parameter", "SetupField", "SetupDomain", "SetupSchema", "build_setup_schema", "setup_fields_for_selection", "ValidationReport", "validate_case",
    "MeshCatalog", "MeshPatch", "PhysicsSpec", "FieldSpec", "physics_spec",
    "default_parameters", "BoundaryFieldSpec", "BoundaryTypeSpec",
    "scalar_fields_for_boundary_type", "boundary_defaults",
    "CFDXMainWindow", "create_application", "launch",
]


def __getattr__(name: str):
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    return getattr(importlib.import_module(module_name), name)


def __dir__() -> list:
    return sorted(__all__)
