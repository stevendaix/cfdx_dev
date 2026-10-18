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

from cfdx.io.gap_analysis import (
    Finding,
    GapAnalysis,
    GapAnalysisReport,
    Severity,
)
from cfdx.io.interfaces import (
    ConversionResult,
    SolverAdapter,
    SourceInfo,
)
from cfdx.io.mapping_rules import MappingRules
from cfdx.io.schema import (
    CFDX_SCHEMA_VERSION,
    BCType,
    BCValueType,
    BoundarySpec,
    CaseSetup,
    InitialCondition,
    MaterialSpec,
    MeshMetadata,
    NumericalScheme,
)

from .boundary_setup import (
    BoundaryFieldSpec,
    BoundaryTypeSpec,
    boundary_defaults,
    scalar_fields_for_boundary_type,
)
from .case import Case, ExecutionConfig
from .case_io import read_case, read_case_with_dat, save_case, save_case_with_dat
from .execution import ExecutionController, ExecutionError
from .mesh_model import MeshCatalog, MeshPatch, read_mesh_catalog
from .metrics import SolverMetrics, SolverMetricsParser
from .physics_setup import FieldSpec, PhysicsSpec, default_parameters, physics_spec
from .project import Project
from .renderer import NullRenderer, Renderer, RenderObject
from .runner import ProcessResult, SolverRunner
from .session import CaseNode, CFDXSession, ChangeImpact, Checkpoint, SimulationState
from .setup_model import (
    MeshSelection,
    Parameter,
    ParameterType,
    SetupDiagnostic,
    typed_parameter,
)
from .setup_schema import (
    SetupDomain,
    SetupField,
    SetupSchema,
    build_setup_schema,
    setup_fields_for_selection,
)
from .tui import TuiRenderer
from .validation import ValidationReport, validate_case

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
    "CFDX_SCHEMA_VERSION",
    "BCType",
    "BCValueType",
    "BoundaryFieldSpec",
    "BoundarySpec",
    "BoundaryTypeSpec",
    "CFDXMainWindow",
    "CFDXSession",
    "Case",
    "CaseNode",
    "CaseSetup",
    "ChangeImpact",
    "Checkpoint",
    "ConversionResult",
    "ExecutionConfig",
    "ExecutionController",
    "ExecutionError",
    "FieldSpec",
    "Finding",
    "GapAnalysis",
    "GapAnalysisReport",
    "InitialCondition",
    "MappingRules",
    "MaterialSpec",
    "MeshCatalog",
    "MeshMetadata",
    "MeshPatch",
    "MeshSelection",
    "NullRenderer",
    "NumericalScheme",
    "Parameter",
    "ParameterType",
    "PhysicsSpec",
    "ProcessResult",
    "Project",
    "RenderObject",
    "Renderer",
    "SetupDiagnostic",
    "SetupDomain",
    "SetupField",
    "SetupSchema",
    "Severity",
    "SimulationState",
    "SolverAdapter",
    "SolverMetrics",
    "SolverMetricsParser",
    "SolverRunner",
    "SourceInfo",
    "TuiRenderer",
    "ValidationReport",
    "boundary_defaults",
    "build_setup_schema",
    "create_application",
    "default_parameters",
    "launch",
    "physics_spec",
    "read_case",
    "read_case_with_dat",
    "save_case",
    "save_case_with_dat",
    "scalar_fields_for_boundary_type",
    "setup_fields_for_selection",
    "typed_parameter",
    "validate_case",
]


def __getattr__(name: str):
    module_name = _LAZY_EXPORTS.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    import importlib

    return getattr(importlib.import_module(module_name), name)


def __dir__() -> list:
    return sorted(__all__)
