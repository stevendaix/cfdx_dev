"""Thin Python orchestration API for CFDX."""
from .case import Case, ExecutionConfig
from .session import CFDXSession, CaseNode, ChangeImpact, Checkpoint, SimulationState
from .tui import TuiRenderer
from .runner import ProcessResult, SolverRunner
from .metrics import SolverMetrics, SolverMetricsParser
from .renderer import Renderer, RenderObject, NullRenderer
from .execution import ExecutionController, ExecutionError
from .case_io import read_case, read_case_with_dat, save_case, save_case_with_dat
from .setup_model import MeshSelection, Parameter, ParameterType, SetupDiagnostic, typed_parameter
from .validation import ValidationReport, validate_case
from .mesh_model import MeshCatalog, MeshPatch, read_mesh_catalog
from .gui import CFDXMainWindow, create_application, launch

__all__ = [
    "Case", "ExecutionConfig", "CFDXSession", "CaseNode", "ChangeImpact",
    "Checkpoint", "SimulationState", "TuiRenderer", "ProcessResult",
    "SolverRunner", "SolverMetrics", "SolverMetricsParser",
    "Renderer", "RenderObject", "NullRenderer",
    "ExecutionController", "ExecutionError", "read_case", "read_case_with_dat",
    "save_case", "save_case_with_dat", "MeshSelection", "Parameter", "ParameterType",
    "SetupDiagnostic", "typed_parameter", "ValidationReport", "validate_case",
    "MeshCatalog", "MeshPatch", "read_mesh_catalog",
    "CFDXMainWindow", "create_application", "launch",
]
