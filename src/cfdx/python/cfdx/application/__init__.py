"""Headless application contracts shared by CFDX frontends."""

from .application import Application
from .commands import (
    Command,
    PauseSolver,
    ResumeSolver,
    RunSolver,
    SetNumericalOption,
    SetProperty,
    StopSolver,
    ValidateCase,
)
from .events import (
    ApplicationEvent,
    ApplicationStateChanged,
    CapabilitiesChanged,
    EventBus,
    ExecutionChanged,
    ProjectChanged,
    ResultsChanged,
    SelectionChanged,
    ValidationChanged,
)
from .state import (
    ApplicationState,
    ExecutionState,
    ProjectState,
    SelectionState,
    build_application_state,
)
from .workflow import WorkflowStatus, WorkflowStepState, build_workflow_state, workflow_children
from .properties import PropertyState, properties_for_selection, set_property
from ..setup_schema import SetupField, SetupDomain, SetupSchema, build_setup_schema, setup_fields_for_selection
from .run_center import RunCenterModel, RunCenterState
from .results import DisplayObject, ResultDataset, ResultFrameState, ResultsState, build_results_state
from .postprocessing import PostProcessingRegistry, PostProcessingSpec, build_post_processing_registry

__all__ = [
    "Application",
    "ApplicationEvent",
    "ApplicationState",
    "ApplicationStateChanged",
    "CapabilitiesChanged",
    "Command",
    "EventBus",
    "ExecutionChanged",
    "ExecutionState",
    "ProjectChanged",
    "ProjectState",
    "ResultsChanged",
    "PauseSolver",
    "ResumeSolver",
    "RunSolver",
    "RunCenterModel",
    "RunCenterState",
    "DisplayObject",
    "ResultDataset",
    "PostProcessingRegistry",
    "PostProcessingSpec",
    "build_post_processing_registry",
    "ResultFrameState",
    "ResultsState",
    "build_results_state",
    "SelectionChanged",
    "SelectionState",
    "SetNumericalOption",
    "SetProperty",
    "StopSolver",
    "ValidateCase",
    "ValidationChanged",
    "WorkflowStatus",
    "WorkflowStepState",
    "build_application_state",
    "build_workflow_state",
    "workflow_children",
    "PropertyState",
    "properties_for_selection",
    "set_property",
    "SetupField",
    "SetupDomain",
    "SetupSchema",
    "build_setup_schema",
    "setup_fields_for_selection",
]
