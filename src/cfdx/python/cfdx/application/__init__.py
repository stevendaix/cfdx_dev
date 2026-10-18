"""Headless application contracts shared by CFDX frontends."""

from ..setup_schema import (
    SetupDomain,
    SetupField,
    SetupSchema,
    build_setup_schema,
    setup_fields_for_selection,
)
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
from .postprocessing import (
    FieldDataset,
    PostProcessingRegistry,
    PostProcessingSpec,
    build_post_processing_registry,
    execute_derived,
)
from .properties import PropertyState, properties_for_selection, set_property
from .results import (
    DisplayObject,
    ResultDataset,
    ResultFrameState,
    ResultsState,
    build_results_state,
)
from .run_center import RunCenterModel, RunCenterState
from .state import (
    ApplicationState,
    ExecutionState,
    ProjectState,
    SelectionState,
    build_application_state,
)
from .workflow import (
    WorkflowStatus,
    WorkflowStepState,
    build_workflow_state,
    workflow_children,
)

__all__ = [
    "Application",
    "ApplicationEvent",
    "ApplicationState",
    "ApplicationStateChanged",
    "CapabilitiesChanged",
    "Command",
    "DisplayObject",
    "EventBus",
    "ExecutionChanged",
    "ExecutionState",
    "FieldDataset",
    "PauseSolver",
    "PostProcessingRegistry",
    "PostProcessingSpec",
    "ProjectChanged",
    "ProjectState",
    "PropertyState",
    "ResultDataset",
    "ResultFrameState",
    "ResultsChanged",
    "ResultsState",
    "ResumeSolver",
    "RunCenterModel",
    "RunCenterState",
    "RunSolver",
    "SelectionChanged",
    "SelectionState",
    "SetNumericalOption",
    "SetProperty",
    "SetupDomain",
    "SetupField",
    "SetupSchema",
    "StopSolver",
    "ValidateCase",
    "ValidationChanged",
    "WorkflowStatus",
    "WorkflowStepState",
    "build_application_state",
    "build_post_processing_registry",
    "build_results_state",
    "build_setup_schema",
    "build_workflow_state",
    "execute_derived",
    "properties_for_selection",
    "set_property",
    "setup_fields_for_selection",
    "workflow_children",
]
