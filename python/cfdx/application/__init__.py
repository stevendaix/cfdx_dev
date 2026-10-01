"""Headless application contracts shared by CFDX frontends."""

from .application import Application
from .commands import (
    Command,
    PauseSolver,
    ResumeSolver,
    RunSolver,
    SetNumericalOption,
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
    "SelectionChanged",
    "SelectionState",
    "SetNumericalOption",
    "StopSolver",
    "ValidateCase",
    "ValidationChanged",
    "build_application_state",
]
