"""Application events independent from Qt."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Generic, TypeVar

from .state import ApplicationState


@dataclass(frozen=True)
class ApplicationEvent:
    state: ApplicationState


@dataclass(frozen=True)
class ApplicationStateChanged(ApplicationEvent):
    pass


@dataclass(frozen=True)
class ProjectChanged(ApplicationEvent):
    pass


@dataclass(frozen=True)
class SelectionChanged(ApplicationEvent):
    pass


@dataclass(frozen=True)
class ValidationChanged(ApplicationEvent):
    pass


@dataclass(frozen=True)
class ExecutionChanged(ApplicationEvent):
    pass


@dataclass(frozen=True)
class ResultsChanged(ApplicationEvent):
    pass


@dataclass(frozen=True)
class CapabilitiesChanged(ApplicationEvent):
    pass


EventT = TypeVar("EventT", bound=ApplicationEvent)
Subscriber = Callable[[EventT], None]


class EventBus:
    """Small synchronous event bus suitable for GUI adapters and tests."""

    def __init__(self) -> None:
        self._subscribers: dict[type[ApplicationEvent], list[Subscriber]] = {}

    def subscribe(self, event_type: type[EventT], callback: Subscriber[EventT]) -> None:
        self._subscribers.setdefault(event_type, []).append(callback)

    def unsubscribe(self, event_type: type[EventT], callback: Subscriber[EventT]) -> None:
        callbacks = self._subscribers.get(event_type, [])
        if callback in callbacks:
            callbacks.remove(callback)

    def publish(self, event: EventT) -> None:
        for callback in tuple(self._subscribers.get(type(event), ())):
            callback(event)
        if type(event) is not ApplicationStateChanged:
            for callback in tuple(self._subscribers.get(ApplicationStateChanged, ())):
                callback(event)
