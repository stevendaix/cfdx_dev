"""Backend-independent rendering contracts for CFDX post-processing."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, Sequence


@dataclass(frozen=True)
class RenderObject:
    """Stable object identity exposed to a view."""

    object_id: str
    label: str


class Renderer(Protocol):
    """Minimal renderer contract usable by GUI and headless clients."""

    def load(self, source: str) -> Sequence[RenderObject]:
        """Load a dataset and return stable selectable objects."""

    def clear(self) -> None:
        """Remove the current scene."""

    def select(self, object_id: str) -> None:
        """Select an object by stable logical identifier."""


class NullRenderer:
    """Headless renderer useful for tests and non-visual execution."""

    def __init__(self) -> None:
        self.objects: tuple[RenderObject, ...] = ()
        self.selected: str | None = None

    def load(self, source: str) -> tuple[RenderObject, ...]:
        self.objects = (RenderObject("dataset", source),)
        return self.objects

    def clear(self) -> None:
        self.objects = ()

    def select(self, object_id: str) -> None:
        if not any(obj.object_id == object_id for obj in self.objects):
            raise KeyError(object_id)
        self.selected = object_id
