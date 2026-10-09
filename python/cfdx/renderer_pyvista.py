"""Optional PyVista renderer adapter for CFDX result datasets."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from .renderer import RenderObject


class PyVistaRenderer:
    """Thin adapter; PyVista remains an optional visualization dependency."""

    def __init__(self) -> None:
        try:
            import pyvista as pv
        except ImportError as exc:
            raise RuntimeError("PyVista is required for 3D rendering") from exc
        self._pv = pv
        self.dataset: Any | None = None
        self.source: str | None = None
        self.selected: str | None = None

    def load(self, source: str) -> tuple[RenderObject, ...]:
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(source)
        self.dataset = self._pv.read(path)
        self.source = str(path)
        return (RenderObject("dataset", path.name),)

    def clear(self) -> None:
        self.dataset = None
        self.source = None
        self.selected = None

    def select(self, object_id: str) -> None:
        if object_id != "dataset" or self.dataset is None:
            raise KeyError(object_id)
        self.selected = object_id

    def slice(self, **kwargs: Any) -> Any:
        if self.dataset is None:
            raise RuntimeError("no dataset loaded")
        return self.dataset.slice(**kwargs)

    def contour(self, isosurfaces: Any, scalars: str) -> Any:
        if self.dataset is None:
            raise RuntimeError("no dataset loaded")
        return self.dataset.contour(isosurfaces=isosurfaces, scalars=scalars)

    def glyph(self, **kwargs: Any) -> Any:
        if self.dataset is None:
            raise RuntimeError("no dataset loaded")
        return self.dataset.glyph(**kwargs)
