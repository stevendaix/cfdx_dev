"""Optional PyVistaQt 3D view for the CFDX GUI."""
from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    from PySide6.QtWidgets import QVBoxLayout, QWidget
except ImportError:
    QWidget = None


class PyVistaQtView(QWidget if QWidget is not None else object):
    """Qt container around pyvistaqt for CFDX result visualization."""

    def __init__(self, parent: Any = None) -> None:
        if QWidget is None:
            raise RuntimeError("PySide6 is required for the 3D GUI")
        super().__init__(parent)
        try:
            import pyvista as pv
            from pyvistaqt import QtInteractor
        except ImportError as exc:
            raise RuntimeError("PyVista and pyvistaqt are required for the 3D GUI") from exc
        self._pv = pv
        self.plotter = QtInteractor(self)
        layout = QVBoxLayout(self)
        layout.addWidget(self.plotter)

    def load(self, source: str) -> None:
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(source)
        self.plotter.clear()
        self.plotter.add_mesh(self._pv.read(path))
        self.plotter.reset_camera()
        self.plotter.render()

    def clear(self) -> None:
        self.plotter.clear()

    def closeEvent(self, event) -> None:
        self.plotter.close()
        super().closeEvent(event)
