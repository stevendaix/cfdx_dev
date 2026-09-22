"""Optional PyVistaQt 3D view for the CFDX GUI."""
from __future__ import annotations

from pathlib import Path
from typing import Any

try:
    from PySide6.QtWidgets import QVBoxLayout, QWidget
except ImportError:
    QWidget = None


class PyVistaQtView(QWidget if QWidget is not None else object):
    """Qt container around pyvistaqt with the CFDX renderer contract."""

    def __init__(self, parent: Any = None) -> None:
        if QWidget is None:
            raise RuntimeError("PySide6 is required for the 3D GUI")
        super().__init__(parent)
        try:
            from pyvistaqt import QtInteractor
        except ImportError as exc:
            raise RuntimeError("pyvistaqt is required for the 3D GUI") from exc
        self.plotter = QtInteractor(self)
        layout = QVBoxLayout(self)
        layout.addWidget(self.plotter)

    def load(self, source: str) -> None:
        path = Path(source)
        if not path.exists():
            raise FileNotFoundError(source)
        self.plotter.clear()
        self.plotter.add_mesh(self.plotter.reader(path).read())
        self.plotter.reset_camera()
        self.plotter.render()

    def clear(self) -> None:
        self.plotter.clear()

    def closeEvent(self, event) -> None:
        self.plotter.close()
        super().closeEvent(event)
