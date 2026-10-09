"""Optional Qt mesh browser backed by the real MeshCatalog model."""
from __future__ import annotations
from PySide6.QtCore import Signal
from PySide6.QtWidgets import QLabel, QListWidget, QVBoxLayout, QWidget
from .mesh_model import MeshCatalog

class MeshBrowserPanel(QWidget):
    selection_changed = Signal(object)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.catalog: MeshCatalog | None = None
        self.status = QLabel("No mesh loaded")
        self.patches = QListWidget()
        self.patches.currentRowChanged.connect(self._select)
        layout=QVBoxLayout(self); layout.addWidget(self.status); layout.addWidget(self.patches)

    def set_catalog(self, catalog: MeshCatalog) -> None:
        self.catalog=catalog; self.patches.clear()
        for patch in catalog.patches:
            self.patches.addItem(f"{patch.name}  [{patch.type}]  faces={patch.face_count}")
        self.status.setText(f"{catalog.n_cells} cells / {catalog.n_faces} faces / {catalog.n_points} points")

    def clear(self) -> None:
        self.catalog=None; self.patches.clear(); self.status.setText("No mesh loaded")

    def _select(self, row: int) -> None:
        if self.catalog is not None and row >= 0:
            self.selection_changed.emit(self.catalog.patch(row).selection)
