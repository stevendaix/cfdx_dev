"""Optional Qt mesh browser backed by the real MeshCatalog model.

The panel is a regions -> patches -> faces tree.  Selections travel as
:class:`~cfdx.setup_model.MeshSelection` objects carrying the canonical
stable ids shared with the 3D view, so both widgets highlight the same
mesh object for one click.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QLabel,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .mesh_model import MeshCatalog, MeshPatch, MeshRegion
from .selection import face_selection_id
from .setup_model import MeshSelection

_FACE_PAGE_SIZE = 250


class MeshBrowserPanel(QWidget):
    """Tree browser over mesh regions, patches and faces."""

    selection_changed = Signal(object)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.catalog: MeshCatalog | None = None
        self.status = QLabel("No mesh loaded")
        self.patches = QTreeWidget()
        self.patches.setObjectName("workbench.mesh_browser.patches")
        self.patches.setHeaderLabel("Mesh")
        self.patches.itemSelectionChanged.connect(self._select)
        self.patches.itemExpanded.connect(self._populate_item_children)
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.patches)

    def set_catalog(self, catalog: MeshCatalog) -> None:
        self.catalog = catalog
        self.patches.clear()
        regions = catalog.regions
        if not regions or (len(regions) == 1 and regions[0].name == "default"):
            for patch in catalog.patches:
                self._add_patch_item(self.patches, patch)
        else:
            for index, region in enumerate(regions):
                self._add_region_item(region, index)
        self.status.setText(
            f"{catalog.n_cells} cells / {catalog.n_faces} faces / "
            f"{catalog.n_points} points"
        )

    def clear(self) -> None:
        self.catalog = None
        self.patches.clear()
        self.status.setText("No mesh loaded")

    def _add_region_item(self, region: MeshRegion, index: int) -> None:
        """Add one region node with its patches as children."""
        label = f"{region.name}  [{region.n_cells} cells]"
        item = QTreeWidgetItem(self.patches, [label])
        item.setData(
            0,
            Qt.ItemDataRole.UserRole,
            MeshSelection("region", index, region.name, region.stable_id),
        )
        for patch in region.patches:
            self._add_patch_item(item, patch)

    def _add_patch_item(
        self, parent: QTreeWidget | QTreeWidgetItem, patch: MeshPatch
    ) -> None:
        """Add one patch node whose face children load lazily."""
        item = QTreeWidgetItem(
            parent, [f"{patch.name}  [{patch.type}]  faces={patch.face_count}"]
        )
        item.setData(0, Qt.ItemDataRole.UserRole, patch.selection)
        item.setData(0, Qt.ItemDataRole.UserRole + 1, patch)
        if patch.face_count:
            item.setChildIndicatorPolicy(
                QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
            )

    def _populate_item_children(self, item: QTreeWidgetItem) -> None:
        """Create one bounded page of face children on demand."""
        page_end = item.data(0, Qt.ItemDataRole.UserRole + 2)
        if isinstance(page_end, int):
            parent = item.parent()
            patch = item.data(0, Qt.ItemDataRole.UserRole + 1)
            if parent is None or patch is None:
                return
            parent.removeChild(item)
            self._append_face_page(parent, patch, page_end)
            return
        patch = item.data(0, Qt.ItemDataRole.UserRole + 1)
        if patch is None or item.childCount():
            return
        self._append_face_page(item, patch, 0)

    def _append_face_page(
        self, item: QTreeWidgetItem, patch: MeshPatch, start: int
    ) -> None:
        end = min(start + _FACE_PAGE_SIZE, patch.face_count)
        for face_index in range(start, end):
            face_id = face_selection_id(patch.stable_id, face_index)
            child = QTreeWidgetItem(item, [f"Face {face_index}"])
            child.setData(
                0,
                Qt.ItemDataRole.UserRole,
                MeshSelection("face", face_index, patch.name, face_id),
            )
        item.setData(0, Qt.ItemDataRole.UserRole + 3, end)
        if end < patch.face_count:
            more = QTreeWidgetItem(
                item,
                [
                    f"Load faces {end}–{min(end + _FACE_PAGE_SIZE, patch.face_count) - 1}…"
                ],
            )
            more.setData(0, Qt.ItemDataRole.UserRole + 1, patch)
            more.setData(0, Qt.ItemDataRole.UserRole + 2, end)
            more.setChildIndicatorPolicy(
                QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
            )

    def select_stable_id(self, stable_id: str) -> MeshSelection | None:
        """Select a currently represented browser node by its canonical ID."""

        def find(item: QTreeWidgetItem) -> QTreeWidgetItem | None:
            selection = item.data(0, Qt.ItemDataRole.UserRole)
            if (
                isinstance(selection, MeshSelection)
                and selection.stable_id == stable_id
            ):
                return item
            for index in range(item.childCount()):
                match = find(item.child(index))
                if match is not None:
                    return match
            return None

        for index in range(self.patches.topLevelItemCount()):
            item = find(self.patches.topLevelItem(index))
            if item is not None:
                selection = item.data(0, Qt.ItemDataRole.UserRole)
                ancestor = item.parent()
                while ancestor is not None:
                    ancestor.setExpanded(True)
                    ancestor = ancestor.parent()
                signals_were_blocked = self.patches.blockSignals(True)
                try:
                    self.patches.setCurrentItem(item)
                    self.patches.scrollToItem(item)
                finally:
                    self.patches.blockSignals(signals_were_blocked)
                return selection
        return None

    def _select(self) -> None:
        items = self.patches.selectedItems()
        if not items:
            return
        selection = items[0].data(0, Qt.ItemDataRole.UserRole)
        if selection is not None:
            self.selection_changed.emit(selection)
