"""Qt panel for metadata-first CFDX HDF5 inspection."""
from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QLabel, QTreeWidget, QTreeWidgetItem, QVBoxLayout, QWidget

from .case_inspector import CaseFileSummary


class CaseFilePanel(QWidget):
    """Display case structure and dataset metadata without loading arrays."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.summary: CaseFileSummary | None = None
        self.status = QLabel("No CFDX file inspected")
        self.tree = QTreeWidget()
        self.tree.setObjectName("workbench.case_file.tree")
        self.tree.setHeaderLabels(["Case file", "Details"])
        self.tree.setColumnWidth(0, 260)
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.tree)

    def clear(self) -> None:
        self.summary = None
        self.tree.clear()
        self.status.setText("No CFDX file inspected")

    def set_summary(self, summary: CaseFileSummary) -> None:
        self.summary = summary
        self.tree.clear()
        state = "CFDX" if summary.is_cfdx else "HDF5 (unrecognised CFDX layout)"
        self.status.setText(f"{summary.path.name} · {state}")
        root = QTreeWidgetItem(self.tree, ["File", str(summary.path)])
        QTreeWidgetItem(root, ["Format", summary.format or "not declared"])
        QTreeWidgetItem(root, ["Schema", summary.schema_version or "not declared"])
        QTreeWidgetItem(
            root,
            [
                "Project configuration",
                "available" if summary.has_project_config else "not present",
            ],
        )
        QTreeWidgetItem(
            root,
            [
                "Imported case metadata",
                "available"
                if summary.has_case_config and not summary.has_project_config
                else "not present",
            ],
        )
        mesh = QTreeWidgetItem(self.tree, ["Mesh", ""])
        if summary.has_mesh:
            for label, value in (
                ("Points", summary.n_points),
                ("Faces", summary.n_faces),
                ("Cells", summary.n_cells),
                ("Regions", len(summary.regions) if summary.regions else 1),
            ):
                QTreeWidgetItem(mesh, [label, str(value) if value is not None else "not declared"])
            for region in summary.regions:
                QTreeWidgetItem(mesh, ["Region", region])
        else:
            QTreeWidgetItem(mesh, ["Status", "No recognised mesh topology"])
        fields = QTreeWidgetItem(self.tree, ["Fields", str(summary.total_field_count)])
        for label, names in (("Scalar", summary.scalar_fields), ("Vector", summary.vector_fields)):
            group = QTreeWidgetItem(fields, [label, str(len(names))])
            for name in names:
                QTreeWidgetItem(group, [name, ""])
        if summary.case_name:
            QTreeWidgetItem(self.tree, ["Case name", summary.case_name])
        for label, names in (
            ("Physics", summary.physics),
            ("Numerics", summary.numerics),
            ("Boundaries", summary.boundaries),
        ):
            group = QTreeWidgetItem(self.tree, [label, str(len(names))])
            for name in names:
                QTreeWidgetItem(group, [name, ""])
        datasets = QTreeWidgetItem(self.tree, ["HDF5 datasets", str(len(summary.datasets))])
        for dataset in summary.datasets:
            item = QTreeWidgetItem(
                datasets,
                [
                    dataset.path,
                    f"{dataset.shape} · {dataset.dtype} · {dataset.size_bytes:,} bytes",
                ],
            )
            item.setToolTip(0, dataset.path)
            item.setData(0, Qt.ItemDataRole.UserRole, dataset.path)
        if summary.datasets_truncated:
            QTreeWidgetItem(
                datasets,
                ["Dataset listing", "Truncated at the configured metadata limit"],
            )
        attributes = QTreeWidgetItem(self.tree, ["Attributes", str(len(summary.attributes))])
        for name, value in summary.attributes:
            item = QTreeWidgetItem(attributes, [name, value])
            item.setToolTip(1, value)
        self.tree.expandToDepth(1)
