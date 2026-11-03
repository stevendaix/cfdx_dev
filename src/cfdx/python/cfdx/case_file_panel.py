"""Qt panel for metadata-first CFDX HDF5 inspection."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QLabel,
    QPlainTextEdit,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .case_inspector import CaseFileSummary, read_hdf5_dataset_sample, read_dataset_preview


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
        self.preview = QPlainTextEdit()
        self.preview.setObjectName("workbench.case_file.dataset_preview")
        self.preview.setReadOnly(True)
        self.preview.setMaximumBlockCount(40)
        self.preview.setPlainText("Select a numeric dataset to inspect a bounded sample.")
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.tree, 3)
        layout.addWidget(self.preview, 1)
        self.tree.itemClicked.connect(self._preview_dataset)

    def clear(self) -> None:
        self.summary = None
        self.tree.clear()
        self.preview.setPlainText("Select a numeric dataset to inspect a bounded sample.")
        self.status.setText("No CFDX file inspected")

    def _preview_dataset(self, item: QTreeWidgetItem, _column: int) -> None:
        """Display a bounded sample only when a dataset row is selected."""
        dataset_path = item.data(0, Qt.ItemDataRole.UserRole)
        if not dataset_path or self.summary is None:
            return
        try:
            preview = read_dataset_preview(self.summary.path, dataset_path)
        except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
            preview = f"Unable to read dataset preview: {exc}"
        self.preview.setPlainText(preview)

    def set_summary(self, summary: CaseFileSummary) -> None:
        self.summary = summary
        self.tree.clear()
        self.sample.clear()
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
                QTreeWidgetItem(
                    mesh,
                    [label, str(value) if value is not None else "not declared"],
                )
            for region in summary.regions:
                QTreeWidgetItem(mesh, ["Region", region])
            if summary.regions_truncated:
                QTreeWidgetItem(mesh, ["Regions", "Listing truncated at 300 entries"])
        else:
            QTreeWidgetItem(mesh, ["Status", "No recognised mesh topology"])
        fields = QTreeWidgetItem(self.tree, ["Fields", str(summary.total_field_count)])
        for label, names in (
            ("Scalar", summary.scalar_fields),
            ("Vector", summary.vector_fields),
        ):
            truncated = (
                summary.scalar_fields_truncated
                if label == "Scalar"
                else summary.vector_fields_truncated
            )
            count_label = f"≥{len(names)}" if truncated else str(len(names))
            group = QTreeWidgetItem(fields, [label, count_label])
            for name in names:
                QTreeWidgetItem(group, [name, ""])
            if truncated:
                QTreeWidgetItem(group, ["Field listing", "Truncated at 300 entries"])
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
        datasets = QTreeWidgetItem(
            self.tree, ["HDF5 datasets", str(len(summary.datasets))]
        )
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
        attributes = QTreeWidgetItem(
            self.tree, ["Attributes", str(len(summary.attributes))]
        )
        for name, value in summary.attributes:
            item = QTreeWidgetItem(attributes, [name, value])
            item.setToolTip(1, value)
        if summary.attributes_truncated:
            QTreeWidgetItem(
                attributes,
                ["Attribute listing", "Truncated at the configured metadata limit"],
            )
        self.tree.expandToDepth(1)

    def _dataset_selected(self, current, _previous) -> None:
        """Read at most 32 values from the selected HDF5 dataset."""
        dataset_path = (
            current.data(0, Qt.ItemDataRole.UserRole) if current is not None else None
        )
        if not dataset_path or self.summary is None:
            return
        try:
            sample = read_hdf5_dataset_sample(self.summary.path, dataset_path)
        except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
            self.sample.setPlainText(f"Sample unavailable: {exc}")
            return
        values = ", ".join(sample.values) if sample.values else "<empty dataset>"
        if len(values) > 4000:
            values = values[:3997] + "..."
        self.sample.setPlainText(
            f"{sample.path}\n"
            f"Shape: {sample.shape} · dtype: {sample.dtype}\n"
            f"Showing {len(sample.values)} of {sample.total_values:,} values\n"
            f"{values}"
        )
