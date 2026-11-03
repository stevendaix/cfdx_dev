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

from .case_inspector import (\n    CaseFileSummary,\n    read_hdf5_dataset_sample,\n)


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
        self.preview.setPlainText("Select a dataset to inspect a bounded sample.")
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.tree, 3)
        layout.addWidget(self.preview, 1)
        self.tree.itemClicked.connect(self._preview_dataset)

    def clear(self) -> None:
        self.summary = None
        self.tree.clear()
        self.preview.setPlainText("Select a dataset to inspect a bounded sample.")
        self.status.setText("No CFDX file inspected")

    def _preview_dataset(self, item: QTreeWidgetItem, _column: int) -> None:
        """Display a bounded sample only when a dataset row is selected."""
        dataset_path = item.data(0, Qt.ItemDataRole.UserRole)
        if not dataset_path or self.summary is None:
            return
        try:
            sample = read_hdf5_dataset_sample(
                self.summary.path, dataset_path, max_values=32
            )
            values = ", ".join(sample.values) or "(empty dataset)"
            truncated = (
                f"\nshowing {len(sample.values)} of {sample.total_values:,} values"
                if sample.total_values > len(sample.values)
                else ""
            )
            preview = (
                f"{sample.path}\nshape: {sample.shape}\ndtype: {sample.dtype}"
                f"\nfirst values: {values}{truncated}"
            )
        except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
            preview = f"Unable to read dataset sample: {exc}"
        self.preview.setPlainText(preview)

    def set_summary(self, summary: CaseFileSummary) -> None:
        self.summary = summary
        self.tree.clear()
        self.preview.setPlainText(
            "Select a dataset to inspect a bounded sample."
        )
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
                display_value = (
                    str(value) if value is not None else "not declared"
                )
                if not summary.mesh_counts_complete and label != "Regions":
                    display_value = f"≥{display_value}"
                if label == "Regions" and summary.regions_truncated:
                    display_value = f"{len(summary.regions)}+"
                QTreeWidgetItem(mesh, [label, display_value])
            if not summary.mesh_counts_complete:
                QTreeWidgetItem(
                    mesh, ["Count status", "Partial: first 300 regions only"]
                )
            for region in summary.regions:
                QTreeWidgetItem(mesh, ["Region", region])
            if summary.regions_truncated:
                QTreeWidgetItem(mesh, ["Regions", "Listing truncated at 300 entries"])
        else:
            QTreeWidgetItem(mesh, ["Status", "No recognised mesh topology"])
        fields_truncated = (
            summary.scalar_fields_truncated or summary.vector_fields_truncated
        )
        field_total = (
            f"{summary.total_field_count}+"
            if fields_truncated
            else str(summary.total_field_count)
        )
        fields = QTreeWidgetItem(self.tree, ["Fields", field_total])
        for label, names in (
            ("Scalar", summary.scalar_fields),
            ("Vector", summary.vector_fields),
        ):
            truncated = (
                summary.scalar_fields_truncated
                if label == "Scalar"
                else summary.vector_fields_truncated
            )
            count_label = f"{len(names)}+" if truncated else str(len(names))
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
            size_label = (
                f"{dataset.size_bytes:,} bytes"
                if dataset.size_bytes is not None
                else "variable/unknown size"
            )
            chunk_label = (
                (
                    f"chunks={dataset.chunks}"
                    if dataset.chunks is not None
                    else "contiguous"
                )
            )
            compression_label = (
                f"compression={dataset.compression}"
                if dataset.compression is not None
                else "uncompressed"
            )
            details = (
                f"{dataset.shape} · {dataset.dtype} · {size_label} · "
                f"{chunk_label} · {compression_label}"
            )
            item = QTreeWidgetItem(datasets, [dataset.path, details])
            item.setToolTip(0, dataset.path)
            item.setToolTip(1, details)
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

