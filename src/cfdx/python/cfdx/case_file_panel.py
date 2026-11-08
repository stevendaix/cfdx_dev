"""Qt panel for metadata-first CFDX HDF5 inspection."""

from __future__ import annotations

import json

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Qt, Signal
from PySide6.QtWidgets import (
    QApplication,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .case_inspector import (
    CaseFileSummary,
    read_hdf5_dataset_sample,
    read_hdf5_dataset_statistics,
)
from .mesh_quality import (
    SURFACE_QUALITY_MAX_FACES,
    analyze_cfdx_surface_file,
)
from .volume_skewness_view import VolumeSkewnessDialog
from .mesh_health import analyze_mesh_health_file, mesh_health_text


class _QualityTaskSignals(QObject):
    ready = Signal(int, object, object)


class _QualityTask(QRunnable):
    """Analyze numeric topology in a worker, never on the Qt GUI thread."""

    def __init__(self, job_id: int, path, signature) -> None:
        super().__init__()
        self.job_id = job_id
        self.path = path
        self.signature = signature
        self.signals = _QualityTaskSignals()

    def run(self) -> None:
        try:
            result = analyze_mesh_health_file(
                self.path, expected_signature=self.signature,
            )
        except (OSError, RuntimeError, ValueError, TypeError, KeyError) as exc:
            self.signals.ready.emit(self.job_id, None, str(exc))
        else:
            self.signals.ready.emit(self.job_id, result, None)


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
        self.filter_input = QLineEdit()
        self.filter_input.setObjectName("workbench.case_file.filter")
        self.filter_input.setPlaceholderText(
            "Filter datasets, fields, regions and metadata…"
        )
        self.filter_input.setClearButtonEnabled(True)
        self.expand_button = QPushButton("Expand all")
        self.expand_button.setObjectName("workbench.case_file.expand_all")
        self.collapse_button = QPushButton("Collapse all")
        self.collapse_button.setObjectName("workbench.case_file.collapse_all")
        self.copy_path_button = QPushButton("Copy dataset path")
        self.copy_path_button.setObjectName("workbench.case_file.copy_dataset_path")
        self.copy_path_button.setEnabled(False)
        self.statistics_button = QPushButton("Compute statistics")
        self.statistics_button.setObjectName("workbench.case_file.compute_statistics")
        self.statistics_button.setEnabled(False)
        self.mesh_quality_button = QPushButton("Analyze surface quality")
        self.mesh_quality_button.setObjectName("workbench.case_file.surface_quality")
        self.mesh_quality_button.setEnabled(False)
        self.volume_quality_button = QPushButton("Inspect volume skewness")
        self.volume_quality_button.setObjectName("workbench.case_file.volume_skewness")
        self.volume_quality_button.setEnabled(False)
        self._volume_quality_dialog = None
        self._mesh_health_snapshot = None
        self._quality_job_id = 0
        self._quality_tasks: dict[int, _QualityTask] = {}
        self.volume_health_button = QPushButton("Mesh health dashboard")
        self.volume_health_button.setObjectName("workbench.case_file.volume_health")
        self.volume_health_button.setEnabled(False)
        self.mesh_health_export = QPushButton("Export mesh audit")
        self.mesh_health_export.setObjectName("workbench.case_file.export_mesh_health")
        self.mesh_health_export.setEnabled(False)
        self.preview = QPlainTextEdit()
        self.preview.setObjectName("workbench.case_file.dataset_preview")
        self.preview.setReadOnly(True)
        self.preview.setMaximumBlockCount(40)
        self.preview.setPlainText("Select a dataset to inspect a bounded sample.")
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.filter_input)
        controls = QHBoxLayout()
        controls.addWidget(self.expand_button)
        controls.addWidget(self.collapse_button)
        controls.addStretch(1)
        controls.addWidget(self.mesh_quality_button)
        controls.addWidget(self.volume_quality_button)
        controls.addWidget(self.volume_health_button)
        controls.addWidget(self.mesh_health_export)
        controls.addWidget(self.statistics_button)
        controls.addWidget(self.copy_path_button)
        layout.addLayout(controls)
        layout.addWidget(self.tree, 3)
        layout.addWidget(self.preview, 1)
        self.tree.itemClicked.connect(self._preview_dataset)
        self.tree.currentItemChanged.connect(self._preview_current_item)
        self.tree.itemSelectionChanged.connect(self._update_selection_actions)
        self.filter_input.textChanged.connect(self._filter_tree)
        self.expand_button.clicked.connect(self.tree.expandAll)
        self.collapse_button.clicked.connect(self.tree.collapseAll)
        self.copy_path_button.clicked.connect(self.copy_selected_dataset_path)
        self.statistics_button.clicked.connect(self.compute_selected_statistics)
        self.mesh_quality_button.clicked.connect(self.compute_surface_quality)
        self.volume_quality_button.clicked.connect(self.open_volume_skewness)
        self.volume_health_button.clicked.connect(self.compute_volume_health)
        self.mesh_health_export.clicked.connect(self.export_mesh_health)

    def clear(self) -> None:
        self.summary = None
        self.filter_input.clear()
        self.tree.clear()
        self.copy_path_button.setEnabled(False)
        self.statistics_button.setEnabled(False)
        self.mesh_quality_button.setEnabled(False)
        self.volume_quality_button.setEnabled(False)
        self.volume_health_button.setEnabled(False)
        self.mesh_health_export.setEnabled(False)
        self._mesh_health_snapshot = None
        self._quality_job_id += 1
        self.preview.setPlainText("Select a dataset to inspect a bounded sample.")
        self.status.setText("No CFDX file inspected")

    def _filter_tree(self, text: str) -> None:
        """Filter the metadata tree while retaining matching ancestor nodes."""
        query = text.strip().casefold()
        for index in range(self.tree.topLevelItemCount()):
            self._filter_item(self.tree.topLevelItem(index), query)

    def _filter_item(self, item: QTreeWidgetItem, query: str) -> bool:
        own_match = not query or any(
            query in item.text(column).casefold() for column in range(2)
        )
        child_match = False
        for index in range(item.childCount()):
            child_match = self._filter_item(item.child(index), query) or child_match
        visible = own_match or child_match
        item.setHidden(not visible)
        return visible

    def _preview_current_item(
        self, current: QTreeWidgetItem | None, _previous: QTreeWidgetItem | None
    ) -> None:
        """Preview datasets selected with either the mouse or keyboard."""
        if current is not None:
            self._preview_dataset(current, 0)

    def _update_selection_actions(self) -> None:
        item = self.tree.currentItem()
        dataset_path = (
            item.data(0, Qt.ItemDataRole.UserRole) if item is not None else None
        )
        self.copy_path_button.setEnabled(bool(dataset_path))
        self.statistics_button.setEnabled(bool(dataset_path))

    def compute_selected_statistics(self) -> None:
        """Compute bounded streaming statistics for the selected dataset."""
        item = self.tree.currentItem()
        dataset_path = (
            item.data(0, Qt.ItemDataRole.UserRole) if item is not None else None
        )
        if not dataset_path or self.summary is None:
            return
        try:
            statistics = read_hdf5_dataset_statistics(
                self.summary.path, str(dataset_path)
            )

            def display(value: float | None) -> str:
                return "n/a" if value is None else f"{value:.8g}"

            report = (
                f"{statistics.path} — streaming statistics\n"
                f"shape: {statistics.shape} · dtype: {statistics.dtype}\n"
                f"values: {statistics.value_count:,} · "
                f"finite: {statistics.finite_count:,} · "
                f"non-finite: {statistics.non_finite_count:,}\n"
                f"min: {display(statistics.minimum)} · "
                f"max: {display(statistics.maximum)}\n"
                f"mean: {display(statistics.mean)} · "
                f"RMS: {display(statistics.rms)}\n"
                f"scanned: {statistics.scanned_bytes:,} bytes "
                "(uniform value statistics; no cell-volume weighting)"
            )
        except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
            report = f"Unable to compute dataset statistics: {exc}"
        self.preview.setPlainText(report)

    def compute_surface_quality(self) -> None:
        """Run explicit, bounded geometric diagnostics on the inspected surface."""
        if self.summary is None or not self.summary.preview_supported:
            return
        try:
            report = analyze_cfdx_surface_file(
                self.summary.path,
                expected_signature=self.summary.file_signature,
            )

            def display(value: float | None) -> str:
                return "n/a" if value is None else f"{value:.8g}"

            worst_faces = (
                ", ".join(str(face_id) for face_id in report.worst_warpage_face_ids)
                or "none"
            )
            report_text = (
                "Surface geometry diagnostics (not volume-cell quality)\n"
                f"faces: {report.face_count:,} · valid: {report.valid_face_count:,} · "
                f"degenerate: {len(report.degenerate_face_ids):,}\n"
                f"area min/max: {display(report.face_area_min)} / "
                f"{display(report.face_area_max)}\n"
                f"edge-length ratio min/p95/max: "
                f"{display(report.edge_length_ratio_min)} / "
                f"{display(report.edge_length_ratio_p95)} / "
                f"{display(report.edge_length_ratio_max)}\n"
                f"warpage ratio p95/max: {display(report.warpage_ratio_p95)} / "
                f"{display(report.warpage_ratio_max)}\n"
                f"zero-length-edge faces: {len(report.zero_length_edge_face_ids):,}\n"
                f"worst warpage face IDs: {worst_faces}\n"
                "Volume-cell skewness, non-orthogonality and Jacobian are not "
                "available from this surface-only topology."
            )
        except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
            report_text = f"Unable to analyze surface quality: {exc}"
        self.preview.setPlainText(report_text)

    def open_volume_skewness(self) -> None:
        """Open the interactive cut-plane view for bounded volume topology."""
        if self.summary is None or not self.summary.preview_supported:
            return
        try:
            dialog = VolumeSkewnessDialog(
                self.summary.path,
                expected_signature=self.summary.file_signature,
                parent=self,
            )
            self._volume_quality_dialog = dialog
            dialog.show()
        except (OSError, RuntimeError, TypeError, ValueError, KeyError) as exc:
            self.preview.setPlainText(f"Unable to open volume skewness view: {exc}")

    def compute_volume_health(self) -> None:
        """Opt-in, bounded computation with stale-selection protection."""
        summary = self.summary
        if summary is None or not self.volume_health_button.isEnabled():
            return
        self._quality_job_id += 1
        job_id = self._quality_job_id
        self._mesh_health_snapshot = None
        self.mesh_health_export.setEnabled(False)
        self.volume_health_button.setEnabled(False)
        self.preview.setPlainText("Computing bounded volume quality in a worker…")
        task = _QualityTask(job_id, summary.path, summary.file_signature)
        self._quality_tasks[job_id] = task
        task.signals.ready.connect(self._volume_health_ready)
        QThreadPool.globalInstance().start(task)

    def _volume_health_ready(self, job_id: int, snapshot, error) -> None:
        self._quality_tasks.pop(job_id, None)
        if job_id != self._quality_job_id or self.summary is None:
            return  # File or selection replaced before this worker completed.
        self.volume_health_button.setEnabled(self.volume_quality_button.isEnabled())
        if error is not None:
            self.preview.setPlainText(f"Mesh health analysis failed: {error}")
            return
        if (snapshot.source_file != str(self.summary.path)
                or snapshot.file_signature != self.summary.file_signature):
            self.preview.setPlainText("Mesh health result ignored: inspected file changed")
            return
        self._mesh_health_snapshot = snapshot
        self.mesh_health_export.setEnabled(True)
        self.preview.setPlainText(mesh_health_text(snapshot))

    def export_mesh_health(self) -> None:
        snapshot = self._mesh_health_snapshot
        if snapshot is None or self.summary is None:
            return
        try:
            source = self.summary.path
            stat = source.stat()
            identity = (stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns)
            if identity != snapshot.file_signature:
                raise RuntimeError("source mesh has changed; recalculate the quality report")
        except (OSError, RuntimeError) as exc:
            self.preview.setPlainText(f"Mesh audit export rejected: {exc}")
            self.mesh_health_export.setEnabled(False)
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Export mesh quality evidence", "cfdx-mesh-health.json", "JSON (*.json)"
        )
        if not path:
            return
        try:
            import json
            from pathlib import Path
            Path(path).write_text(
                json.dumps(snapshot.to_dict(), indent=2, sort_keys=True, allow_nan=False) + "\n",
                encoding="utf-8",
            )
        except (OSError, ValueError, TypeError) as exc:
            self.preview.setPlainText(f"Mesh audit export failed: {exc}")

    def copy_selected_dataset_path(self) -> None:
        """Copy the selected HDF5 dataset path without reading its payload."""
        item = self.tree.currentItem()
        dataset_path = (
            item.data(0, Qt.ItemDataRole.UserRole) if item is not None else None
        )
        if dataset_path:
            QApplication.clipboard().setText(str(dataset_path))

    def _preview_dataset(self, item: QTreeWidgetItem, _column: int) -> None:
        """Display a bounded sample only when a dataset row is selected."""
        dataset_path = item.data(0, Qt.ItemDataRole.UserRole)
        if not dataset_path or self.summary is None:
            self.preview.setPlainText("Select a dataset to inspect a bounded sample.")
            self._update_selection_actions()
            return
        try:
            sample = read_hdf5_dataset_sample(
                self.summary.path, dataset_path, max_values=32
            )
            values = ", ".join(sample.values) or "(empty dataset)"
            if len(values) > 4000:
                values = values[:3997] + "..."
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
        self._quality_job_id += 1
        self._mesh_health_snapshot = None
        self.mesh_health_export.setEnabled(False)
        self.filter_input.clear()
        self.tree.clear()
        self.copy_path_button.setEnabled(False)
        self.statistics_button.setEnabled(False)
        self.mesh_quality_button.setEnabled(
            summary.preview_supported
            and summary.n_faces is not None
            and summary.n_faces <= SURFACE_QUALITY_MAX_FACES
        )
        dataset_paths = {dataset.path for dataset in summary.datasets}
        self.volume_quality_button.setEnabled(
            summary.preview_supported
            and {
                "/" + name
                for name in (
                    "cell_faces",
                    "cell_offsets",
                    "owner",
                    "neighbour",
                    "face_vertices",
                    "face_offsets",
                    "points",
                )
            }.issubset(dataset_paths)
        )
        self.volume_health_button.setEnabled(self.volume_quality_button.isEnabled())
        self.preview.setPlainText("Select a dataset to inspect a bounded sample.")
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
                display_value = str(value) if value is not None else "not declared"
                if (
                    summary.regions_truncated
                    and label != "Regions"
                    and value is not None
                ):
                    display_value = f"≥{display_value}"
                if label == "Regions" and summary.regions_truncated:
                    display_value = (
                        f"{len(summary.regions)}+"
                        if len(summary.regions) == 300
                        else "incomplete"
                    )
                QTreeWidgetItem(mesh, [label, display_value])
            if not summary.mesh_counts_complete:
                count_note = (
                    "Region catalog incomplete; totals may be lower bounds"
                    if summary.regions_truncated
                    else "Incomplete topology count metadata"
                )
                QTreeWidgetItem(mesh, ["Count status", count_note])
            for region in summary.regions:
                QTreeWidgetItem(mesh, ["Region", region])
            if summary.regions_truncated:
                QTreeWidgetItem(
                    mesh,
                    [
                        "Regions",
                        "Listing incomplete (safety cap or non-internal links)",
                    ],
                )
        else:
            QTreeWidgetItem(mesh, ["Status", "No recognised mesh topology"])
        fields_truncated = (
            summary.scalar_fields_truncated or summary.vector_fields_truncated
        )
        field_count_is_lower_bound = (
            summary.scalar_fields_truncated and len(summary.scalar_fields) == 300
        ) or (summary.vector_fields_truncated and len(summary.vector_fields) == 300)
        field_total = (
            f"{summary.total_field_count}+"
            if fields_truncated and field_count_is_lower_bound
            else "incomplete"
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
            count_label = (
                f"{len(names)}+"
                if truncated and len(names) == 300
                else "incomplete"
                if truncated
                else str(len(names))
            )
            group = QTreeWidgetItem(fields, [label, count_label])
            for name in names:
                QTreeWidgetItem(group, [name, ""])
            if truncated:
                QTreeWidgetItem(
                    group,
                    [
                        "Field listing",
                        "Listing incomplete (safety cap or non-internal links)",
                    ],
                )
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
        dataset_groups: dict[str, QTreeWidgetItem] = {"": datasets}
        for dataset in summary.datasets:
            size_label = (
                f"{dataset.size_bytes:,} bytes"
                if dataset.size_bytes is not None
                else "variable/unknown size"
            )
            chunk_label = (
                f"chunks={dataset.chunks}"
                if dataset.chunks is not None
                else "contiguous"
            )
            allocated_label = (
                f"allocated={dataset.storage_bytes:,} bytes"
                if dataset.storage_bytes is not None
                else "allocated=unknown"
            )
            ratio = dataset.storage_ratio
            storage_label = (
                f"storage={ratio:.1%} of logical"
                if ratio is not None
                else "storage ratio=n/a"
            )
            compression_label = (
                f"compression={dataset.compression}"
                if dataset.compression is not None
                else "uncompressed"
            )
            details = (
                f"{dataset.shape} · {dataset.dtype} · {size_label} · "
                f"{allocated_label} · {storage_label} · "
                f"{chunk_label} · {compression_label}"
            )
            parts = dataset.path.strip("/").split("/")
            parent_path = ""
            parent = datasets
            for group_name in parts[:-1]:
                parent_path = f"{parent_path}/{group_name}"
                group_item = dataset_groups.get(parent_path)
                if group_item is None:
                    group_item = QTreeWidgetItem(parent, [group_name, "HDF5 group"])
                    dataset_groups[parent_path] = group_item
                parent = group_item
            item = QTreeWidgetItem(parent, [dataset.path, details])
            item.setToolTip(0, dataset.path)
            item.setToolTip(1, details)
            item.setData(0, Qt.ItemDataRole.UserRole, dataset.path)
        if summary.datasets_truncated:
            QTreeWidgetItem(
                datasets,
                ["Dataset listing", "Truncated at the configured metadata limit"],
            )
        if summary.tree_truncated:
            QTreeWidgetItem(
                datasets,
                ["HDF5 object tree", "Traversal stopped at the safety limit"],
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
        self.tree.expandToDepth(2)
        self._update_selection_actions()
