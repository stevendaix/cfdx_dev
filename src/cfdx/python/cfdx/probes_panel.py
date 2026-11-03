"""Qt adapter for the application-level point-probe catalogue."""

from __future__ import annotations

from collections import deque
from pathlib import Path
from weakref import ref

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QFileDialog,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .application import Application, ApplicationStateChanged, ProbeSampleReceived
from .application.commands import ConfigureProbes
from .probe import Probe, ProbeCatalog, ProbeField
from .probe_plot import create_probe_plot, update_probe_plot
from .probe_validation import ProbeSample, ProbeSeries
from .session import SimulationState

_FIELD_CHOICES = [field.cli for field in ProbeField]
_MAX_LIVE_SAMPLES_PER_PROBE = 2000


class ProbesPanel(QWidget):
    """Editor for the case's point-probe catalogue.

    The catalogue is committed to the case through the application's
    :class:`~cfdx.application.commands.ConfigureProbes` command; the execution
    backend then forwards probes to the solver as ``--probe``/``--probe-csv``
    arguments, producing the probe CSV at run time.
    """

    _state_event = Signal(object)
    _probe_event = Signal(object)

    def __init__(self, application: Application, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.application = application
        self.setObjectName("workbench.probes_panel")

        self._table = QTableWidget(0, 5)
        self._table.setHorizontalHeaderLabels(["name", "x", "y", "z", "field"])
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self._table.itemChanged.connect(self._on_item_changed)

        self._add_button = QPushButton("Add")
        self._add_button.clicked.connect(self._add_row)
        self._remove_button = QPushButton("Remove")
        self._remove_button.clicked.connect(self._remove_selected)
        self._remove_button.setEnabled(False)
        self._table.selectionModel().selectionChanged.connect(
            self._on_selection_changed
        )

        self._apply_button = QPushButton("Apply to case")
        self._apply_button.setEnabled(False)
        self._apply_button.clicked.connect(self._apply)

        self._load_results_button = QPushButton("Load results CSV…")
        self._load_results_button.clicked.connect(self._load_results)
        self._export_results_button = QPushButton("Export results CSV…")
        self._export_results_button.setEnabled(False)
        self._export_results_button.clicked.connect(self._export_results)
        self._export_live_button = QPushButton("Export live stream…")
        self._export_live_button.setEnabled(False)
        self._export_live_button.setToolTip(
            "Export all samples retained in the live buffer, not only visible curves"
        )
        self._export_live_button.clicked.connect(self._export_live_stream)
        self._results_summary = QLabel("No probe results loaded")
        self._results_summary.setObjectName("workbench.probes.results_summary")
        self._results_summary.setProperty("role", "secondary")
        self._loaded_probe_series: tuple[ProbeSeries, ...] = ()
        self._display_probe_series: tuple[ProbeSeries, ...] = ()
        self._live_samples: dict[str, deque[ProbeSample]] = {}
        self._live_metadata: dict[str, tuple[tuple[float, ...], str, str]] = {}
        self._live_run_active = False
        self._live_dirty = False
        self._rendered_probe_catalog: tuple[Probe, ...] = ()
        self._plot_field_selector = QComboBox()
        self._plot_field_selector.setObjectName("workbench.probes.plot_field")
        self._plot_field_selector.setEnabled(False)
        self._plot_axis_selector = QComboBox()
        self._plot_axis_selector.setObjectName("workbench.probes.plot_axis")
        self._plot_axis_selector.addItem("Physical time", "time")
        self._plot_axis_selector.addItem("Iteration", "iteration")
        self._plot_axis_selector.currentIndexChanged.connect(self._redraw_probe_plot)
        self._plot_probe_list = QListWidget()
        self._plot_probe_list.setObjectName("workbench.probes.plot_series")
        self._plot_probe_list.setMaximumHeight(100)
        self._plot_probe_list.itemChanged.connect(self._redraw_probe_plot)
        try:
            self._probe_plot = create_probe_plot()
        except RuntimeError as exc:
            self._probe_plot = None
            self._probe_plot_unavailable = QLabel(str(exc))
            self._probe_plot_unavailable.setWordWrap(True)
        else:
            self._probe_plot_unavailable = None

        self._dirty = False

        controls = QHBoxLayout()
        for button in (self._add_button, self._remove_button, self._apply_button):
            controls.addWidget(button)

        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        results_controls = QHBoxLayout()
        results_controls.addWidget(self._load_results_button)
        results_controls.addWidget(self._export_results_button)
        results_controls.addWidget(self._export_live_button)
        results_controls.addWidget(self._results_summary, 1)
        layout.addLayout(results_controls)
        plot_controls = QHBoxLayout()
        plot_controls.addWidget(QLabel("Field / unit"))
        plot_controls.addWidget(self._plot_field_selector, 1)
        plot_controls.addWidget(QLabel("X axis"))
        plot_controls.addWidget(self._plot_axis_selector)
        self._plot_group = QGroupBox("Interactive probe history")
        plot_layout = QVBoxLayout(self._plot_group)
        plot_layout.addLayout(plot_controls)
        plot_layout.addWidget(self._plot_probe_list)
        if self._probe_plot is not None:
            plot_layout.addWidget(self._probe_plot, 1)
        else:
            plot_layout.addWidget(self._probe_plot_unavailable)
        layout.addWidget(self._plot_group, 2)
        layout.addWidget(self._table)

        self._plot_field_selector.currentIndexChanged.connect(
            self._populate_probe_choices
        )

        self._state_event.connect(self._refresh)
        self._probe_event.connect(self._on_probe_sample)
        panel_ref = ref(self)

        def queue_state_event(event: object) -> None:
            panel = panel_ref()
            if panel is not None:
                panel._state_event.emit(event)

        def queue_probe_event(event: object) -> None:
            panel = panel_ref()
            if panel is not None:
                panel._probe_event.emit(event)

        self.application.events.subscribe(ApplicationStateChanged, queue_state_event)
        self.application.events.subscribe(ProbeSampleReceived, queue_probe_event)
        event_bus = self.application.events
        self.destroyed.connect(
            lambda *_args: event_bus.unsubscribe(
                ApplicationStateChanged, queue_state_event
            )
        )
        self.destroyed.connect(
            lambda *_args: event_bus.unsubscribe(ProbeSampleReceived, queue_probe_event)
        )
        self._live_timer = QTimer(self)
        self._live_timer.setInterval(100)
        self._live_timer.timeout.connect(self._flush_live_plot)
        self._live_timer.start()
        self._refresh()

    def _refresh(self, *_args: object) -> None:
        """Rebuild the table from the case's current probe catalogue."""
        catalog = tuple(self.application.session.case.probes)
        if catalog == self._rendered_probe_catalog:
            self._update_live_run_state()
            return
        self._table.blockSignals(True)
        try:
            self._table.setRowCount(0)
            for probe in catalog:
                self._add_row(probe)
        finally:
            self._table.blockSignals(False)
        self._rendered_probe_catalog = catalog
        self._set_dirty(False)
        self._update_live_run_state()

    def _update_live_run_state(self) -> None:
        running = self.application.session.state in {
            SimulationState.RUNNING,
            SimulationState.PAUSED,
        }
        if running and not self._live_run_active:
            self._begin_live_run()
        elif not running and self._live_run_active:
            self._live_run_active = False

    def _begin_live_run(self) -> None:
        self._live_run_active = True
        self._live_samples.clear()
        self._live_metadata.clear()
        self._display_probe_series = ()
        self._live_dirty = False
        self._export_live_button.setEnabled(False)
        self._results_summary.setText("Waiting for live probe samples…")
        iteration_index = self._plot_axis_selector.findData("iteration")
        if iteration_index >= 0:
            self._plot_axis_selector.setCurrentIndex(iteration_index)
        self._populate_probe_plot_controls()

    def _on_probe_sample(self, event: ProbeSampleReceived) -> None:
        sample = event.sample
        if not self._live_run_active:
            self._begin_live_run()
        definition = next(
            (
                probe
                for probe in self.application.session.case.probes
                if probe.name == sample.name
            ),
            None,
        )
        if definition is not None and (
            definition.field.cli != sample.field or definition.field.unit != sample.unit
        ):
            self._show_results_message(
                f"Live probe metadata mismatch for {sample.name!r}"
            )
            return
        point = (
            (definition.x, definition.y, definition.z) if definition is not None else ()
        )
        metadata = (point, sample.field, sample.unit)
        previous_metadata = self._live_metadata.get(sample.name)
        if previous_metadata is not None and previous_metadata != metadata:
            self._show_results_message(
                f"Live probe metadata changed during run for {sample.name!r}"
            )
            return
        history = self._live_samples.get(sample.name)
        if history is None:
            history = deque(maxlen=_MAX_LIVE_SAMPLES_PER_PROBE)
            self._live_samples[sample.name] = history
            self._live_metadata[sample.name] = metadata
        if history and (
            sample.iteration <= history[-1].iteration or sample.time < history[-1].time
        ):
            self._show_results_message(
                f"Out-of-order live sample for probe {sample.name!r} was ignored"
            )
            return
        history.append(ProbeSample(sample.iteration, sample.time, sample.value))
        self._live_dirty = True
        self._export_live_button.setEnabled(True)

    def _flush_live_plot(self) -> None:
        if not self._live_dirty:
            return
        series: list[ProbeSeries] = []
        for name in sorted(self._live_samples):
            samples = self._live_samples[name]
            point, field, unit = self._live_metadata[name]
            series.append(
                ProbeSeries(name, point, tuple(samples), field=field, unit=unit)
            )
        self._display_probe_series = tuple(series)
        sample_count = sum(len(probe.samples) for probe in series)
        iterations = [sample.iteration for probe in series for sample in probe.samples]
        interval = f"iterations {min(iterations)}–{max(iterations)}"
        probe_label = "probe" if len(series) == 1 else "probes"
        self._results_summary.setText(
            f"Live · {len(series)} {probe_label} · {sample_count} buffered samples · {interval}"
        )
        self._populate_probe_plot_controls(preserve_selection=True)
        self._live_dirty = False

    def _add_row(self, probe: Probe | None = None) -> None:
        row = self._table.rowCount()
        self._table.insertRow(row)
        self._table.setItem(row, 0, QTableWidgetItem(probe.name if probe else ""))
        self._table.setItem(row, 1, QTableWidgetItem(repr(probe.x) if probe else "0"))
        self._table.setItem(row, 2, QTableWidgetItem(repr(probe.y) if probe else "0"))
        self._table.setItem(row, 3, QTableWidgetItem(repr(probe.z) if probe else "0"))
        combo = QComboBox()
        combo.addItems(_FIELD_CHOICES)
        combo.setCurrentText(probe.field.cli if probe else _FIELD_CHOICES[0])
        self._table.setCellWidget(row, 4, combo)

    def _remove_selected(self) -> None:
        for index in sorted(
            {index.row() for index in self._table.selectionModel().selectedRows()},
            reverse=True,
        ):
            self._table.removeRow(index)
        self._set_dirty(True)

    def _on_item_changed(self, *_args: object) -> None:
        self._set_dirty(True)

    def _on_selection_changed(self) -> None:
        self._remove_button.setEnabled(
            bool(self._table.selectionModel().selectedRows())
        )

    def _set_dirty(self, dirty: bool) -> None:
        self._dirty = dirty
        self._apply_button.setEnabled(dirty)

    def _apply(self) -> None:
        catalog = ProbeCatalog()
        for row in range(self._table.rowCount()):
            name = (self._table.item(row, 0).text() or "").strip()
            try:
                x = float(self._table.item(row, 1).text())
                y = float(self._table.item(row, 2).text())
                z = float(self._table.item(row, 3).text())
            except (TypeError, ValueError, AttributeError):
                continue
            combo = self._table.cellWidget(row, 4)
            field_token = combo.currentText() if combo is not None else ""
            if not name or field_token not in _FIELD_CHOICES:
                continue
            catalog.add(
                Probe(name=name, x=x, y=y, z=z, field=ProbeField.from_cli(field_token))
            )
        try:
            self.application.execute(ConfigureProbes(catalog))
        except ValueError as exc:
            parent = self.window()
            if parent is not None and hasattr(parent, "statusBar"):
                parent.statusBar().showMessage(f"Probes: {exc}")
        self._set_dirty(False)

    def _load_results(self) -> None:
        """Load a solver probe CSV through the shared application contract."""
        source, _selected_filter = QFileDialog.getOpenFileName(
            self,
            "Load probe results",
            "",
            "CFDX probe CSV (*.csv);;All files (*)",
        )
        if not source:
            return
        self._loaded_probe_series = ()
        self._display_probe_series = ()
        self._export_results_button.setEnabled(False)
        self._results_summary.setText("No probe results loaded")
        self._populate_probe_plot_controls()
        try:
            series = self.application.read_probe_results(source)
        except (OSError, ValueError) as exc:
            self._show_results_message(f"Could not load probe results: {exc}")
            return
        self._loaded_probe_series = series
        self._display_probe_series = series
        self._populate_probe_plot_controls()
        sample_count = sum(len(probe.samples) for probe in series)
        iterations = [sample.iteration for probe in series for sample in probe.samples]
        interval = (
            f"iterations {min(iterations)}–{max(iterations)}"
            if iterations
            else "no samples"
        )
        probe_label = "probe" if len(series) == 1 else "probes"
        self._results_summary.setText(
            f"{len(series)} {probe_label} · {sample_count} samples · {interval}"
        )
        self._export_results_button.setEnabled(bool(series))
        self._show_results_message(f"Loaded probe results: {source}")

    def _probe_group(self, probe: ProbeSeries) -> tuple[str, str]:
        return (probe.field or "value", probe.unit or "")

    def _populate_probe_plot_controls(
        self, *, preserve_selection: bool = False
    ) -> None:
        group_data = self._plot_field_selector.currentData()
        current_group = tuple(group_data) if group_data is not None else None
        selected_names = (
            {
                self._plot_probe_list.item(row).data(Qt.ItemDataRole.UserRole)
                for row in range(self._plot_probe_list.count())
                if self._plot_probe_list.item(row).checkState() == Qt.CheckState.Checked
            }
            if preserve_selection and current_group is not None
            else None
        )
        groups = sorted(
            {self._probe_group(probe) for probe in self._display_probe_series}
        )
        self._plot_field_selector.blockSignals(True)
        self._plot_field_selector.clear()
        for field, unit in groups:
            label = f"{field} ({unit})" if unit else field
            self._plot_field_selector.addItem(label, (field, unit))
        self._plot_field_selector.setEnabled(
            bool(groups) and self._probe_plot is not None
        )
        if preserve_selection and current_group in groups:
            matching_index = next(
                (
                    index
                    for index in range(self._plot_field_selector.count())
                    if tuple(self._plot_field_selector.itemData(index)) == current_group
                ),
                -1,
            )
            if matching_index >= 0:
                self._plot_field_selector.setCurrentIndex(matching_index)
        elif preserve_selection and current_group is not None:
            selected_names = None
        self._plot_field_selector.blockSignals(False)
        self._populate_probe_choices(
            selected_names=selected_names if preserve_selection else None
        )

    def _populate_probe_choices(
        self, *_args: object, selected_names: set[str] | None = None
    ) -> None:
        group = self._plot_field_selector.currentData()
        self._plot_probe_list.blockSignals(True)
        self._plot_probe_list.clear()
        if group is not None:
            for probe in self._display_probe_series:
                if self._probe_group(probe) != tuple(group):
                    continue
                item = QListWidgetItem(probe.name)
                item.setData(Qt.ItemDataRole.UserRole, probe.name)
                item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
                item.setCheckState(
                    Qt.CheckState.Checked
                    if selected_names is None or probe.name in selected_names
                    else Qt.CheckState.Unchecked
                )
                self._plot_probe_list.addItem(item)
        self._plot_probe_list.setEnabled(self._plot_probe_list.count() > 0)
        self._plot_probe_list.blockSignals(False)
        self._redraw_probe_plot()

    def _redraw_probe_plot(self, *_args: object) -> None:
        if self._probe_plot is None:
            return
        group = self._plot_field_selector.currentData()
        selected = {
            self._plot_probe_list.item(row).data(Qt.ItemDataRole.UserRole)
            for row in range(self._plot_probe_list.count())
            if self._plot_probe_list.item(row).checkState() == Qt.CheckState.Checked
        }
        visible = tuple(
            probe
            for probe in self._display_probe_series
            if group is not None
            and self._probe_group(probe) == tuple(group)
            and probe.name in selected
        )
        update_probe_plot(
            self._probe_plot,
            visible,
            x_axis=self._plot_axis_selector.currentData() or "time",
        )

    def _export_results(self) -> None:
        """Write the loaded history using the canonical native CSV schema."""
        if not self._loaded_probe_series:
            self._show_results_message("Load probe results before exporting")
            return
        destination, _selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export probe results",
            "cfdx-probes.csv",
            "CFDX probe CSV (*.csv);;All files (*)",
        )
        if not destination:
            return
        try:
            exported = self.application.export_probe_results(
                destination, self._loaded_probe_series
            )
        except (OSError, ValueError) as exc:
            self._show_results_message(f"Could not export probe results: {exc}")
            return
        self._show_results_message(f"Exported probe results: {exported}")

    def _live_probe_series(self) -> tuple[ProbeSeries, ...]:
        """Return every sample retained by the live buffers, independent of plot filtering."""
        return tuple(
            ProbeSeries(
                name,
                self._live_metadata[name][0],
                tuple(self._live_samples[name]),
                field=self._live_metadata[name][1],
                unit=self._live_metadata[name][2],
            )
            for name in sorted(self._live_samples)
            if self._live_samples[name]
        )

    def _export_live_stream(self) -> None:
        """Export the bounded live-stream snapshot as canonical CSV or versioned JSON."""
        series = self._live_probe_series()
        if not series:
            self._show_results_message("No live probe samples are available to export")
            self._export_live_button.setEnabled(False)
            return

        destination, selected_filter = QFileDialog.getSaveFileName(
            self,
            "Export live probe stream",
            "cfdx-live-probes.csv",
            "CFDX probe CSV (*.csv);;CFDX probe JSON (*.json)",
        )
        if not destination:
            return
        export_json = "JSON" in selected_filter.upper()
        target = Path(destination)
        expected_suffix = ".json" if export_json else ".csv"
        if target.suffix.lower() != expected_suffix:
            target = target.with_suffix(expected_suffix)

        try:
            if export_json:
                exported = self.application.export_probe_results_json(
                    target,
                    series,
                    max_samples_per_probe=_MAX_LIVE_SAMPLES_PER_PROBE,
                )
            else:
                exported = self.application.export_probe_results(target, series)
        except (OSError, ValueError) as exc:
            self._show_results_message(f"Could not export live probe stream: {exc}")
            return
        sample_count = sum(len(probe.samples) for probe in series)
        self._show_results_message(
            f"Exported {sample_count} buffered live probe samples: {exported}"
        )

    def _show_results_message(self, message: str) -> None:
        parent = self.window()
        if parent is not None and hasattr(parent, "statusBar"):
            parent.statusBar().showMessage(message)
