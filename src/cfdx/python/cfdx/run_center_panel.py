"""Qt solver cockpit with bounded logs and evidence-first convergence plots."""

from __future__ import annotations

import json
from pathlib import Path
from weakref import ref

from PySide6.QtCore import QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from .application import Application, ApplicationStateChanged
from .cell_hotspots_panel import CellHotspotsPanel
from .convergence_diagnostics import analyze_monitor_series, metric_group
from .theme import DEFAULT_THEME, style_scientific_plot


class RunCenterPanel(QWidget):
    """Live diagnostics, never a surrogate for the solver convergence gate."""

    _state_event = Signal(object)
    cell_frame_selected = Signal(object, object, str)
    full_cell_map_selected = Signal(object)
    cell_hotspot_selected = Signal(object, object)
    _GROUPS = (
        "All reported",
        "Nonlinear / equation",
        "Conservation",
        "Linear solve",
        "State updates",
        "Other reported",
    )

    def __init__(self, application: Application, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.application = application
        self.setObjectName("workbench.run_center")
        self.status = QLabel("No execution attached")
        self.metrics = QLabel("Iteration: — | Time: — | CFL: —")
        self.verdict = QLabel("Solver convergence evidence: not available")
        self.verdict.setObjectName("workbench.run.verdict")
        self.verdict.setWordWrap(True)
        self.output = QTextEdit()
        self.output.setObjectName("workbench.run.logs")
        self.output.setReadOnly(True)
        self.output.setPlaceholderText("Solver output will appear here")
        self._visible_output: tuple[str, ...] = ()
        self.run_button = QPushButton("Run")
        self.pause_button = QPushButton("Pause")
        self.resume_button = QPushButton("Resume")
        self.stop_button = QPushButton("Stop")
        self.spatial_diagnostics = QCheckBox("Cell hotspots (opt-in)", self)
        self.full_cell_map = QCheckBox("Full per-cell map (opt-in, 512 MiB cap)", self)
        self.full_cell_map.setObjectName("workbench.run.full_cell_map_opt_in")
        self.full_cell_map.setToolTip(
            "Observational full-field map, no physics acceptance change; max 512 MiB"
        )
        self.full_cell_map.toggled.connect(
            lambda enabled: (
                self.spatial_diagnostics.setChecked(True) if enabled else None
            )
        )
        self.spatial_diagnostics.setObjectName("workbench.run.cell_hotspots_opt_in")
        self.spatial_diagnostics.setToolTip(
            "Add --cell-hotspots to the native CFDX solver command for the next run. "
            "Creates a bounded, live per-cell diagnostic without changing solver tolerances."
        )
        controls = QHBoxLayout()
        for button in (
            self.run_button,
            self.pause_button,
            self.resume_button,
            self.stop_button,
        ):
            controls.addWidget(button)
        layout = QVBoxLayout(self)
        layout.addWidget(self.status)
        layout.addWidget(self.metrics)
        layout.addWidget(self.verdict)
        layout.addLayout(controls)
        layout.addWidget(self.spatial_diagnostics)
        layout.addWidget(self.full_cell_map)

        self.tabs = QTabWidget(self)
        self.tabs.setObjectName("workbench.run.tabs")
        convergence = QWidget(self)
        convergence_layout = QVBoxLayout(convergence)
        self.convergence_notice = QLabel(
            "Diagnostic only: slope/oscillation and plots are computed from actual "
            "solver telemetry; they are not scientific convergence criteria."
        )
        self.convergence_notice.setWordWrap(True)
        convergence_layout.addWidget(self.convergence_notice)
        self.metric_filter = QComboBox(self)
        self.metric_filter.setObjectName("workbench.run.metric_filter")
        self.metric_filter.addItems(self._GROUPS)
        self.metric_filter.currentTextChanged.connect(self._refresh_monitor_details)
        convergence_layout.addWidget(self.metric_filter)
        self.table = QTableWidget(0, 6, self)
        self.table.setObjectName("workbench.run.metric_table")
        self.table.setHorizontalHeaderLabels(
            (
                "Group",
                "Quantity",
                "Latest",
                "Trend",
                "Decades / iteration",
                "History",
            )
        )
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setStretchLastSection(True)
        convergence_layout.addWidget(self.table)
        self.plot = None
        self._plot_curves: dict[str, object] = {}
        try:
            import pyqtgraph as pg
        except ImportError:
            pg = None
        if pg is not None:
            plot = pg.PlotWidget(self)
            style_scientific_plot(plot, DEFAULT_THEME)
            plot.setObjectName("workbench.run.residual_plot")
            plot.setLogMode(x=False, y=True)
            plot.setLabel("bottom", "Iteration")
            plot.setLabel("left", "Reported magnitude (log scale)")
            plot.showGrid(x=True, y=True, alpha=0.18)
            plot.addLegend()
            convergence_layout.addWidget(plot)
            self.plot = plot
        self.tabs.addTab(convergence, "Convergence")
        # Native, persisted evidence is deliberately separate from transient
        # console telemetry. A filled Evidence tab requires *both* fresh
        # CFDX-EXECUTION and CFDX-CONVERGENCE v1 artifacts.
        evidence_tab = QWidget(self)
        evidence_layout = QVBoxLayout(evidence_tab)
        self.evidence_summary = QLabel(
            "No verified native execution artifacts for this run"
        )
        self.evidence_summary.setObjectName("workbench.run.evidence_summary")
        self.evidence_summary.setWordWrap(True)
        evidence_layout.addWidget(self.evidence_summary)
        self.evidence_table = QTableWidget(0, 3, self)
        self.evidence_table.setObjectName("workbench.run.evidence_table")
        self.evidence_table.setHorizontalHeaderLabels(
            ("Group", "Final metric", "Reported value")
        )
        self.evidence_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers
        )
        self.evidence_table.horizontalHeader().setStretchLastSection(True)
        evidence_layout.addWidget(self.evidence_table)
        self.export_evidence = QPushButton("Export verified evidence (JSON)")
        self.export_evidence.setObjectName("workbench.run.export_evidence")
        self.export_evidence.setEnabled(False)
        self.export_evidence.clicked.connect(self._export_evidence)
        evidence_layout.addWidget(self.export_evidence)
        self.tabs.addTab(evidence_tab, "Evidence")
        self.cell_hotspots = CellHotspotsPanel(self)
        self.tabs.addTab(self.cell_hotspots, "Cell hotspots")
        self.cell_hotspots.frame_selected.connect(self.cell_frame_selected.emit)
        self.cell_hotspots.full_map_selected.connect(self.full_cell_map_selected.emit)
        self.cell_hotspots.cell_selected.connect(self.cell_hotspot_selected.emit)
        self.tabs.addTab(self.output, "Logs")
        layout.addWidget(self.tabs)
        self._trends = ()
        self._last_series = None
        self._last_evidence = None
        self._last_evidence_error = None
        self.run_button.clicked.connect(self._run_with_cell_diagnostics)
        self.pause_button.clicked.connect(lambda: self._invoke("pause"))
        self.resume_button.clicked.connect(lambda: self._invoke("resume"))
        self.stop_button.clicked.connect(lambda: self._invoke("stop"))
        self._state_event.connect(self._state_changed)
        self._refresh_timer = QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.setInterval(120)
        self._refresh_timer.timeout.connect(self._flush_state)
        self._pending_state = None
        panel_ref = ref(self)

        def queue_state_event(event: object) -> None:
            panel = panel_ref()
            if panel is not None:
                panel._state_event.emit(event)

        application.events.subscribe(ApplicationStateChanged, queue_state_event)
        event_bus = application.events
        self.destroyed.connect(
            lambda *_args: event_bus.unsubscribe(
                ApplicationStateChanged, queue_state_event
            )
        )
        self._refresh(application.state)

    def apply_theme(self, theme: str) -> None:
        """Refresh only the plot chrome; solver records remain unchanged."""
        if self.plot is not None:
            style_scientific_plot(self.plot, theme)

    def _run_with_cell_diagnostics(self) -> None:
        """Never enable an unsupported telemetry flag on an unrelated runner."""
        center = self.application.run_center
        if center is None:
            self.status.setText("No execution controller attached")
            return
        runner = center.controller.runner
        command = tuple(runner.command)
        if self.spatial_diagnostics.isChecked() or self.full_cell_map.isChecked():
            executable = Path(command[0]).name.casefold()
            if executable not in (
                "cfdx_production_solver",
                "cfdx_production_solver.exe",
            ):
                self.status.setText(
                    "Cell hotspot opt-in requires the native cfdx_production_solver CLI. "
                    "Configure that executable or turn off the checkbox."
                )
                return
            extra = tuple(
                flag
                for flag in ("--cell-hotspots", "--cell-map-full")
                if flag not in command
                and (flag != "--cell-map-full" or self.full_cell_map.isChecked())
            )
            if extra:
                runner.command = command + extra
                self._spatial_flag_added = True
            elif (
                not self.full_cell_map.isChecked()
                and "--cell-map-full" in command
                and self._spatial_flag_added
            ):
                runner.command = tuple(
                    arg for arg in command if arg != "--cell-map-full"
                )
        elif getattr(self, "_spatial_flag_added", False):
            runner.command = tuple(
                arg
                for arg in command
                if arg not in {"--cell-hotspots", "--cell-map-full"}
            )
            self._spatial_flag_added = False
        self._invoke("run")

    def _invoke(self, action: str) -> None:
        try:
            getattr(self.application, action)()
        except (ValueError, RuntimeError, OSError, NotImplementedError) as exc:
            self.status.setText(f"{action.title()} unavailable: {exc}")

    def _state_changed(self, event) -> None:
        # Batch high-frequency stdout/stderr events; never update Qt widgets
        # directly from a solver thread.
        self._pending_state = event.state
        if not self._refresh_timer.isActive():
            self._refresh_timer.start()

    def _flush_state(self) -> None:
        if self._pending_state is not None:
            self._refresh(self._pending_state)
            self._pending_state = None

    def _refresh(self, _state) -> None:
        center = self.application.run_center
        if center is None:
            self.status.setText("No execution controller attached")
            self.metrics.setText("Iteration: — | Time: — | CFL: —")
            self.verdict.setText("Solver convergence evidence: not available")
            self.output.clear()
            self._visible_output = ()
            self.table.setRowCount(0)
            self._trends = ()
            self._last_series = None
            self._refresh_evidence(None, "no execution controller attached")
            # A different CFDX project must never inherit old spatial samples.
            self.cell_hotspots.set_source(None)
            for button in (
                self.run_button,
                self.pause_button,
                self.resume_button,
                self.stop_button,
            ):
                button.setEnabled(False)
            return
        state = center.state
        self.status.setText(f"Process state: {state.simulation_state.value}")
        cfl = state.latest_metrics.cfl if state.latest_metrics else None
        cfl_text = f"{cfl:g}" if cfl is not None else "—"
        self.metrics.setText(
            f"Iteration: {state.iteration} | Time: {state.time:g} | CFL: {cfl_text}"
        )
        if state.convergence_verdict is True:
            verdict_text = "Native execution verdict: converged (two fresh artifacts). Independent physical qualification remains required."
        elif state.convergence_verdict is False:
            verdict_text = "Native execution verdict: NOT converged (overall model including requested physics)."
        else:
            verdict_text = "NOT VERIFIED — no fresh compatible execution.json + convergence.json. A zero process exit is not convergence."
        self.verdict.setText(verdict_text)
        self._refresh_evidence(
            state.convergence_evidence, state.convergence_evidence_error
        )
        # Link the spatial trace only to runs launched with the native opt-in.
        # Never auto-load an old file from a directory used by another run.
        controller = center.controller
        command = tuple(controller.runner.command)
        output_dir = getattr(controller, "_run_output_dir", None)
        if "--cell-hotspots" in command and output_dir is not None:
            self.cell_hotspots.set_source(
                Path(output_dir) / "cell_convergence.jsonl", automatic=True
            )
        elif self.cell_hotspots.source is not None and getattr(
            self, "_spatial_flag_added", False
        ):
            self.cell_hotspots.set_source(None)

        output = state.output
        scroll = self.output.verticalScrollBar()
        follow = scroll.value() >= scroll.maximum() - 10
        if output[: len(self._visible_output)] == self._visible_output:
            for line in output[len(self._visible_output) :]:
                self.output.append(line)
        else:
            # A new run or bounded rolling history invalidates the prior view.
            self.output.setPlainText("\n".join(output))
        self._visible_output = output
        if follow:
            scroll.setValue(scroll.maximum())

        series = getattr(center.controller, "monitor_series", None)
        self._last_series = series
        try:
            self._trends = analyze_monitor_series(series)
        except ValueError as exc:
            self._trends = ()
            self.convergence_notice.setText(f"Invalid solver telemetry: {exc}")
        else:
            self.convergence_notice.setText(
                "Reported solver measurements only. Trend is a diagnostic, not a convergence certificate."
            )
        self._refresh_monitor_details()
        running = state.simulation_state.value == "RUNNING"
        paused = state.simulation_state.value == "PAUSED"
        self.run_button.setEnabled(not running and not paused)
        supports_pause = bool(
            getattr(center.controller.runner, "supports_pause", False)
        )
        self.pause_button.setEnabled(running and supports_pause)
        self.resume_button.setEnabled(paused and supports_pause)
        self.stop_button.setEnabled(running or paused)

    def _refresh_evidence(self, evidence, error: str | None) -> None:
        if evidence is self._last_evidence and error == self._last_evidence_error:
            return
        self._last_evidence = evidence
        self._last_evidence_error = error
        self.evidence_table.setRowCount(0)
        self.export_evidence.setEnabled(evidence is not None)
        if evidence is None:
            message = "No fresh, native execution certificate available."
            if error:
                message += f" Details: {error}"
            self.evidence_summary.setText(message)
            return
        self.evidence_summary.setText(
            f"Native solver verdict: {'CONVERGED' if evidence.overall_converged else 'NOT CONVERGED'}"
            f" | flow={evidence.flow_converged} | iterations={evidence.iterations}"
            f" | flow history={evidence.flow_history_count} | reason: {evidence.reason}\n"
            "SHA-256 integrity recorded for execution.json and convergence.json. "
            "This is a solver verdict, NOT independent verification/validation."
        )
        self.evidence_table.setRowCount(len(evidence.final_metrics))
        for row, (name, value) in enumerate(evidence.final_metrics):
            for col, label in enumerate(
                (
                    metric_group(name),
                    name,
                    f"{value:.9g}",
                )
            ):
                self.evidence_table.setItem(row, col, QTableWidgetItem(label))

    def _export_evidence(self) -> None:
        evidence = self._last_evidence
        if evidence is None:
            return
        path, _ = QFileDialog.getSaveFileName(
            self,
            "Export verified run evidence",
            "cfdx-run-evidence.json",
            "JSON (*.json)",
        )
        if not path:
            return
        try:
            Path(path).write_text(
                json.dumps(
                    evidence.to_dict(), indent=2, sort_keys=True, allow_nan=False
                )
                + "\n",
                encoding="utf-8",
            )
        except (OSError, TypeError, ValueError) as exc:
            self.evidence_summary.setText(f"Evidence export failed: {exc}")

    def _refresh_monitor_details(self, *_args) -> None:
        group = self.metric_filter.currentText()
        trends = [
            item
            for item in self._trends
            if group == "All reported" or item.group == group
        ]
        self.table.setSortingEnabled(False)
        self.table.setRowCount(len(trends))
        for row, trend in enumerate(trends):
            columns = (
                trend.group,
                trend.name,
                f"{trend.last:.5g}",
                trend.trend,
                f"{trend.decades_per_iteration:+.3g}"
                if trend.decades_per_iteration is not None
                else "—",
                trend.sparkline,
            )
            for col, value in enumerate(columns):
                self.table.setItem(row, col, QTableWidgetItem(value))
        if self.plot is not None:
            for curve in self._plot_curves.values():
                self.plot.removeItem(curve)
            self._plot_curves = {}
            series = self._last_series
            if series is not None:
                import pyqtgraph as pg

                for index, trend in enumerate(trends[:8]):
                    samples = [
                        sample
                        for sample in series.samples[-512:]
                        if trend.name in sample.values
                    ]
                    x = [sample.iteration for sample in samples]
                    y = [
                        max(abs(sample.values[trend.name]), 1e-300)
                        for sample in samples
                    ]
                    curve = self.plot.plot(
                        x,
                        y,
                        name=trend.name,
                        pen=pg.mkPen(pg.intColor(index, hues=8), width=1.5),
                    )
                    self._plot_curves[trend.name] = curve
        self.table.setSortingEnabled(True)
