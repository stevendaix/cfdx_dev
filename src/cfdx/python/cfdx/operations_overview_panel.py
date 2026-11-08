"""Connected, read-only CFDX Mission Control panel.

Only the application/event layer owns scientific state. The panel never
synthesizes solver data, launches arbitrary files or certifies a CFD result.
"""
from __future__ import annotations

import json
import os
from pathlib import Path
import tempfile

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import (
    QAbstractItemView, QFileDialog, QFrame, QGridLayout, QHBoxLayout,
    QLabel, QListWidget, QListWidgetItem, QPushButton, QVBoxLayout, QWidget,
)

from .application.state import ApplicationState
from .operations_overview import OperationsOverview, build_operations_overview
from .theme import DEFAULT_THEME, theme_tokens


class OperationsOverviewPanel(QWidget):
    """Project snapshot, check provenance and verified native-solver evidence."""

    workspace_requested = Signal(str)
    check_requested = Signal()
    open_mesh_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("workbench.operations_overview")
        self.snapshot: OperationsOverview | None = None
        self._theme = DEFAULT_THEME

        heading = QLabel("MISSION CONTROL", self)
        heading.setObjectName("workbench.overview.heading")
        heading.setProperty("role", "eyebrow")
        subtitle = QLabel("State from CFDX • configuration is not convergence", self)
        subtitle.setProperty("role", "secondary")
        subtitle.setWordWrap(True)
        self.project_label = QLabel("No project snapshot", self)
        self.project_label.setObjectName("workbench.overview.project")
        self.project_label.setWordWrap(True)
        self.project_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        cards = QWidget(self)
        card_layout = QGridLayout(cards)
        card_layout.setSpacing(8)
        self.configuration_card, self.configuration_value = self._create_card(
            "Configuration", "—", "workbench.overview.configuration")
        self.checks_card, self.checks_value = self._create_card(
            "Setup checks", "Unverified", "workbench.overview.checks")
        self.evidence_card, self.evidence_value = self._create_card(
            "Solver evidence", "No certificate", "workbench.overview.evidence")
        for index, card in enumerate((self.configuration_card, self.checks_card, self.evidence_card)):
            card_layout.addWidget(card, 0, index)
        self.detail = QLabel("No case evidence loaded", self)
        self.detail.setObjectName("workbench.overview.detail")
        self.detail.setProperty("role", "secondary")
        self.detail.setWordWrap(True)
        self.finding_title = QLabel("Next actions and diagnostics", self)
        self.findings = QListWidget(self)
        self.findings.setObjectName("workbench.overview.findings")
        self.findings.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.findings.itemDoubleClicked.connect(self._navigate_item)
        self.findings.setMinimumHeight(145)
        self.info_label = QLabel(
            "Workflow COMPLETE means configuration recorded, not a scientifically "
            "qualified simulation. Physical validation remains independent.", self
        )
        self.info_label.setWordWrap(True)
        self.info_label.setProperty("role", "secondary")

        self.check_button = QPushButton("Check case", self)
        self.check_button.setObjectName("workbench.overview.check")
        self.check_button.setToolTip("Run the current CFDX setup validator; this does not execute the solver")
        self.check_button.clicked.connect(lambda _checked=False: self.check_requested.emit())
        self.mesh_button = QPushButton("Inspect mesh", self)
        self.mesh_button.setObjectName("workbench.overview.mesh")
        self.mesh_button.clicked.connect(lambda _checked=False: self.open_mesh_requested.emit())
        self.navigate_button = QPushButton("Open selected", self)
        self.navigate_button.setObjectName("workbench.overview.navigate")
        self.navigate_button.clicked.connect(self._navigate_selected)
        self.export_button = QPushButton("Export snapshot", self)
        self.export_button.setObjectName("workbench.overview.export")
        self.export_button.clicked.connect(self._export_snapshot)
        self.export_button.setEnabled(False)
        actions = QHBoxLayout()
        for button in (self.check_button, self.mesh_button, self.navigate_button, self.export_button):
            actions.addWidget(button)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(9)
        layout.addWidget(heading)
        layout.addWidget(subtitle)
        layout.addWidget(self.project_label)
        layout.addWidget(cards)
        layout.addWidget(self.detail)
        layout.addWidget(self.finding_title)
        layout.addWidget(self.findings, 1)
        layout.addLayout(actions)
        layout.addWidget(self.info_label)
        self.apply_theme(DEFAULT_THEME)

    def _create_card(self, title: str, initial: str, object_name: str) -> tuple[QFrame, QLabel]:
        frame = QFrame(self)
        frame.setObjectName(object_name)
        frame.setProperty("role", "dashboardCard")
        frame_layout = QVBoxLayout(frame)
        label = QLabel(title, frame)
        label.setProperty("role", "secondary")
        value = QLabel(initial, frame)
        value.setObjectName(f"{object_name}.value")
        value.setProperty("role", "dashboardValue")
        value.setWordWrap(True)
        frame_layout.addWidget(label)
        frame_layout.addWidget(value)
        return frame, value

    def apply_theme(self, theme: str) -> None:
        """Refresh diagnostic text highlights for a theme switch."""
        self._theme = theme
        if self.snapshot is not None:
            self._fill_findings(self.snapshot)

    def update_state(self, state: ApplicationState, *, run=None, mesh_catalog_loaded: bool = False) -> None:
        snapshot = build_operations_overview(
            state, run=run, mesh_catalog_loaded=mesh_catalog_loaded
        )
        self.snapshot = snapshot
        self.project_label.setText(
            snapshot.project + ("  •  modified" if snapshot.unsaved_changes else "")
        )
        self.project_label.setToolTip(snapshot.project_path or "Unsaved in-memory case")
        self.configuration_value.setText(f"{snapshot.configured_sections}/{snapshot.section_count} recorded")
        self.checks_value.setText({
            "passed": "Passed (setup)", "warnings": "Warnings", "errors": "Errors",
            "stale": "Check outdated", "not_checked": "Not checked",
        }[snapshot.check_state])
        self.evidence_value.setText({
            "converged": "Converged (solver)",
            "not_converged": "Not converged",
            "missing_or_unverified": "Not verified",
        }[snapshot.solver_evidence])
        self.detail.setText(
            f"{snapshot.execution_state}  •  iteration {snapshot.iteration:,}  •  "
            f"time {snapshot.physical_time:g}  •  {snapshot.results_frames} result frame(s)  •  "
            f"check scope: {snapshot.check_scope.replace('_', ' ')}"
        )
        self._fill_findings(snapshot)
        self.export_button.setEnabled(True)

    def _fill_findings(self, snapshot: OperationsOverview) -> None:
        previous = self.findings.currentItem()
        previous_key = previous.data(Qt.ItemDataRole.UserRole + 1) if previous else None
        self.findings.clear()
        tokens = theme_tokens(self._theme)
        palette = {
            "error": tokens["status.error"],
            "warning": tokens["status.warning"],
            "action": tokens["accent.primary"],
            "information": tokens["text.secondary"],
        }
        for finding in snapshot.findings:
            item = QListWidgetItem(f"{finding.title}  •  {finding.workspace.title()}")
            item.setToolTip(finding.detail)
            item.setData(Qt.ItemDataRole.UserRole, finding.workspace)
            item.setData(Qt.ItemDataRole.UserRole + 1, finding.key)
            item.setForeground(QBrush(QColor(palette[finding.level])))
            self.findings.addItem(item)
            if finding.key == previous_key:
                self.findings.setCurrentItem(item)
        if not snapshot.findings:
            item = QListWidgetItem("No pending overview actions. Consult independent V&V before qualification.")
            item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsSelectable)
            self.findings.addItem(item)
        elif self.findings.currentRow() < 0:
            self.findings.setCurrentRow(0)
        self.navigate_button.setEnabled(bool(snapshot.findings))

    def _navigate_item(self, item: QListWidgetItem) -> None:
        workspace = item.data(Qt.ItemDataRole.UserRole)
        if workspace in {"setup", "mesh", "run", "analyze", "validate"}:
            self.workspace_requested.emit(workspace)

    def _navigate_selected(self) -> None:
        if self.findings.currentItem() is not None:
            self._navigate_item(self.findings.currentItem())

    def _export_snapshot(self) -> None:
        if self.snapshot is None:
            return
        filename, _ = QFileDialog.getSaveFileName(
            self, "Export CFDX Mission Control", "cfdx_overview.json", "JSON (*.json)"
        )
        if not filename:
            return
        dest = Path(filename)
        temporary: Path | None = None
        try:
            payload = json.dumps(self.snapshot.to_dict(), ensure_ascii=False,
                                 sort_keys=True, indent=2, allow_nan=False) + "\n"
            # Use a sibling temporary file + atomic rename: partial reports
            # must never be mistaken for a valid evidence bundle.
            with tempfile.NamedTemporaryFile(
                mode="w", encoding="utf-8", prefix=".cfdx-overview-",
                suffix=".tmp", dir=dest.parent, delete=False
            ) as handle:
                temporary = Path(handle.name)
                handle.write(payload)
            os.replace(temporary, dest)
            self.info_label.setText(f"Snapshot exported: {dest.name} (not a physics qualification)")
        except (OSError, ValueError) as exc:
            self.info_label.setText(f"Could not export snapshot: {exc}")
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)
