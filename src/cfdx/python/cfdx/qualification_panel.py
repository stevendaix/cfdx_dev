"""Qt view of the read-only, commit-unverified CFDX qualification register."""

from __future__ import annotations

import json
import re
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPlainTextEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .qualification_center import (
    STATUSES,
    QualificationCase,
    QualificationRegistry,
    discover_qualification_registry,
    read_qualification_registry,
)


class QualificationCenterPanel(QWidget):
    """View many qualification cases without claiming new HEAD evidence."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("workbench.qualification.panel")
        self.registry: QualificationRegistry | None = None
        self._visible_cases: dict[str, QualificationCase] = {}
        self._proofs: dict[str, str] = {}
        self._proof_commit: str | None = None

        layout = QVBoxLayout(self)
        self.caution = QLabel(
            "Registry declarations only — PASS does not certify this checkout. "
            "Current-commit evidence must be checked independently."
        )
        self.caution.setObjectName("workbench.qualification.provenance")
        self.caution.setWordWrap(True)
        layout.addWidget(self.caution)
        self.path_label = QLabel("No qualification registry loaded")
        self.path_label.setObjectName("workbench.qualification.path")
        self.path_label.setWordWrap(True)
        layout.addWidget(self.path_label)

        actions = QHBoxLayout()
        self.open_button = QPushButton("Open registry…")
        self.open_button.setObjectName("workbench.qualification.open")
        self.open_button.clicked.connect(self.open_registry_dialog)
        actions.addWidget(self.open_button)
        self.reload_button = QPushButton("Reload")
        self.reload_button.setObjectName("workbench.qualification.reload")
        self.reload_button.setEnabled(False)
        self.reload_button.clicked.connect(self.reload)
        actions.addWidget(self.reload_button)
        self.proof_button = QPushButton("Open CI evidence manifest…")
        self.proof_button.setObjectName("workbench.qualification.evidence")
        self.proof_button.clicked.connect(self.open_evidence_dialog)
        actions.addWidget(self.proof_button)
        layout.addLayout(actions)

        search_row = QHBoxLayout()
        self.search = QLineEdit()
        self.search.setObjectName("workbench.qualification.search")
        self.search.setPlaceholderText("Find case, domain or ID…")
        self.search.textChanged.connect(self.refresh)
        search_row.addWidget(self.search)
        self.status_filter = QComboBox()
        self.status_filter.setObjectName("workbench.qualification.status_filter")
        self.status_filter.addItem("All statuses", "")
        for status in STATUSES:
            self.status_filter.addItem(status, status)
        self.status_filter.currentIndexChanged.connect(self.refresh)
        search_row.addWidget(self.status_filter)
        self.domain_filter = QComboBox()
        self.domain_filter.setObjectName("workbench.qualification.domain_filter")
        self.domain_filter.addItem("All domains", "")
        self.domain_filter.currentIndexChanged.connect(self.refresh)
        search_row.addWidget(self.domain_filter)
        layout.addLayout(search_row)

        self.summary = QLabel("0 cases — no registry loaded")
        self.summary.setObjectName("workbench.qualification.summary")
        layout.addWidget(self.summary)
        self.table = QTableWidget(0, 5, self)
        self.table.setObjectName("workbench.qualification.table")
        self.table.setHorizontalHeaderLabels(
            ("Case ID", "Case", "Domain", "Declared status", "CI proof")
        )
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.horizontalHeader().setStretchLastSection(True)
        self.table.setSortingEnabled(True)
        self.table.itemSelectionChanged.connect(self._show_details)
        layout.addWidget(self.table, 3)
        self.details = QPlainTextEdit()
        self.details.setObjectName("workbench.qualification.details")
        self.details.setReadOnly(True)
        self.details.setPlaceholderText(
            "Select a case to inspect its recorded criteria."
        )
        layout.addWidget(self.details, 2)
        self.message = QLabel("")
        self.message.setObjectName("workbench.qualification.message")
        self.message.setWordWrap(True)
        layout.addWidget(self.message)

        candidate = discover_qualification_registry()
        if candidate is not None:
            try:
                self.load_registry(candidate)
            except (OSError, UnicodeError, ValueError, RuntimeError) as exc:
                self.message.setText(f"Registry could not be loaded: {exc}")

    def load_registry(self, path: str | Path) -> None:
        """Parse before replacing displayed state; failure is transactional."""
        new_registry = read_qualification_registry(path)
        self.registry = new_registry
        self.path_label.setText(f"Source: {new_registry.source}")
        self.reload_button.setEnabled(True)
        self.message.clear()
        self.domain_filter.blockSignals(True)
        self.domain_filter.clear()
        self.domain_filter.addItem("All domains", "")
        for domain in sorted({case.domain for case in new_registry.cases}):
            self.domain_filter.addItem(domain, domain)
        self.domain_filter.blockSignals(False)
        self.refresh()

    def open_evidence_dialog(self) -> None:
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open commit-bound CFDX evidence manifest", "", "JSON (*.json)"
        )
        if filename:
            try:
                self.load_evidence_manifest(filename)
            except (
                OSError,
                ValueError,
                TypeError,
                KeyError,
                json.JSONDecodeError,
            ) as exc:
                self.message.setText(f"CI evidence rejected: {exc}")

    def load_evidence_manifest(self, source: str | Path) -> None:
        if self.registry is None:
            raise ValueError("load the case registry before evidence")
        source = Path(source)
        if (
            source.is_symlink()
            or not source.is_file()
            or source.stat().st_size > 2 * 1024 * 1024
        ):
            raise ValueError("unsafe or oversized CI manifest")
        data = json.loads(source.read_text(encoding="utf-8"))
        sha = data.get("commit_sha")
        if (
            data.get("format") != "CFDX-CI-QUALIFICATION-EVIDENCE"
            or data.get("schema_version") != 1
            or not isinstance(sha, str)
            or not re.fullmatch(r"[0-9a-f]{40}", sha)
        ):
            raise ValueError("invalid CI manifest schema or commit")
        entries = data.get("cases")
        if not isinstance(entries, list):
            raise ValueError("missing case evidence list")
        proofs = {}
        allowed = {case.id for case in self.registry.cases}
        for item in entries:
            if not isinstance(item, dict) or item.get("case_id") not in allowed:
                raise ValueError("unknown qualification case")
            key = item["case_id"]
            if key in proofs or item.get("proof_status") not in (
                "PASS_EVIDENCED",
                "UNVERIFIED",
            ):
                raise ValueError("duplicate or invalid case proof")
            proofs[key] = item["proof_status"]
        self._proofs = proofs
        self._proof_commit = sha
        self.message.setText(
            f"External CI manifest: {sha[:12]} — not proof for this checkout unless the HEAD matches and original artifacts are trusted"
        )
        self.refresh()

    def open_registry_dialog(self) -> None:
        filename, _ = QFileDialog.getOpenFileName(
            self, "Open CFDX Qualification Registry", "", "JSON files (*.json)"
        )
        if filename:
            self._try_load(filename)

    def reload(self) -> None:
        if self.registry is not None:
            self._try_load(self.registry.source)

    def _try_load(self, path: str | Path) -> None:
        try:
            self.load_registry(path)
        except (OSError, UnicodeError, ValueError, RuntimeError) as exc:
            self.message.setText(f"Reload rejected; previous snapshot retained: {exc}")

    def selected_case_id(self) -> str | None:
        row = self.table.currentRow()
        if row < 0:
            return None
        item = self.table.item(row, 0)
        return item.data(Qt.ItemDataRole.UserRole) if item is not None else None

    def refresh(self, *_args) -> None:
        previous = self.selected_case_id()
        self.table.setSortingEnabled(False)
        self.table.setRowCount(0)
        if self.registry is None:
            self._visible_cases = {}
            self.summary.setText("0 cases — no registry loaded")
            self.details.clear()
            return
        cases = self.registry.filter(
            self.search.text(),
            domain=self.domain_filter.currentData() or "",
            status=self.status_filter.currentData() or "",
        )
        self._visible_cases = {case.id: case for case in cases}
        self.table.setRowCount(len(cases))
        for row, case in enumerate(cases):
            proof = self._proofs.get(case.id, "UNVERIFIED")
            for col, value in enumerate(
                (case.id, case.name, case.domain, case.status, proof)
            ):
                item = QTableWidgetItem(value)
                if col == 0:
                    item.setData(Qt.ItemDataRole.UserRole, case.id)
                if col == 3:
                    item.setToolTip(
                        "Status from registry only; never a current-HEAD certificate"
                    )
                if col == 4:
                    item.setToolTip(
                        f"Claimed external CI proof at {self._proof_commit or 'unknown SHA'}; reverify artifacts on exact HEAD"
                    )
                self.table.setItem(row, col, item)
        self.table.setSortingEnabled(True)
        counts = self.registry.counts()
        count_text = " · ".join(
            f"{key}: {count}" for key, count in counts.items() if count
        )
        self.summary.setText(
            f"Showing {len(cases)} / {len(self.registry.cases)} cases · "
            f"Registry declarations: {count_text}"
        )
        target_row = next(
            (
                row
                for row in range(self.table.rowCount())
                if self.table.item(row, 0).data(Qt.ItemDataRole.UserRole) == previous
            ),
            0 if cases else -1,
        )
        if target_row >= 0:
            self.table.setCurrentCell(target_row, 0)
            self.table.selectRow(target_row)
        else:
            self.details.clear()

    def _show_details(self) -> None:
        case_id = self.selected_case_id()
        case = self._visible_cases.get(case_id) if case_id else None
        if case is None or self.registry is None:
            self.details.clear()
            return
        lines = [
            f"{case.id} — {case.name}",
            f"Domain: {case.domain}",
            f"Registry status: {case.status} (DECLARED, NOT CURRENT-HEAD CERTIFIED)",
            f"Definition: {self.registry.policy_for(case.status)}",
            "",
        ]
        lines.extend(
            f"{key.replace('_', ' ').capitalize()}: {value}"
            for key, value in case.details
        )
        self.details.setPlainText("\n".join(lines))
