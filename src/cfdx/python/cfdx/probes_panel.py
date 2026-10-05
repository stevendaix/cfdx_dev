"""Qt adapter for the application-level point-probe catalogue."""
from __future__ import annotations

from .application import Application, ApplicationStateChanged
from .application.commands import ConfigureProbes
from .probe import Probe, ProbeCatalog, ProbeField

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QComboBox,
    QHeaderView,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

_FIELD_CHOICES = [field.cli for field in ProbeField]


class ProbesPanel(QWidget):
    """Editor for the case's point-probe catalogue.

    The catalogue is committed to the case through the application's
    :class:`~cfdx.application.commands.ConfigureProbes` command; the execution
    backend then forwards probes to the solver as ``--probe``/``--probe-csv``
    arguments, producing the probe CSV at run time.
    """

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
        self._table.selectionModel().selectionChanged.connect(self._on_selection_changed)

        self._apply_button = QPushButton("Apply to case")
        self._apply_button.setEnabled(False)
        self._apply_button.clicked.connect(self._apply)

        self._dirty = False

        controls = QHBoxLayout()
        for button in (self._add_button, self._remove_button, self._apply_button):
            controls.addWidget(button)

        layout = QVBoxLayout(self)
        layout.addLayout(controls)
        layout.addWidget(self._table)

        self.application.events.subscribe(ApplicationStateChanged, self._refresh)
        self._refresh()

    def _refresh(self, *_args: object) -> None:
        """Rebuild the table from the case's current probe catalogue."""
        self._table.blockSignals(True)
        try:
            self._table.setRowCount(0)
            for probe in self.application.session.case.probes:
                self._add_row(probe)
        finally:
            self._table.blockSignals(False)
        self._set_dirty(False)

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
        self._remove_button.setEnabled(bool(self._table.selectionModel().selectedRows()))

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
            catalog.add(Probe(name=name, x=x, y=y, z=z, field=ProbeField.from_cli(field_token)))
        try:
            self.application.execute(ConfigureProbes(catalog))
        except ValueError as exc:
            parent = self.window()
            if parent is not None and hasattr(parent, "statusBar"):
                parent.statusBar().showMessage(f"Probes: {exc}")
        self._set_dirty(False)
