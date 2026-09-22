"""Reusable Qt setup controls for CFDX case physics and boundaries."""
from __future__ import annotations

from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from .case import Case


class CaseSetupPanel(QWidget):
    """Edit the public Python Case physics and boundary dictionaries."""

    def __init__(self, case: Case, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.case = case

        self.physics_model = QComboBox()
        self.physics_model.addItems(["energy", "laminar", "k-epsilon", "k-omega"])
        self.physics_enabled = QCheckBox("Enabled")
        self.physics_button = QPushButton("Apply physics")

        physics_form = QFormLayout()
        physics_form.addRow("Model", self.physics_model)
        physics_form.addRow("", self.physics_enabled)
        physics_form.addRow("", self.physics_button)

        self.boundary_name = QLineEdit()
        self.boundary_type = QComboBox()
        self.boundary_type.addItems(["inlet", "outlet", "wall", "symmetry"])
        self.boundary_value = QLineEdit()
        self.boundary_button = QPushButton("Apply boundary")

        boundary_form = QFormLayout()
        boundary_form.addRow("Name", self.boundary_name)
        boundary_form.addRow("Type", self.boundary_type)
        boundary_form.addRow("Value", self.boundary_value)
        boundary_form.addRow("", self.boundary_button)

        layout = QVBoxLayout(self)
        layout.addLayout(physics_form)
        layout.addLayout(boundary_form)

        self.physics_button.clicked.connect(self._apply_physics)
        self.boundary_button.clicked.connect(self._apply_boundary)

    def _apply_physics(self) -> None:
        model = self.physics_model.currentText()
        self.case.physics[model] = {"enabled": self.physics_enabled.isChecked()}

    def _apply_boundary(self) -> None:
        name = self.boundary_name.text().strip()
        if not name:
            return
        self.case.set_boundary(
            name,
            type=self.boundary_type.currentText(),
            value=self.boundary_value.text(),
        )
