"""Contract-driven Qt setup editor for CFDX cases."""
from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPushButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .boundary_setup import FIELD_SPECS, SCALAR_TYPES, VELOCITY_TYPES
from .case import Case
from .initialization import InitializationMode, InitializationSpec
from .materials import MaterialSpec
from .physics_setup import PHYSICS_SPECS, TURBULENCE_MODELS
from .setup_model import ParameterType


class CaseSetupPanel(QWidget):
    """Thin Qt adapter over the shared headless setup contracts."""

    changed = Signal()

    def __init__(self, case: Case, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.case = case
        tabs = QTabWidget(self)
        tabs.addTab(self._physics_tab(), "Physics")
        tabs.addTab(self._material_tab(), "Materials")
        tabs.addTab(self._boundary_tab(), "Boundaries")
        tabs.addTab(self._initialization_tab(), "Initialization")
        layout = QVBoxLayout(self)
        layout.addWidget(tabs)

    @staticmethod
    def _real(value: float) -> QDoubleSpinBox:
        box = QDoubleSpinBox()
        box.setRange(-1e300, 1e300)
        box.setDecimals(10)
        box.setValue(value)
        return box

    def _physics_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.physics_model = QComboBox()
        self.physics_model.addItems([spec.key for spec in PHYSICS_SPECS])
        self.physics_enabled = QCheckBox("Enabled")
        self.physics_button = QPushButton("Apply physics")
        self.physics_fields_form = QFormLayout()
        form.addRow("Model", self.physics_model)
        form.addRow("", self.physics_enabled)
        form.addRow(self.physics_fields_form)
        form.addRow("", self.physics_button)
        self.physics_fields: dict[str, QWidget] = {}

        self.physics_model.currentTextChanged.connect(self._rebuild_physics_fields)
        self.physics_enabled.toggled.connect(self._update_physics_dependencies)
        self.physics_button.clicked.connect(self._apply_physics)
        self._rebuild_physics_fields(self.physics_model.currentText())
        return w

    def _clear_form(self, form: QFormLayout) -> None:
        while form.count():
            item = form.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    def _rebuild_physics_fields(self, model: str) -> None:
        self._clear_form(self.physics_fields_form)
        self.physics_fields = {}
        try:
            spec = next(item for item in PHYSICS_SPECS if item.key == model)
        except StopIteration:
            return
        current = self.case.physics.get(model, {})
        for field in spec.fields:
            if field.name == "model":
                box = QComboBox()
                box.addItems(TURBULENCE_MODELS)
                box.setCurrentText(str(current.get(field.name, field.default)))
                widget: QWidget = box
            elif field.kind is ParameterType.REAL:
                widget = self._real(float(current.get(field.name, field.default)))
            else:
                box = QComboBox()
                box.addItem(str(field.default))
                box.setCurrentText(str(current.get(field.name, field.default)))
                widget = box
            self.physics_fields[field.name] = widget
            label = field.name.replace("_", " ").title()
            if field.unit:
                label += f" [{field.unit}]"
            self.physics_fields_form.addRow(label, widget)
        self._update_physics_dependencies()

    def _update_physics_dependencies(self) -> None:
        model = self.physics_model.currentText()
        enabled_models = {
            key for key, value in self.case.physics.items()
            if isinstance(value, dict) and bool(value.get("enabled", True))
        }
        for spec in PHYSICS_SPECS:
            if spec.key == model:
                dependencies_ok = all(
                    dependency in enabled_models or dependency == model
                    for dependency in spec.requires
                )
                self.physics_enabled.setEnabled(dependencies_ok)
                if not dependencies_ok:
                    self.physics_enabled.setToolTip(
                        "Enable required models first: " + ", ".join(spec.requires)
                    )
                else:
                    self.physics_enabled.setToolTip("")
                return

    def _apply_physics(self) -> None:
        model = self.physics_model.currentText()
        values: dict[str, object] = {"enabled": self.physics_enabled.isChecked()}
        for spec in PHYSICS_SPECS:
            if spec.key != model:
                continue
            for field in spec.fields:
                widget = self.physics_fields.get(field.name)
                if widget is None:
                    continue
                if isinstance(widget, QDoubleSpinBox):
                    values[field.name] = widget.value()
                elif isinstance(widget, QComboBox):
                    values[field.name] = widget.currentText()
            break
        self.case.physics[model] = values
        self._update_physics_dependencies()
        self.changed.emit()

    def _material_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.material_name = QLineEdit()
        self.material_density = self._real(1.0)
        self.material_viscosity = self._real(1e-3)
        self.material_cp = self._real(1000.0)
        self.material_conductivity = self._real(0.0)
        self.material_button = QPushButton("Apply material")
        form.addRow("Name", self.material_name)
        form.addRow("Density [kg/m³]", self.material_density)
        form.addRow("Dynamic viscosity [Pa·s]", self.material_viscosity)
        form.addRow("Cp [J/(kg·K)]", self.material_cp)
        form.addRow("Conductivity [W/(m·K)]", self.material_conductivity)
        form.addRow("", self.material_button)
        self.material_button.clicked.connect(self._apply_material)
        return w

    def _boundary_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.boundary_list = QListWidget()
        self.boundary_name = QLineEdit()
        self.boundary_type = QComboBox()
        self.boundary_type.addItems(
            ["inlet", "outlet", "wall", "symmetry", "periodic", "interface", "empty"]
        )
        self.boundary_field = QComboBox()
        self.boundary_field.addItems(list(FIELD_SPECS))
        self.boundary_scalar_type = QComboBox()
        self.boundary_scalar_type.addItems(SCALAR_TYPES)
        self.boundary_value = self._real(0.0)
        self.boundary_gradient = self._real(0.0)
        self.boundary_x = self._real(0.0)
        self.boundary_y = self._real(0.0)
        self.boundary_z = self._real(0.0)
        self.boundary_button = QPushButton("Apply / Update")
        self.boundary_remove_button = QPushButton("Remove selected")
        form.addRow("Existing patches", self.boundary_list)
        form.addRow("Name", self.boundary_name)
        form.addRow("Patch type", self.boundary_type)
        form.addRow("Field", self.boundary_field)
        form.addRow("Scalar condition", self.boundary_scalar_type)
        form.addRow("Value", self.boundary_value)
        form.addRow("Gradient", self.boundary_gradient)
        form.addRow("Vector X [m/s]", self.boundary_x)
        form.addRow("Vector Y [m/s]", self.boundary_y)
        form.addRow("Vector Z [m/s]", self.boundary_z)
        row = QHBoxLayout()
        row.addWidget(self.boundary_button)
        row.addWidget(self.boundary_remove_button)
        form.addRow("", row)
        self.boundary_button.clicked.connect(self._apply_boundary)
        self.boundary_remove_button.clicked.connect(self._remove_boundary)
        self.boundary_list.currentTextChanged.connect(self._load_boundary)
        self.boundary_field.currentTextChanged.connect(self._update_boundary_field_controls)
        self._update_boundary_field_controls(self.boundary_field.currentText())
        self._refresh_boundaries()
        return w

    def _initialization_tab(self) -> QWidget:
        w = QWidget()
        form = QFormLayout(w)
        self.initialization_mode = QComboBox()
        self.initialization_mode.addItems([m.value for m in InitializationMode])
        self.initialization_field = QLineEdit("U")
        self.initialization_value = self._real(0.0)
        self.initialization_button = QPushButton("Apply initialization")
        form.addRow("Mode", self.initialization_mode)
        form.addRow("Field", self.initialization_field)
        form.addRow("Uniform value", self.initialization_value)
        form.addRow("", self.initialization_button)
        self.initialization_button.clicked.connect(self._apply_initialization)
        self.initialization_mode.currentTextChanged.connect(
            lambda mode: self.initialization_value.setEnabled(mode == "uniform")
        )
        return w

    def set_case(self, case: Case) -> None:
        self.case = case
        self._rebuild_physics_fields(self.physics_model.currentText())
        self._refresh_boundaries()

    def _apply_material(self) -> None:
        name = self.material_name.text().strip()
        if not name:
            QMessageBox.warning(self, "Material", "Material name is required")
            return
        material = MaterialSpec(
            name,
            self.material_density.value(),
            self.material_viscosity.value(),
            self.material_cp.value(),
            self.material_conductivity.value(),
        )
        try:
            material.validate()
        except ValueError as exc:
            QMessageBox.warning(self, "Material validation", str(exc))
            return
        self.case.materials[name] = {
            "density": material.density,
            "dynamic_viscosity": material.dynamic_viscosity,
            "cp": material.cp,
            "conductivity": material.conductivity,
        }
        self.changed.emit()

    def _update_boundary_field_controls(self, field: str) -> None:
        vector = field == "velocity"
        for box in (self.boundary_x, self.boundary_y, self.boundary_z):
            box.setVisible(vector)
        self.boundary_value.setVisible(not vector)
        self.boundary_gradient.setVisible(not vector)
        self.boundary_scalar_type.clear()
        self.boundary_scalar_type.addItems(VELOCITY_TYPES if vector else SCALAR_TYPES)

    def _apply_boundary(self) -> None:
        name = self.boundary_name.text().strip()
        if not name:
            QMessageBox.warning(self, "Boundary", "Boundary name is required")
            return
        values: dict[str, object] = {"type": self.boundary_type.currentText()}
        if values["type"] not in {"periodic", "interface", "empty"}:
            field = self.boundary_field.currentText()
            values["fields"] = [field]
            if field == "velocity":
                values.update(
                    {
                        "velocity_type": self.boundary_scalar_type.currentText(),
                        "velocity_value": [
                            self.boundary_x.value(),
                            self.boundary_y.value(),
                            self.boundary_z.value(),
                        ],
                    }
                )
            else:
                values.update(
                    {
                        "scalar_type": self.boundary_scalar_type.currentText(),
                        "value": self.boundary_value.value(),
                        "gradient": self.boundary_gradient.value(),
                    }
                )
        self.case.set_boundary(name, **values)
        self._refresh_boundaries()
        self.changed.emit()

    def _remove_boundary(self) -> None:
        item = self.boundary_list.currentItem()
        if item is not None:
            self.case.boundaries.pop(item.text(), None)
            self._refresh_boundaries()
            self.changed.emit()

    def _refresh_boundaries(self) -> None:
        self.boundary_list.blockSignals(True)
        self.boundary_list.clear()
        self.boundary_list.addItems(list(self.case.boundaries))
        self.boundary_list.blockSignals(False)

    def _load_boundary(self, name: str) -> None:
        if not name:
            return
        boundary = self.case.boundaries.get(name, {})
        index = self.boundary_type.findText(str(boundary.get("type", "wall")))
        self.boundary_type.setCurrentIndex(max(0, index))
        field = str((boundary.get("fields") or ["pressure"])[0])
        index = self.boundary_field.findText(field)
        self.boundary_field.setCurrentIndex(max(0, index))
        index = self.boundary_scalar_type.findText(
            str(boundary.get("scalar_type", boundary.get("velocity_type", "ZERO_GRADIENT")))
        )
        self.boundary_scalar_type.setCurrentIndex(max(0, index))
        vector = boundary.get("velocity_value", [0.0, 0.0, 0.0])
        self.boundary_x.setValue(float(vector[0]))
        self.boundary_y.setValue(float(vector[1]))
        self.boundary_z.setValue(float(vector[2]))
        self.boundary_value.setValue(float(boundary.get("value", 0.0)))
        self.boundary_gradient.setValue(float(boundary.get("gradient", 0.0)))
        self.boundary_name.setText(name)

    def _apply_initialization(self) -> None:
        mode = InitializationMode(self.initialization_mode.currentText())
        field = self.initialization_field.text().strip()
        value = self.initialization_value.value() if mode is InitializationMode.UNIFORM else None
        spec = InitializationSpec(mode, field, value)
        try:
            spec.validate()
        except ValueError as exc:
            QMessageBox.warning(self, "Initialization validation", str(exc))
            return
        self.case.physics["initialization"] = {
            "mode": mode.value,
            "field": spec.field,
            "value": spec.value,
        }
        self.changed.emit()
