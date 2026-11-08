"""Contract-driven Qt setup editor for CFDX cases."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QMessageBox,
    QPushButton,
    QTableWidget, QTableWidgetItem, QAbstractItemView,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from .boundary_setup import FIELD_SPECS, SCALAR_TYPES, VELOCITY_TYPES
from .boundary_matrix import FIELD_UNITS, field_rows, remove_field_condition, setup_advice
from .boundary_editor import (
    BOUNDARY_PRESETS, apply_boundary_condition, boundary_runtime_support,
    condition_for_field,
)
from .case import Case
from .initialization import InitializationMode, InitializationSpec
from .materials import MaterialSpec
from .physics_setup import (
    PHYSICS_SPECS,
    TURBULENCE_CATALOG,
    turbulence_model_from_case,
    validate_turbulence_selection,
)
from .setup_model import ParameterType


class CaseSetupPanel(QWidget):
    """Thin Qt adapter over the shared headless setup contracts."""

    changed = Signal()

    def __init__(self, case: Case, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.case = case
        self.tabs = QTabWidget(self)
        self.tabs.addTab(self._physics_tab(), "Physics")
        self.tabs.addTab(self._material_tab(), "Materials")
        self.tabs.addTab(self._boundary_tab(), "Boundaries")
        self.tabs.addTab(self._initialization_tab(), "Initialization")
        layout = QVBoxLayout(self)
        self.validation_summary = QLabel("No validation diagnostics")
        self.validation_summary.setWordWrap(True)
        layout.addWidget(self.validation_summary)
        layout.addWidget(self.tabs)

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
        if model == "turbulence":
            normalized = turbulence_model_from_case(current)
            self.turbulence_family = QComboBox()
            families = sorted(
                {item.family.value.upper() for item in TURBULENCE_CATALOG}
            )
            self.turbulence_family.addItems(families)
            self.turbulence_family.setCurrentText(str(normalized.get("family", "RANS")))
            self.physics_fields["family"] = self.turbulence_family
            self.physics_fields_form.addRow("Turbulence family", self.turbulence_family)

            self.turbulence_wall = QComboBox()
            self.turbulence_wall.addItems(
                ["resolved", "wall_function", "all_y_plus", "none"]
            )
            self.turbulence_wall.setCurrentText(
                str(normalized.get("wall_treatment", "resolved"))
            )
            self.physics_fields["wall_treatment"] = self.turbulence_wall
            self.physics_fields_form.addRow("Wall treatment", self.turbulence_wall)

            self.turbulence_family.currentTextChanged.connect(
                self._filter_turbulence_models
            )
        for field in spec.fields:
            if field.name == "model":
                box = QComboBox()
                self.physics_fields["model"] = box
                self._filter_turbulence_models(
                    str(normalized.get("model", field.default))
                    if model == "turbulence"
                    else str(current.get(field.name, field.default))
                )
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

    def _filter_turbulence_models(self, preferred: str | None = None) -> None:
        box = self.physics_fields.get("model")
        family_box = self.physics_fields.get("family")
        if not isinstance(box, QComboBox) or not isinstance(family_box, QComboBox):
            return
        family = family_box.currentText().lower()
        choices = [
            item.key
            for item in TURBULENCE_CATALOG
            if item.status.value == "solver_ready" and item.family.value == family
        ]
        current = preferred or box.currentText()
        box.blockSignals(True)
        box.clear()
        box.addItems(choices)
        if current in choices:
            box.setCurrentText(current)
        box.blockSignals(False)

    def _update_physics_dependencies(self) -> None:
        model = self.physics_model.currentText()
        enabled_models = {
            key
            for key, value in self.case.physics.items()
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
        if model == "turbulence":
            try:
                values = validate_turbulence_selection(values)
            except ValueError as exc:
                QMessageBox.warning(self, "Turbulence selection", str(exc))
                return
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
        self.boundary_search = QLineEdit()
        self.boundary_search.setObjectName("workbench.boundaries.search")
        self.boundary_search.setPlaceholderText("Filter boundary patches...")
        self.boundary_preset = QComboBox()
        self.boundary_preset.setObjectName("workbench.boundaries.preset")
        self.boundary_preset.addItem("Choose a suggested setup (never auto-applied)", -1)
        for index, preset in enumerate(BOUNDARY_PRESETS):
            self.boundary_preset.addItem(preset.name, index)
        self.boundary_preset_button = QPushButton("Use suggested setup")
        self.boundary_preset_button.setObjectName("workbench.boundaries.use_preset")
        self.boundary_note = QLabel()
        self.boundary_note.setWordWrap(True)
        self.boundary_note.setObjectName("workbench.boundaries.note")
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
        self.boundary_field_remove_button = QPushButton("Remove selected field")
        self.boundary_field_remove_button.setObjectName("workbench.boundaries.remove_field")
        self.boundary_matrix = QTableWidget(0, 4, w)
        self.boundary_matrix.setObjectName("workbench.boundaries.field_matrix")
        self.boundary_matrix.setHorizontalHeaderLabels(("Field", "Condition", "Value", "Execution status"))
        self.boundary_matrix.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.boundary_matrix.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.boundary_matrix.horizontalHeader().setStretchLastSection(True)
        self.boundary_advice = QLabel("Choose a patch to inspect its fields")
        self.boundary_advice.setObjectName("workbench.boundaries.advice")
        self.boundary_advice.setWordWrap(True)
        self.boundary_remove_button = QPushButton("Remove selected")
        form.addRow("Suggestions", self.boundary_preset)
        form.addRow("", self.boundary_preset_button)
        form.addRow("Find patch", self.boundary_search)
        form.addRow("Existing patches", self.boundary_list)
        form.addRow("Name", self.boundary_name)
        form.addRow("Patch type", self.boundary_type)
        form.addRow("Field", self.boundary_field)
        form.addRow("Scalar condition", self.boundary_scalar_type)
        self.boundary_value_label = QLabel("Value")
        self.boundary_gradient_label = QLabel("Gradient")
        form.addRow(self.boundary_value_label, self.boundary_value)
        form.addRow(self.boundary_gradient_label, self.boundary_gradient)
        form.addRow("Vector X [m/s]", self.boundary_x)
        form.addRow("Vector Y [m/s]", self.boundary_y)
        form.addRow("Vector Z [m/s]", self.boundary_z)
        form.addRow("Configured fields", self.boundary_matrix)
        form.addRow("", self.boundary_field_remove_button)
        form.addRow("Field / runtime support", self.boundary_note)
        form.addRow("Engineering checks", self.boundary_advice)
        row = QHBoxLayout()
        row.addWidget(self.boundary_button)
        row.addWidget(self.boundary_remove_button)
        form.addRow("", row)
        self.boundary_button.clicked.connect(self._apply_boundary)
        self.boundary_field_remove_button.clicked.connect(self._remove_boundary_field)
        self.boundary_matrix.itemSelectionChanged.connect(self._selected_field_in_matrix)
        self.boundary_remove_button.clicked.connect(self._remove_boundary)
        self.boundary_list.currentTextChanged.connect(self._load_boundary)
        self.boundary_search.textChanged.connect(self._filter_boundary_list)
        self.boundary_preset_button.clicked.connect(self._select_boundary_preset)
        self.boundary_type.currentTextChanged.connect(self._update_boundary_guidance)
        self.boundary_field.currentTextChanged.connect(
            self._update_boundary_field_controls
        )
        self._update_boundary_field_controls(self.boundary_field.currentText())
        self._update_boundary_guidance()
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
        self._refresh_boundary_matrix()

    @staticmethod
    def _set_validation_state(widget: QWidget, state: str) -> None:
        """Apply a semantic validation state through the shared Qt stylesheet."""
        widget.setProperty("validationState", state)
        style = widget.style()
        style.unpolish(widget)
        style.polish(widget)
        widget.update()

    def set_diagnostics(self, diagnostics) -> None:
        """Show feedback, select the relevant tab, and highlight its editor."""
        diagnostics = tuple(diagnostics)
        widgets = (
            *self.physics_fields.values(),
            self.material_name,
            self.material_density,
            self.material_viscosity,
            self.material_cp,
            self.material_conductivity,
            self.boundary_name,
            self.boundary_type,
            self.boundary_field,
            self.boundary_scalar_type,
            self.boundary_value,
            self.boundary_gradient,
            self.boundary_x,
            self.boundary_y,
            self.boundary_z,
            self.initialization_mode,
            self.initialization_field,
            self.initialization_value,
        )
        self._set_validation_state(self.validation_summary, "")
        for widget in widgets:
            self._set_validation_state(widget, "")
            widget.setToolTip("")

        if not diagnostics:
            self.validation_summary.setText("No validation diagnostics")
            return

        errors = sum(item.severity == "error" for item in diagnostics)
        self.validation_summary.setText(
            f"Validation: {errors} error(s), {len(diagnostics) - errors} warning(s)"
            + " — "
            + "; ".join(item.message for item in diagnostics[:3])
        )
        self._set_validation_state(
            self.validation_summary, "error" if errors else "warning"
        )
        for item in diagnostics:
            target, tab_index = self._diagnostic_target(item.path)
            if target is not None:
                self.tabs.setCurrentIndex(tab_index)
                target.setToolTip(f"{item.code}: {item.message}")
                self._set_validation_state(
                    target, "error" if item.severity == "error" else "warning"
                )

    def _diagnostic_target(self, path: str):
        """Resolve a validation path to a widget and its Setup tab."""
        parts = path.split(".")
        if parts[:1] == ["physics"]:
            if len(parts) == 2 and parts[1] == "initialization":
                return self.initialization_mode, 3
            if len(parts) >= 3:
                model, field = parts[1], parts[2]
                if model != self.physics_model.currentText():
                    self.physics_model.setCurrentText(model)
                return self.physics_fields.get(field), 0
            return None, 0
        if parts[:1] == ["materials"]:
            field = parts[2] if len(parts) >= 3 else None
            material_widgets = {
                "density": self.material_density,
                "dynamic_viscosity": self.material_viscosity,
                "cp": self.material_cp,
                "conductivity": self.material_conductivity,
            }
            return material_widgets.get(field, self.material_name), 1
        if parts[:1] == ["boundaries"]:
            if len(parts) >= 2:
                name = parts[1]
                row = self.boundary_list.findItems(name, Qt.MatchFlag.MatchExactly)
                if row:
                    self.boundary_list.setCurrentItem(row[0])
            field = parts[2] if len(parts) >= 3 else None
            boundary_widgets = {
                "type": self.boundary_type,
                "fields": self.boundary_field,
                "scalar_type": self.boundary_scalar_type,
                "value": self.boundary_value,
                "gradient": self.boundary_gradient,
                "velocity_type": self.boundary_scalar_type,
            }
            return boundary_widgets.get(field, self.boundary_name), 2
        return None, 0

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
        self.boundary_scalar_type.blockSignals(True)
        previous = self.boundary_scalar_type.currentText()
        self.boundary_scalar_type.clear()
        choices = (VELOCITY_TYPES if vector else SCALAR_TYPES +
                   (("WALL_K_ZERO",) if field == "k" else
                    ("WALL_OMEGA_SST",) if field == "omega" else ()))
        self.boundary_scalar_type.addItems(choices)
        if previous in choices:
            self.boundary_scalar_type.setCurrentText(previous)
        self.boundary_scalar_type.blockSignals(False)
        unit = FIELD_UNITS.get(field, "")
        self.boundary_value_label.setText(f"Value [{unit}]")
        self.boundary_gradient_label.setText(f"Normal gradient [{unit}/m]")
        self._update_boundary_guidance()
        self._load_selected_field_condition()

    def _update_boundary_guidance(self, *_args) -> None:
        kind = self.boundary_type.currentText()
        is_geometric = kind in {"symmetry", "periodic", "interface", "empty"}
        self.boundary_field.setEnabled(not is_geometric)
        self.boundary_scalar_type.setEnabled(not is_geometric)
        for editor in (
            self.boundary_value, self.boundary_gradient,
            self.boundary_x, self.boundary_y, self.boundary_z,
        ):
            editor.setEnabled(not is_geometric)
        selected = self.boundary_name.text().strip()
        existing = self.case.boundaries.get(selected, {})
        if existing:
            supported, message = boundary_runtime_support(existing)
            self.boundary_note.setText(
                ("Editor schema accepted; native mapping unverified: "
                 if supported else "Native execution blocked: ") + message
            )
            self._set_validation_state(self.boundary_note, "" if supported else "warning")
        else:
            self.boundary_note.setText(
                "Geometry-only patch: no field conditions" if is_geometric else
                "Single-field CFDX solver contract. Additional fields can be staged "
                "without data loss, but will block execution until natively mapped."
            )
            self._set_validation_state(self.boundary_note, "")

    def _refresh_boundary_matrix(self) -> None:
        name = self.boundary_name.text().strip()
        self.boundary_matrix.setRowCount(0)
        boundary = self.case.boundaries.get(name)
        if isinstance(boundary, dict):
            try:
                rows = field_rows(name, boundary)
            except ValueError as exc:
                self.boundary_advice.setText(f"Invalid stored boundary fields: {exc}")
                self.boundary_field_remove_button.setEnabled(False)
                return
            self.boundary_matrix.setRowCount(len(rows))
            for row, item in enumerate(rows):
                for column, value in enumerate((
                    item.field, item.mode, f"{item.value} {item.units}",
                    "editor-valid / native unverified" if item.native_supported
                    else "staged / BLOCKED",
                )):
                    self.boundary_matrix.setItem(row, column, QTableWidgetItem(value))
        self.boundary_field_remove_button.setEnabled(False)
        advice = setup_advice(self.case.boundaries)
        self.boundary_advice.setText("; ".join(advice) if advice else "No additional editor warnings")

    def _selected_field_in_matrix(self) -> None:
        row = self.boundary_matrix.currentRow()
        self.boundary_field_remove_button.setEnabled(row >= 0)
        if row >= 0 and self.boundary_matrix.item(row, 0) is not None:
            field = self.boundary_matrix.item(row, 0).text()
            if field != self.boundary_field.currentText():
                self.boundary_field.setCurrentText(field)

    def _remove_boundary_field(self) -> None:
        name = self.boundary_name.text().strip()
        row = self.boundary_matrix.currentRow()
        if row < 0 or self.boundary_matrix.item(row, 0) is None:
            return
        field = self.boundary_matrix.item(row, 0).text()
        response = QMessageBox.question(
            self, "Remove field condition", f"Remove {field!r} from patch {name!r}?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if response != QMessageBox.StandardButton.Yes:
            return
        try:
            new_boundary = remove_field_condition(self.case.boundaries[name], field)
        except (KeyError, ValueError) as exc:
            QMessageBox.warning(self, "Boundary validation", str(exc))
            return
        self.case.set_boundary(name, **new_boundary)
        self._load_boundary(name)
        self._refresh_boundary_matrix()
        self.changed.emit()

    def _filter_boundary_list(self, query: str) -> None:
        needle = query.casefold().strip()
        for row in range(self.boundary_list.count()):
            item = self.boundary_list.item(row)
            item.setHidden(needle not in item.text().casefold())

    def _select_boundary_preset(self) -> None:
        index = self.boundary_preset.currentData()
        if not isinstance(index, int) or index < 0:
            return
        preset = BOUNDARY_PRESETS[index]
        self.boundary_type.setCurrentText(preset.patch_type)
        if preset.field is not None:
            self.boundary_field.setCurrentText(preset.field)
            self.boundary_scalar_type.setCurrentText(preset.condition or "ZERO_GRADIENT")
        if preset.name == "No-slip wall":
            for box in (self.boundary_x, self.boundary_y, self.boundary_z):
                box.setValue(0.0)
        self.boundary_note.setText(preset.explanation + " Click Apply / Update to save.")

    def _apply_boundary(self) -> None:
        name = self.boundary_name.text().strip()
        if not name:
            QMessageBox.warning(self, "Boundary", "Boundary name is required")
            return
        existing = self.case.boundaries.get(name, {})
        kind = self.boundary_type.currentText()
        changed_type = bool(existing and existing.get("type") != kind)
        if changed_type:
            response = QMessageBox.question(
                self, "Change boundary kind",
                "Changing this patch kind will discard the previous field conditions. "
                "Continue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if response != QMessageBox.StandardButton.Yes:
                return
        try:
            values = apply_boundary_condition(
                existing,
                patch_type=kind,
                field=self.boundary_field.currentText(),
                mode=self.boundary_scalar_type.currentText(),
                value=self.boundary_value.value(),
                gradient=self.boundary_gradient.value(),
                velocity=(self.boundary_x.value(), self.boundary_y.value(), self.boundary_z.value()),
                allow_type_change=changed_type,
            )
        except ValueError as exc:
            QMessageBox.warning(self, "Boundary validation", str(exc))
            return
        # No partial update: the validated copy is committed only once.
        self.case.set_boundary(name, **values)
        self._refresh_boundaries(select=name)
        self._update_boundary_guidance()
        self._refresh_boundary_matrix()
        self.changed.emit()

    def _remove_boundary(self) -> None:
        item = self.boundary_list.currentItem()
        if item is None:
            return
        response = QMessageBox.question(
            self, "Remove boundary", f"Remove patch configuration {item.text()!r}?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if response != QMessageBox.StandardButton.Yes:
            return
        self.case.boundaries.pop(item.text(), None)
        self._refresh_boundaries()
        self._refresh_boundary_matrix()
        self.changed.emit()

    def _refresh_boundaries(self, *, select: str | None = None) -> None:
        selected = select or (
            self.boundary_list.currentItem().text()
            if self.boundary_list.currentItem() is not None else None
        )
        self.boundary_list.blockSignals(True)
        self.boundary_list.clear()
        self.boundary_list.addItems(list(self.case.boundaries))
        self.boundary_list.blockSignals(False)
        self._filter_boundary_list(self.boundary_search.text())
        if selected is not None:
            rows = self.boundary_list.findItems(selected, Qt.MatchFlag.MatchExactly)
            if rows:
                self.boundary_list.setCurrentItem(rows[0])

    def _load_selected_field_condition(self) -> None:
        name = self.boundary_name.text().strip()
        boundary = self.case.boundaries.get(name)
        if not isinstance(boundary, dict):
            return
        field = self.boundary_field.currentText()
        condition = condition_for_field(boundary, field)
        # Switching to a not-yet-configured field must not inherit the last
        # field's mode. In particular, an inlet pressure preset must not leak
        # FIXED_VALUE into a newly selected temperature condition.
        mode = condition.get("velocity_type" if field == "velocity" else "scalar_type", "ZERO_GRADIENT")
        if isinstance(mode, str):
            self.boundary_scalar_type.setCurrentText(mode)
        vector = condition.get("velocity_value", (0.0, 0.0, 0.0))
        if isinstance(vector, (list, tuple)) and len(vector) == 3:
            for spin, value in zip((self.boundary_x, self.boundary_y, self.boundary_z), vector):
                spin.setValue(float(value))
        self.boundary_value.setValue(float(condition.get("value", 0.0)))
        self.boundary_gradient.setValue(float(condition.get("gradient", 0.0)))

    def _load_boundary(self, name: str) -> None:
        if not name:
            return
        boundary = self.case.boundaries.get(name, {})
        self.boundary_name.setText(name)
        self.boundary_type.setCurrentText(str(boundary.get("type", "wall")))
        fields = boundary.get("fields") or ("pressure",)
        self.boundary_field.setCurrentText(str(fields[0]))
        self._load_selected_field_condition()
        self._update_boundary_guidance()
        self._refresh_boundary_matrix()

    def _apply_initialization(self) -> None:
        mode = InitializationMode(self.initialization_mode.currentText())
        field = self.initialization_field.text().strip()
        value = (
            self.initialization_value.value()
            if mode is InitializationMode.UNIFORM
            else None
        )
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
