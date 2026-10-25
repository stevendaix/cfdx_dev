"""Typed contextual property inspector for the Workbench.

The panel is a thin Qt adapter over the shared PropertyState contracts:
each row's declared type picks its editor widget (real, integer, boolean
or string), the label carries the unit when the contract provides one,
and read-only rows can never offer an edit the domain does not support.
"""

from __future__ import annotations

from collections.abc import Callable

from PySide6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFormLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from .application.properties import PropertyState
from .setup_model import ChangeImpact, ParameterType

_REAL_MIN = -1e300
_REAL_MAX = 1e300
_INTEGER_MIN = -2_000_000_000
_INTEGER_MAX = 2_000_000_000


class PropertiesPanel(QWidget):
    """Typed editors for the properties of the current selection."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("workbench.properties")
        self.on_apply: Callable[[str, object], None] | None = None
        self._form = QFormLayout()
        layout = QVBoxLayout(self)
        layout.addLayout(self._form)
        self._form.addRow(QLabel("Select a workflow item"))

    def set_properties(self, properties: tuple[PropertyState, ...]) -> None:
        """Rebuild one typed editor row per declared property."""
        self._clear_form()
        if not properties:
            self._form.addRow(QLabel("No editable properties for this selection"))
            return
        for prop in properties:
            self._add_property(prop)

    def _clear_form(self) -> None:
        while self._form.count():
            item = self._form.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()

    def _add_property(self, prop: PropertyState) -> None:
        """Add the typed editor, unit label and Apply button of one row."""
        editor, read_value = self._editor_for(prop)
        editor.setObjectName(f"workbench.property.{prop.key}")
        self._form.addRow(self._label_for(prop), editor)
        self._describe(prop, editor)
        if not prop.editable:
            editor.setEnabled(False)
            return
        button = QPushButton("Apply")
        button.setObjectName(f"workbench.property.apply.{prop.key}")
        button.clicked.connect(self._make_apply_handler(prop.key, read_value))
        self._form.addRow("", button)

    def _make_apply_handler(
        self, key: str, read_value: Callable[[], object]
    ) -> Callable[[], None]:
        """Bind one property's typed reader to its Apply button."""

        def apply_current() -> None:
            if self.on_apply is not None:
                self.on_apply(key, read_value())

        return apply_current

    def _editor_for(self, prop: PropertyState) -> tuple[QWidget, Callable[[], object]]:
        """Build the typed editor for one declared property."""
        if prop.kind is ParameterType.REAL:
            editor = QDoubleSpinBox()
            editor.setRange(_REAL_MIN, _REAL_MAX)
            editor.setDecimals(10)
            editor.setValue(float(prop.value))
            return editor, editor.value
        if prop.kind is ParameterType.INTEGER:
            editor = QSpinBox()
            editor.setRange(_INTEGER_MIN, _INTEGER_MAX)
            editor.setValue(int(prop.value))
            return editor, editor.value
        if prop.kind is ParameterType.BOOL:
            editor = QCheckBox()
            editor.setChecked(bool(prop.value))
            return editor, editor.isChecked
        line = QLineEdit(str(prop.value))
        if prop.kind is not ParameterType.STRING:
            line.setReadOnly(True)
        return line, line.text

    def _describe(self, prop: PropertyState, editor: QWidget) -> None:
        """Surface the declared impact and message as the editor tooltip."""
        if prop.impact is ChangeImpact.REQUIRES_RESTART:
            editor.setToolTip("Applying this change requires a solver restart")
        if prop.message is not None:
            editor.setToolTip(prop.message)

    @staticmethod
    def _label_for(prop: PropertyState) -> str:
        if prop.unit is None:
            return prop.label
        return f"{prop.label} [{prop.unit}]"
