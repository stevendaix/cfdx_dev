"""Read-only numerical method explorer; never changes production solver setup."""

from __future__ import annotations

from PySide6.QtWidgets import (
    QComboBox, QFormLayout, QHBoxLayout, QLabel, QLineEdit, QPlainTextEdit,
    QPushButton, QVBoxLayout, QWidget,
)

from .modern_methods import (
    METHODS, compare_methods, method_by_key, method_families, search_methods,
)


class MethodsPanel(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("workbench.methods.panel")
        layout = QVBoxLayout(self)
        caution = QLabel(
            "Technique explorer — read-only. An adapter mapping is NOT a "
            "certified production implementation. No solver setting is changed."
        )
        caution.setWordWrap(True)
        caution.setObjectName("workbench.methods.provenance")
        layout.addWidget(caution)
        form = QFormLayout()
        self.search = QLineEdit()
        self.search.setObjectName("workbench.methods.search")
        self.search.setPlaceholderText("Search method, preconditioner or qualification...")
        layout.addWidget(self.search)
        self.integration = QComboBox()
        self.integration.setObjectName("workbench.methods.integration")
        self.integration.addItem("All integration levels", "all")
        self.integration.addItem("Adapter mapped (not certified)", "mapped")
        self.integration.addItem("Research / experimental", "research")
        self.family = QComboBox()
        self.family.setObjectName("workbench.methods.family")
        self.family.addItem("All families", "")
        for family in method_families():
            self.family.addItem(family, family)
        self.method = QComboBox()
        self.method.setObjectName("workbench.methods.method")
        form.addRow("Integration level", self.integration)
        form.addRow("Method family", self.family)
        form.addRow("Method", self.method)
        layout.addLayout(form)
        self.details = QPlainTextEdit()
        self.details.setObjectName("workbench.methods.details")
        self.details.setReadOnly(True)
        layout.addWidget(self.details)
        # Side-by-side textual comparison is read-only and never selects a backend.
        compare_row = QHBoxLayout()
        self.compare_with = QComboBox()
        self.compare_with.setObjectName("workbench.methods.compare_with")
        self.compare_with.addItem("Compare with another method...", "")
        for item in METHODS:
            self.compare_with.addItem(item.label, item.key)
        self.compare_button = QPushButton("Compare methods")
        self.compare_button.setObjectName("workbench.methods.compare")
        self.compare_button.setEnabled(False)
        compare_row.addWidget(self.compare_with, 1)
        compare_row.addWidget(self.compare_button)
        layout.addLayout(compare_row)
        self.comparison = QPlainTextEdit()
        self.comparison.setObjectName("workbench.methods.comparison")
        self.comparison.setReadOnly(True)
        self.comparison.setVisible(False)
        layout.addWidget(self.comparison)
        self.family.currentTextChanged.connect(self._populate)
        self.integration.currentIndexChanged.connect(self._populate)
        self.search.textChanged.connect(self._populate)
        self.method.currentIndexChanged.connect(self._update)
        self.compare_with.currentIndexChanged.connect(self._update_comparison_button)
        self.compare_button.clicked.connect(self._compare)
        self._populate(self.family.currentText())

    def _populate(self, *_args) -> None:
        previous_key = self.method.currentData()
        self.method.blockSignals(True)
        self.method.clear()
        filtered = search_methods(
            self.search.text(),
            family=self.family.currentData() or "",
            integration=self.integration.currentData() or "all",
        )
        for item in filtered:
            self.method.addItem(item.label, item.key)
        if isinstance(previous_key, str):
            index = self.method.findData(previous_key)
            if index >= 0:
                self.method.setCurrentIndex(index)
        self.method.blockSignals(False)
        self._update()

    def _update(self, *_args) -> None:
        key = self.method.currentData()
        if not isinstance(key, str):
            self.details.setPlainText("No methods match the current filters.")
            self._update_comparison_button()
            return
        method = method_by_key(key)
        self.details.setPlainText(
            f"{method.label}\n\nIntegration status: {method.integration}\n\n"
            f"Purpose\n{method.purpose}\n\n"
            f"How to qualify\n{method.qualification}\n\n"
            f"Reference\n{method.reference}"
        )

        self._update_comparison_button()

    def _update_comparison_button(self, *_args) -> None:
        left, right = self.method.currentData(), self.compare_with.currentData()
        self.compare_button.setEnabled(
            isinstance(left, str) and isinstance(right, str)
            and bool(right) and left != right
        )

    def _compare(self) -> None:
        if not self.compare_button.isEnabled():
            return
        self.comparison.setPlainText(
            compare_methods(self.method.currentData(), self.compare_with.currentData())
        )
        self.comparison.setVisible(True)
