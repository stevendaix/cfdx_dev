"""Semantic visual tokens and Qt application for the CFDX Scientific Workbench.

This module deliberately has no Qt import at module import time: the GUI remains
an optional dependency for core CFDX users and command-line workflows.
"""

from __future__ import annotations

# Values are presentation tokens, not scientific colormaps or solver state.
TOKENS: dict[str, str] = {
    "surface.canvas": "#10151B",
    "surface.panel": "#171E26",
    "surface.raised": "#202A35",
    "surface.input": "#0D1217",
    "border.subtle": "#2B3744",
    "text.primary": "#E7EDF3",
    "text.secondary": "#A7B4C2",
    "text.muted": "#758495",
    "accent.primary": "#39C6D8",
    "accent.secondary": "#65D6B0",
    "status.success": "#65C995",
    "status.warning": "#E8B45B",
    "status.error": "#F07878",
    "status.info": "#7DAEFF",
}

WORKBENCH_STYLESHEET = """
QMainWindow, QWidget {
    background-color: %(surface.canvas)s;
    color: %(text.primary)s;
}
QToolBar {
    background-color: %(surface.panel)s;
    border: 0;
    border-bottom: 1px solid %(border.subtle)s;
    spacing: 6px;
    padding: 5px;
}
QToolBar QToolButton {
    background: transparent;
    color: %(text.primary)s;
    border: 1px solid transparent;
    border-radius: 5px;
    padding: 6px 9px;
}
QToolBar QToolButton:hover {
    background: %(surface.raised)s;
    border-color: %(border.subtle)s;
}
QToolBar QToolButton:checked {
    background: %(surface.raised)s;
    border-color: %(accent.primary)s;
}
QDockWidget {
    color: %(text.primary)s;
    font-weight: 600;
}
QDockWidget::title {
    background: %(surface.panel)s;
    border-bottom: 1px solid %(border.subtle)s;
    padding: 7px 9px;
}
QStatusBar {
    background: %(surface.panel)s;
    color: %(text.secondary)s;
    border-top: 1px solid %(border.subtle)s;
}
QTreeWidget, QTableView, QTableWidget, QPlainTextEdit, QTextEdit, QLineEdit,
QComboBox, QSpinBox, QDoubleSpinBox {
    background: %(surface.input)s;
    color: %(text.primary)s;
    border: 1px solid %(border.subtle)s;
    border-radius: 4px;
    selection-background-color: %(accent.primary)s;
    selection-color: %(surface.canvas)s;
}
QTreeWidget::item, QTableView::item, QTableWidget::item {
    padding: 4px;
}
QPushButton {
    background: %(surface.raised)s;
    color: %(text.primary)s;
    border: 1px solid %(border.subtle)s;
    border-radius: 5px;
    padding: 7px 12px;
}
QPushButton:hover {
    border-color: %(accent.primary)s;
}
QPushButton:focus, QLineEdit:focus, QComboBox:focus, QSpinBox:focus,
QDoubleSpinBox:focus, QTreeWidget:focus, QTableView:focus {
    border: 1px solid %(accent.primary)s;
}
QPushButton:disabled {
    color: %(text.muted)s;
    background: %(surface.panel)s;
}
QPushButton[emphasis="primary"] {
    background: %(accent.primary)s;
    color: %(surface.canvas)s;
    border-color: %(accent.primary)s;
    font-weight: 600;
}
QPushButton[emphasis="primary"]:hover {
    background: %(accent.secondary)s;
    border-color: %(accent.secondary)s;
}
QLabel[role="secondary"] {
    color: %(text.secondary)s;
}
QMenu {
    background: %(surface.raised)s;
    color: %(text.primary)s;
    border: 1px solid %(border.subtle)s;
}
QMenu::item:selected {
    background: %(accent.primary)s;
    color: %(surface.canvas)s;
}
QToolTip {
    background: %(surface.raised)s;
    color: %(text.primary)s;
    border: 1px solid %(border.subtle)s;
    padding: 4px;
}
""" % TOKENS


def apply_workbench_theme(window: object) -> None:
    """Apply the semantic palette and shared shell stylesheet to a Qt window.

    Qt is imported lazily so importing this module never makes PySide6 a
    mandatory dependency. The caller is expected to pass a QWidget/QMainWindow.
    """
    try:
        from PySide6.QtGui import QColor, QPalette
    except ImportError as exc:  # pragma: no cover - depends on optional extra
        raise RuntimeError("PySide6 is required to apply the CFDX GUI theme") from exc

    palette = QPalette()
    roles = {
        QPalette.ColorRole.Window: "surface.canvas",
        QPalette.ColorRole.WindowText: "text.primary",
        QPalette.ColorRole.Base: "surface.input",
        QPalette.ColorRole.AlternateBase: "surface.panel",
        QPalette.ColorRole.Text: "text.primary",
        QPalette.ColorRole.Button: "surface.raised",
        QPalette.ColorRole.ButtonText: "text.primary",
        QPalette.ColorRole.Highlight: "accent.primary",
        QPalette.ColorRole.HighlightedText: "surface.canvas",
        QPalette.ColorRole.ToolTipBase: "surface.raised",
        QPalette.ColorRole.ToolTipText: "text.primary",
        QPalette.ColorRole.PlaceholderText: "text.muted",
    }
    for role, token in roles.items():
        palette.setColor(role, QColor(TOKENS[token]))

    window.setPalette(palette)
    window.setStyleSheet(WORKBENCH_STYLESHEET)
