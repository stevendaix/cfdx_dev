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
    "text.muted": "#7A8999",
    "accent.primary": "#39C6D8",
    "accent.secondary": "#65D6B0",
    "status.success": "#65C995",
    "status.warning": "#E8B45B",
    "status.error": "#F07878",
    "status.info": "#7DAEFF",
}

WORKBENCH_STYLESHEET = """
QMainWindow, QWidget {
    background-color: @surface.canvas@;
    color: @text.primary@;
}
QToolBar {
    background-color: @surface.panel@;
    border: 0;
    border-bottom: 1px solid @border.subtle@;
    spacing: 6px;
    padding: 5px;
}
QToolBar QToolButton {
    background: transparent;
    color: @text.primary@;
    border: 1px solid transparent;
    border-radius: 5px;
    padding: 6px 9px;
}
QToolBar QToolButton:hover {
    background: @surface.raised@;
    border-color: @border.subtle@;
}
QToolBar QToolButton:checked {
    background: @surface.raised@;
    border-color: @accent.primary@;
}
QDockWidget {
    color: @text.primary@;
    font-weight: 600;
}
QDockWidget::title {
    background: @surface.panel@;
    border-bottom: 1px solid @border.subtle@;
    padding: 7px 9px;
}
QStatusBar {
    background: @surface.panel@;
    color: @text.secondary@;
    border-top: 1px solid @border.subtle@;
}
QTreeWidget, QTableView, QTableWidget, QPlainTextEdit, QTextEdit, QLineEdit,
QComboBox, QSpinBox, QDoubleSpinBox {
    background: @surface.input@;
    color: @text.primary@;
    border: 1px solid @border.subtle@;
    border-radius: 4px;
    selection-background-color: @accent.primary@;
    selection-color: @surface.canvas@;
}
QTreeWidget::item, QTableView::item, QTableWidget::item {
    padding: 4px;
}
QPushButton {
    background: @surface.raised@;
    color: @text.primary@;
    border: 1px solid @border.subtle@;
    border-radius: 5px;
    padding: 7px 12px;
}
QPushButton:hover {
    border-color: @accent.primary@;
}
QPushButton:focus, QLineEdit:focus, QComboBox:focus, QSpinBox:focus,
QDoubleSpinBox:focus, QTreeWidget:focus, QTableView:focus {
    border: 1px solid @accent.primary@;
}
QPushButton:disabled {
    color: @text.muted@;
    background: @surface.panel@;
}
QPushButton[emphasis="primary"] {
    background: @accent.primary@;
    color: @surface.canvas@;
    border-color: @accent.primary@;
    font-weight: 600;
}
QPushButton[emphasis="primary"]:hover {
    background: @accent.secondary@;
    border-color: @accent.secondary@;
}
QLabel[role="secondary"] {
    color: @text.secondary@;
}
QLabel[validationState="error"] {
    color: @status.error@;
}
QLabel[validationState="warning"] {
    color: @status.warning@;
}
*[validationState="error"] {
    border: 1px solid @status.error@;
}
*[validationState="warning"] {
    border: 1px solid @status.warning@;
}
QMenu {
    background: @surface.raised@;
    color: @text.primary@;
    border: 1px solid @border.subtle@;
}
QMenu::item:selected {
    background: @accent.primary@;
    color: @surface.canvas@;
}
QToolTip {
    background: @surface.raised@;
    color: @text.primary@;
    border: 1px solid @border.subtle@;
    padding: 4px;
}
"""
for _token_name, _token_value in TOKENS.items():
    WORKBENCH_STYLESHEET = WORKBENCH_STYLESHEET.replace(
        f"@{_token_name}@", _token_value
    )


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
