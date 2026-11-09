"""Semantic visual tokens and Qt application for the CFDX Scientific Workbench.

This module deliberately has no Qt import at module import time: the GUI remains
an optional dependency for core CFDX users and command-line workflows.
"""

from __future__ import annotations

# Presentation colors are not scientific color maps or numerical settings.
DEFAULT_THEME = "black"
THEME_LABELS = {"black": "OLED Black", "graphite": "Graphite", "light": "Light"}
THEMES: dict[str, dict[str, str]] = {
    "black": {
        "surface.canvas": "#000000",
        "surface.panel": "#101114",
        "surface.raised": "#1B1E24",
        "surface.input": "#08090C",
        "border.subtle": "#707985",
        "text.primary": "#F2F5FA",
        "text.secondary": "#B7C2D0",
        "text.muted": "#9AAAB9",
        "accent.primary": "#79DFFC",
        "accent.secondary": "#91E4BF",
        "status.success": "#7BDDB2",
        "status.warning": "#F2C174",
        "status.error": "#FF9898",
        "status.info": "#A1BEFF",
    },
    "graphite": {
        "surface.canvas": "#10151B",
        "surface.panel": "#171E26",
        "surface.raised": "#202A35",
        "surface.input": "#0D1217",
        "border.subtle": "#748293",
        "text.primary": "#E7EDF3",
        "text.secondary": "#BCC8D4",
        "text.muted": "#A3B0BF",
        "accent.primary": "#65D9EA",
        "accent.secondary": "#77E2B8",
        "status.success": "#74D7A5",
        "status.warning": "#F4C57F",
        "status.error": "#F9A3A3",
        "status.info": "#9DBFFF",
    },
    "light": {
        "surface.canvas": "#F5F7FA",
        "surface.panel": "#FFFFFF",
        "surface.raised": "#E8EDF3",
        "surface.input": "#FFFFFF",
        "border.subtle": "#5F6F80",
        "text.primary": "#142334",
        "text.secondary": "#35485E",
        "text.muted": "#4D6175",
        "accent.primary": "#006E85",
        "accent.secondary": "#166F4B",
        "status.success": "#136B43",
        "status.warning": "#875400",
        "status.error": "#A12E36",
        "status.info": "#2457A3",
    },
}
TOKENS: dict[str, str] = THEMES[DEFAULT_THEME]


def theme_tokens(theme: str = DEFAULT_THEME) -> dict[str, str]:
    """Return tokens for an explicitly supported theme (never silent fallback)."""
    if theme not in THEMES:
        raise ValueError(f"unsupported CFDX theme: {theme!r}")
    return THEMES[theme]


_STYLE_TEMPLATE = """
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
QAbstractItemView {
    alternate-background-color: @surface.panel@;
    show-decoration-selected: 1;
    outline: 0;
}
QTreeView::item, QTableView::item, QTableWidget::item {
    padding: 4px 6px;
}
QTreeView::item:hover, QTableView::item:hover, QTableWidget::item:hover {
    background: @surface.raised@;
}
QHeaderView::section {
    background: @surface.panel@;
    color: @text.secondary@;
    border: 0;
    border-right: 1px solid @border.subtle@;
    border-bottom: 1px solid @border.subtle@;
    padding: 6px 8px;
    font-weight: 600;
}
QHeaderView::section:checked {
    color: @text.primary@;
    background: @surface.raised@;
}
QAbstractItemView:disabled {
    color: @text.muted@;
    background: @surface.panel@;
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
QDoubleSpinBox:focus, QPlainTextEdit:focus, QTextEdit:focus,
QTreeWidget:focus, QTreeView:focus, QTableView:focus, QTableWidget:focus {
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
QMenuBar {
    background: @surface.panel@;
    color: @text.primary@;
    border-bottom: 1px solid @border.subtle@;
    padding: 3px 8px;
}
QMenuBar::item { padding: 6px 11px; border-radius: 6px; }
QMenuBar::item:selected { background: @surface.raised@; }
QTabWidget::pane { border: 1px solid @border.subtle@; border-radius: 6px; }
QTabBar::tab {
    background: @surface.panel@;
    color: @text.secondary@;
    border: 1px solid @border.subtle@;
    border-bottom: 0;
    padding: 9px 16px;
    margin-right: 3px;
    min-width: 68px;
}
QTabBar::tab:selected {
    background: @surface.raised@;
    color: @text.primary@;
    border-top: 2px solid @accent.primary@;
}
QTabBar::tab:hover:!selected { background: @surface.raised@; }
QGroupBox {
    border: 1px solid @border.subtle@;
    border-radius: 7px;
    margin-top: 12px;
    padding: 12px 9px 8px 9px;
    font-weight: 600;
}
QGroupBox::title { subcontrol-origin: margin; left: 12px; padding: 0 5px; }
QComboBox { padding: 7px 23px 7px 9px; }
QComboBox QAbstractItemView {
    background: @surface.raised@;
    color: @text.primary@;
    selection-background-color: @accent.primary@;
    selection-color: @surface.canvas@;
}
QCheckBox, QRadioButton { background: transparent; spacing: 8px; }
QCheckBox::indicator, QRadioButton::indicator {
    width: 15px;
    height: 15px;
}
QScrollBar:vertical, QScrollBar:horizontal {
    background: @surface.input@;
    border: 0;
    width: 12px;
    height: 12px;
}
QScrollBar::handle:vertical, QScrollBar::handle:horizontal {
    background: @border.subtle@;
    border: 3px solid @surface.input@;
    border-radius: 5px;
    min-height: 30px;
    min-width: 30px;
}
QScrollBar::add-line, QScrollBar::sub-line {
    width: 0px;
    height: 0px;
}
QSplitter::handle { background: @border.subtle@; }
QProgressBar {
    color: @text.primary@;
    background: @surface.input@;
    border: 1px solid @border.subtle@;
    border-radius: 5px;
    text-align: center;
}
QProgressBar::chunk { background: @accent.primary@; border-radius: 4px; }
QLabel[role="eyebrow"] {
    color: @accent.primary@;
    font-weight: 600;
}
QLabel[role="muted"] { color: @text.muted@; }
QFrame[role="dashboardCard"] {
    background: @surface.panel@;
    border: 1px solid @border.subtle@;
    border-radius: 10px;
}
QLabel[role="dashboardValue"] { font-weight: 700; color: @text.primary@; }
QToolBar QToolButton[workbenchRun="true"] {
    background: @accent.primary@;
    color: @surface.canvas@;
    border: 1px solid @accent.primary@;
    border-radius: 6px;
    font-weight: 700;
}
QToolBar QToolButton[workbenchRun="true"]:disabled {
    background: @surface.panel@;
    color: @text.muted@;
    border: 1px solid @border.subtle@;
}
QToolBar QToolButton[workbenchNav="true"] {
    padding: 8px 17px;
    margin: 2px 3px;
    border-radius: 6px;
}
QToolBar QToolButton[workbenchNav="true"]:checked {
    color: @accent.primary@;
    border: 1px solid @accent.primary@;
    background: @surface.raised@;
}
QToolBar QToolButton[workbenchTheme="true"] {
    border: 1px solid @border.subtle@;
    border-radius: 6px;
    padding: 7px 11px;
}
"""


def workbench_stylesheet(theme: str = DEFAULT_THEME) -> str:
    """Compile semantic stylesheet without depending on a Qt installation."""
    result = _STYLE_TEMPLATE
    for name, value in theme_tokens(theme).items():
        result = result.replace(f"@{name}@", value)
    return result


WORKBENCH_STYLESHEET = workbench_stylesheet()


def apply_workbench_theme(window: object, theme: str = DEFAULT_THEME) -> None:
    """Apply one named palette, including disabled state, to an existing Qt widget.

    A missing optional PySide6 installation does not affect headless CFDX.
    Unsupported themes are rejected before mutating window state.
    """
    tokens = theme_tokens(theme)
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
        palette.setColor(role, QColor(tokens[token]))
    for role in (
        QPalette.ColorRole.Text,
        QPalette.ColorRole.WindowText,
        QPalette.ColorRole.ButtonText,
    ):
        palette.setColor(
            QPalette.ColorGroup.Disabled, role, QColor(tokens["text.muted"])
        )

    window.setPalette(palette)
    window.setStyleSheet(workbench_stylesheet(theme))
    window.setProperty("cfdxTheme", theme)


def style_scientific_plot(plot: object, theme: str = DEFAULT_THEME) -> None:
    """Style pyqtgraph axes/background only; do not touch data or scalings."""
    tokens = theme_tokens(theme)
    plot.setBackground(tokens["surface.canvas"])
    axes = plot.getPlotItem()
    for direction in ("left", "bottom"):
        axis = axes.getAxis(direction)
        axis.setPen(tokens["border.subtle"])
        axis.setTextPen(tokens["text.secondary"])
