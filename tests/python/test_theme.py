import importlib.util

import pytest
from cfdx.theme import (
    DEFAULT_THEME, THEMES, TOKENS, WORKBENCH_STYLESHEET,
    apply_workbench_theme, theme_tokens, workbench_stylesheet, style_scientific_plot,
)


def test_theme_tokens_are_semantic_and_complete() -> None:
    expected = {
        "surface.canvas",
        "surface.panel",
        "surface.raised",
        "surface.input",
        "border.subtle",
        "text.primary",
        "text.secondary",
        "text.muted",
        "accent.primary",
        "accent.secondary",
        "status.success",
        "status.warning",
        "status.error",
        "status.info",
    }
    assert set(TOKENS) == expected
    assert all(value.startswith("#") and len(value) == 7 for value in TOKENS.values())
    assert "%(" not in WORKBENCH_STYLESHEET
    assert TOKENS["accent.primary"] in WORKBENCH_STYLESHEET
    assert TOKENS["status.error"] in WORKBENCH_STYLESHEET
    assert 'QLabel[validationState="error"]' in WORKBENCH_STYLESHEET
    assert '*[validationState="warning"]' in WORKBENCH_STYLESHEET


def test_theme_module_imports_without_qt() -> None:
    # Importing semantic tokens must remain safe for base installations.
    assert TOKENS["surface.canvas"] == "#000000"
    if importlib.util.find_spec("PySide6") is None:
        with pytest.raises(RuntimeError, match="PySide6 is required"):
            apply_workbench_theme(object())


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_theme_application_sets_palette_and_stylesheet() -> None:
    from PySide6.QtGui import QPalette
    from PySide6.QtWidgets import QApplication, QMainWindow

    app = QApplication.instance() or QApplication(["cfdx-theme-test"])
    window = QMainWindow()
    apply_workbench_theme(window)
    assert window.styleSheet() == WORKBENCH_STYLESHEET
    assert (
        window.palette().color(QPalette.ColorRole.Window).name()
        == TOKENS["surface.canvas"].lower()
    )
    window.close()
    app.quit()


def _relative_luminance(hex_color: str) -> float:
    channels = [int(hex_color[index : index + 2], 16) / 255 for index in (1, 3, 5)]
    linear = [
        channel / 12.92 if channel <= 0.04045 else ((channel + 0.055) / 1.055) ** 2.4
        for channel in channels
    ]
    return 0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2]


def _contrast_ratio(foreground: str, background: str) -> float:
    luminances = sorted(
        (_relative_luminance(foreground), _relative_luminance(background)),
        reverse=True,
    )
    return (luminances[0] + 0.05) / (luminances[1] + 0.05)


@pytest.mark.parametrize("text_token", ["text.primary", "text.secondary", "text.muted"])
@pytest.mark.parametrize(
    "surface_token", ["surface.canvas", "surface.panel", "surface.input"]
)
def test_text_token_contrast_meets_wcag_aa_for_normal_text(
    text_token: str, surface_token: str
) -> None:
    assert _contrast_ratio(TOKENS[text_token], TOKENS[surface_token]) >= 4.5, (
        f"{text_token} on {surface_token} does not meet the 4.5:1 WCAG AA "
        "contrast threshold for normal text"
    )


@pytest.mark.parametrize(
    "status_token",
    ["status.success", "status.warning", "status.error", "status.info"],
)
def test_status_colors_have_nontext_contrast_on_panel(status_token: str) -> None:
    # WCAG 2.2 non-text contrast: meaningful status indicators need 3:1.
    assert _contrast_ratio(TOKENS[status_token], TOKENS["surface.panel"]) >= 3.0, (
        f"{status_token} on surface.panel does not meet the 3:1 non-text "
        "contrast threshold"
    )


@pytest.mark.parametrize(
    "surface_token", ["surface.canvas", "surface.panel", "surface.input"]
)
def test_focus_indicator_has_nontext_contrast(surface_token: str) -> None:
    # The primary accent is used for keyboard focus borders and must remain visible.
    assert _contrast_ratio(TOKENS["accent.primary"], TOKENS[surface_token]) >= 3.0, (
        f"accent.primary on {surface_token} does not meet the 3:1 non-text "
        "contrast threshold"
    )


def test_focus_styles_cover_multiline_editors_and_table_widgets() -> None:
    focus_rule = WORKBENCH_STYLESHEET.split("QPushButton:focus", maxsplit=1)[1].split(
        "{", maxsplit=1
    )[0]
    for selector in ("QPlainTextEdit:focus", "QTextEdit:focus", "QTableWidget:focus"):
        assert selector in focus_rule
    assert "border: 1px solid " + TOKENS["accent.primary"] in WORKBENCH_STYLESHEET


@pytest.mark.parametrize(
    "status_token",
    ["status.success", "status.warning", "status.error", "status.info"],
)
@pytest.mark.parametrize("surface_token", ["surface.canvas", "surface.panel"])
def test_status_text_meets_wcag_aa_for_normal_text(
    status_token: str, surface_token: str
) -> None:
    # Validation labels use status colors as text, so the stricter text threshold applies.
    assert _contrast_ratio(TOKENS[status_token], TOKENS[surface_token]) >= 4.5, (
        f"{status_token} on {surface_token} does not meet the 4.5:1 WCAG AA "
        "contrast threshold for normal text"
    )


@pytest.mark.parametrize(
    "surface_token", ["surface.canvas", "surface.panel", "surface.input"]
)
def test_semantic_border_has_nontext_contrast(surface_token: str) -> None:
    # Shared control boundaries must remain distinguishable from adjacent surfaces.
    assert _contrast_ratio(TOKENS["border.subtle"], TOKENS[surface_token]) >= 3.0, (
        f"border.subtle on {surface_token} does not meet the 3:1 non-text "
        "contrast threshold"
    )


def test_dense_data_views_use_semantic_table_chrome() -> None:
    # Keep data navigation readable without changing scientific field colormaps.
    alternate_rows = "alternate-background-color: " + TOKENS["surface.panel"]
    assert alternate_rows in WORKBENCH_STYLESHEET
    assert "QHeaderView::section" in WORKBENCH_STYLESHEET
    assert "color: " + TOKENS["text.secondary"] in WORKBENCH_STYLESHEET
    assert "QTreeView::item:hover" in WORKBENCH_STYLESHEET
    assert "QTableView::item:hover" in WORKBENCH_STYLESHEET
    assert "QAbstractItemView:disabled" in WORKBENCH_STYLESHEET
    assert "show-decoration-selected: 1" in WORKBENCH_STYLESHEET
    assert "QTreeView:focus" in WORKBENCH_STYLESHEET


def test_black_is_the_first_run_default_and_compiled_styles_are_complete() -> None:
    assert DEFAULT_THEME == "black"
    assert theme_tokens() is TOKENS
    assert TOKENS["surface.canvas"] == "#000000"
    for name in ("black", "graphite", "light"):
        assert set(THEMES[name]) == set(TOKENS)
        css = workbench_stylesheet(name)
        assert "@surface." not in css
        assert "@status." not in css
        assert "QTabWidget::pane" in css
        assert "QScrollBar:vertical" in css
        assert THEMES[name]["surface.canvas"] in css
        assert "QLabel[role=\"eyebrow\"]" in css


def test_theme_rejects_unsupported_choices_before_mutating_window() -> None:
    with pytest.raises(ValueError, match="unsupported CFDX theme"):
        theme_tokens("fluorescent")
    with pytest.raises(ValueError, match="unsupported CFDX theme"):
        apply_workbench_theme(object(), "unexpected")


@pytest.mark.parametrize("name", list(THEMES))
@pytest.mark.parametrize("text_token", ["text.primary", "text.secondary", "text.muted"])
@pytest.mark.parametrize("surface", ["surface.canvas", "surface.panel", "surface.input", "surface.raised"])
def test_all_theme_text_contrasts(name: str, text_token: str, surface: str) -> None:
    colors = THEMES[name]
    assert _contrast_ratio(colors[text_token], colors[surface]) >= 4.5


@pytest.mark.parametrize("name", list(THEMES))
@pytest.mark.parametrize("highlight", ["border.subtle", "accent.primary"])
@pytest.mark.parametrize("surface", ["surface.canvas", "surface.panel", "surface.input", "surface.raised"])
def test_all_theme_control_indicators_are_legible(name: str, highlight: str, surface: str) -> None:
    colors = THEMES[name]
    assert _contrast_ratio(colors[highlight], colors[surface]) >= 3.0


@pytest.mark.parametrize("name", list(THEMES))
@pytest.mark.parametrize("status", ["status.success", "status.warning", "status.error", "status.info"])
def test_theme_status_labels_are_legible(name: str, status: str) -> None:
    colors = THEMES[name]
    assert _contrast_ratio(colors[status], colors["surface.panel"]) >= 4.5


def test_plot_styling_changes_only_background_and_axis_chrome() -> None:
    class Axis:
        def setPen(self, color): self.pen = color
        def setTextPen(self, color): self.text = color
    class FakePlot:
        def __init__(self): self.axes = {k: Axis() for k in ("left", "bottom")}
        def setBackground(self, color): self.background = color
        def getPlotItem(self): return self
        def getAxis(self, key): return self.axes[key]
    plot = FakePlot()
    style_scientific_plot(plot, "light")
    assert plot.background == THEMES["light"]["surface.canvas"]
    assert plot.axes["left"].text == THEMES["light"]["text.secondary"]
    assert plot.axes["bottom"].pen == THEMES["light"]["border.subtle"]
