import importlib.util

import pytest
from cfdx.theme import TOKENS, WORKBENCH_STYLESHEET, apply_workbench_theme


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
    assert TOKENS["surface.canvas"] == "#10151B"
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
        == TOKENS["surface.canvas"]
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
@pytest.mark.parametrize("surface_token", ["surface.canvas", "surface.panel", "surface.input"])
def test_text_token_contrast_meets_wcag_aa_for_normal_text(
    text_token: str, surface_token: str
) -> None:
    assert _contrast_ratio(TOKENS[text_token], TOKENS[surface_token]) >= 4.5, (
        f"{text_token} on {surface_token} does not meet the 4.5:1 WCAG AA "
        "contrast threshold for normal text"
    )
