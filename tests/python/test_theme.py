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
    assert window.palette().color(QPalette.ColorRole.Window).name() == TOKENS[
        "surface.canvas"
    ]
    window.close()
    app.quit()
