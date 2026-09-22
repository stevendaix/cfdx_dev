from __future__ import annotations

import importlib.util

import pytest

from cfdx import CFDXSession
from cfdx.gui import CFDXMainWindow, launch


def test_gui_module_is_importable() -> None:
    assert CFDXMainWindow is not None


def test_gui_without_qt_has_clear_error() -> None:
    if importlib.util.find_spec("PySide6") is not None:
        pytest.skip("optional GUI dependency is installed")
    with pytest.raises(RuntimeError, match="PySide6 is required"):
        CFDXMainWindow(CFDXSession())
    with pytest.raises(RuntimeError, match="PySide6 is required"):
        launch(CFDXSession(), [])


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_model_tree_and_tui() -> None:
    from PySide6.QtCore import Qt
    from cfdx.gui import create_application

    app = create_application(["cfdx-test"])
    window = CFDXMainWindow(CFDXSession())
    nodes = window.session.case_tree()
    assert window.tree.topLevelItemCount() == len(nodes)
    assert window.tree.topLevelItem(0).data(0, Qt.ItemDataRole.UserRole) == nodes[0].id
    assert "State: CREATED" in window.render_tui()
    window.close()
    app.quit()
