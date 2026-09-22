import importlib.util
import pytest

from cfdx import CFDXSession, SimulationState
from cfdx.gui import CFDXMainWindow


def test_gui_module_is_importable() -> None:
    assert CFDXMainWindow is not None


def test_gui_without_qt_has_clear_error() -> None:
    if importlib.util.find_spec("PySide6") is not None:
        pytest.skip("PySide6 is installed")
    with pytest.raises(RuntimeError, match="PySide6 is required"):
        CFDXMainWindow(CFDXSession())


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_controls_and_parameters() -> None:
    from PySide6.QtCore import Qt
    from cfdx.gui import create_application
    app = create_application(["cfdx-test"])
    session = CFDXSession()
    window = CFDXMainWindow(session)
    nodes = session.case_tree()
    assert window.tree.topLevelItemCount() == len(nodes)
    assert window.tree.topLevelItem(0).data(0, Qt.ItemDataRole.UserRole) == nodes[0].id
    window.cfl.setValue(12.5)
    assert session.case.numerics["cfl"] == 12.5
    window.run_button.click()
    assert session.state is SimulationState.RUNNING
    window.pause_button.click()
    assert session.state is SimulationState.PAUSED
    window.stop_button.click()
    assert session.state is SimulationState.STOPPED
    window.close()
    app.quit()
