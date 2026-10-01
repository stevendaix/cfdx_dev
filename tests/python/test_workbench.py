import importlib.util

import pytest

from cfdx import CFDXSession
from cfdx.workbench import CFDXWorkbenchWindow


def test_workbench_is_importable_without_qt() -> None:
    assert CFDXWorkbenchWindow is not None
    if importlib.util.find_spec("PySide6") is not None:
        pytest.skip("PySide6 is installed")
    with pytest.raises(RuntimeError, match="PySide6 is required"):
        CFDXWorkbenchWindow(CFDXSession())


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_workbench_has_stable_docks() -> None:
    from PySide6.QtWidgets import QDockWidget
    from cfdx.gui import create_application
    from cfdx.session import SimulationState

    app = create_application(["cfdx-workbench-test"])
    session = CFDXSession()
    session.state = SimulationState.READY
    window = CFDXWorkbenchWindow(session)
    assert window.objectName() == ""
    assert window.workflow_tree.objectName() == "workbench.workflow_tree"
    assert window.workflow_tree.topLevelItem(1).child(2).text(0).endswith("[OK]")
    assert window.workflow_tree.topLevelItem(0).child(4).text(0).endswith("[ ]")
    assert window.findChild(QDockWidget, "workbench.dock.workflow") is not None
    assert window.findChild(QDockWidget, "workbench.dock.properties") is not None
    assert window.findChild(QDockWidget, "workbench.dock.monitor") is not None
    assert window.centralWidget().objectName() == "workbench.viewport"
    window.workflow_tree.setCurrentItem(window.workflow_tree.topLevelItem(0).child(4))
    assert window.application.state.selection.stable_id == "boundaries"
    window.close()
    app.quit()
