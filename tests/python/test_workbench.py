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

    app = create_application(["cfdx-workbench-test"])
    window = CFDXWorkbenchWindow(CFDXSession())
    assert window.objectName() == ""
    assert window.workflow_tree.objectName() == "workbench.workflow_tree"
    assert window.findChild(QDockWidget, "workbench.dock.workflow") is not None
    assert window.findChild(QDockWidget, "workbench.dock.properties") is not None
    assert window.findChild(QDockWidget, "workbench.dock.monitor") is not None
    assert window.centralWidget().objectName() == "workbench.viewport"
    window.close()
    app.quit()
