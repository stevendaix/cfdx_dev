import importlib.util
import pytest

from cfdx import CFDXSession
from cfdx.gui import CFDXMainWindow


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_controls_drive_shared_session() -> None:
    from cfdx.gui import create_application
    app = create_application(["cfdx-test"])
    session = CFDXSession()
    window = CFDXMainWindow(session)
    window.run_button.click()
    assert session.state.value == "RUNNING"
    window.pause_button.click()
    assert session.state.value == "PAUSED"
    window.run_button.click()
    window.stop_button.click()
    assert session.state.value == "STOPPED"
    window.close()
    app.quit()
