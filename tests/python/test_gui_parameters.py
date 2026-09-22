import importlib.util
import pytest

from cfdx import CFDXSession
from cfdx.gui import CFDXMainWindow


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_parameter_form_updates_session() -> None:
    from cfdx.gui import create_application
    app = create_application(["cfdx-test"])
    session = CFDXSession()
    window = CFDXMainWindow(session)
    window.cfl.setValue(12.5)
    assert session.case.numerics["cfl"] == 12.5
    window.close()
    app.quit()
