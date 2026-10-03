import importlib.util
import os

import pytest


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None
    or importlib.util.find_spec("pyvistaqt") is None,
    reason="optional Qt/PyVista stack",
)
@pytest.mark.skipif(
    os.environ.get("CFDX_RUN_PYVISTA_QT_SMOKE") != "1",
    reason="PyVistaQt OpenGL smoke test is opt-in and requires a stable display",
)
def test_pyvista_qt_view_smoke() -> None:
    from cfdx.gui_3d import PyVistaQtView, QWidget
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    view = PyVistaQtView()
    assert view.plotter is not None
    view.close()
    app.processEvents()
