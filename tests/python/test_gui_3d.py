import importlib.util

import pytest


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None
    or importlib.util.find_spec("pyvistaqt") is None,
    reason="optional Qt/PyVista stack",
)
def test_pyvista_qt_view_smoke() -> None:
    from cfdx.gui_3d import PyVistaQtView, QWidget

    view = PyVistaQtView()
    assert view.plotter is not None
    view.close()
