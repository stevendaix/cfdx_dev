import importlib.util
import sys
from pathlib import Path

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
    window.setup_panel.physics_model.setCurrentText("energy")
    window.setup_panel.physics_enabled.setChecked(True)
    window.setup_panel.physics_button.click()
    assert session.case.physics["energy"]["enabled"] is True
    window.setup_panel.boundary_name.setText("inlet")
    window.setup_panel.boundary_type.setCurrentText("inlet")
    window.setup_panel.boundary_value.setText("1.0")
    window.setup_panel.boundary_button.click()
    assert session.case.boundaries["inlet"] == {"type": "inlet", "value": "1.0"}
    window.run_button.click()
    assert session.state is SimulationState.RUNNING
    window.pause_button.click()
    assert session.state is SimulationState.PAUSED
    window.stop_button.click()
    assert session.state is SimulationState.STOPPED
    window.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_file_actions_roundtrip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from PySide6.QtWidgets import QFileDialog
    from cfdx.case_io import save_case_with_dat
    from cfdx.gui import create_application

    app = create_application(["cfdx-file-test"])
    session = CFDXSession()
    session.case.name = "channel"
    session.case.execution.solver = sys.executable
    window = CFDXMainWindow(session)

    case_path = tmp_path / "channel.cfdx.h5"
    dat_source = tmp_path / "solver.dat"
    dat_source.write_text("restart\n", encoding="utf-8")

    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        staticmethod(lambda *args, **kwargs: (str(case_path), "CFDX Case (*.cfdx.h5 *.h5)")),
    )
    assert window._save_case_as()
    assert case_path.is_file()

    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *args, **kwargs: (str(dat_source), "Solver restart (*.dat)")),
    )
    assert window._save_case_with_dat()
    assert (tmp_path / "channel.dat").is_file()

    window.close()
    app.quit()
