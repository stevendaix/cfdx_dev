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
def test_gui_controls_and_parameters(tmp_path: Path) -> None:
    from PySide6.QtCore import Qt
    from cfdx.gui import create_application
    app = create_application(["cfdx-test"])
    session = CFDXSession()
    window = CFDXMainWindow(session)
    nodes = session.case_tree()
    assert window.tree.topLevelItemCount() == len(nodes)
    assert window.open_result_button is not None
    assert window.tree.topLevelItem(0).data(0, Qt.ItemDataRole.UserRole) == nodes[0].id
    window.cfl.setValue(12.5)
    assert session.case.numerics["cfl"] == 12.5
    window.setup_panel.physics_model.setCurrentText("incompressible")
    window.setup_panel.physics_enabled.setChecked(True)
    window.setup_panel.physics_button.click()
    window.setup_panel.physics_model.setCurrentText("energy")
    window.setup_panel.physics_enabled.setChecked(True)
    window.setup_panel.physics_button.click()
    assert session.case.physics["energy"]["enabled"] is True
    window.setup_panel.material_name.setText("water")
    window.setup_panel.material_button.click()
    window.setup_panel.boundary_name.setText("inlet")
    window.setup_panel.boundary_type.setCurrentText("inlet")
    window.setup_panel.boundary_value.setText("1.0")
    window.setup_panel.boundary_button.click()
    assert session.case.boundaries["inlet"] == {"type": "inlet", "value": "1.0"}
    session.case.boundaries["inlet"]["vendor_extension"] = {"keep": True}
    window.setup_panel.boundary_value.setValue(2.0)
    window.setup_panel.boundary_button.click()
    assert session.case.boundaries["inlet"]["vendor_extension"] == {"keep": True}
    class FakeRunner:
        running = False

    class FakeController:
        def __init__(self, session):
            self.session = session
            self.runner = FakeRunner()
            self.on_output = None
            self.on_metrics = None
            self.on_complete = None

        def start(self):
            self.session.run()

        def pause(self):
            self.session.pause()

        def resume(self):
            self.session.run()

        def stop(self):
            self.session.stop()

    window.controller = FakeController(session)
    window._case_path = tmp_path / "channel.cfdx.h5"
    window._dirty = False
    window.run_button.click()
    assert session.state is SimulationState.RUNNING
    window.pause_button.click()
    assert session.state is SimulationState.PAUSED
    window.resume_button.click()
    assert session.state is SimulationState.RUNNING
    window.stop_button.click()
    assert session.state is SimulationState.STOPPED
    window._dirty = False
    window.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_file_actions_roundtrip(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from PySide6.QtWidgets import QFileDialog
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
        staticmethod(lambda *args, **kwargs: (str(case_path), "CFDX Case (*.cfdx.h5)")),
    )
    assert window._save_case_as()
    assert case_path.is_file()

    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *args, **kwargs: (str(dat_source), "Solver restart (*.dat)")),
    )
    assert window._save_case_with_dat()
    assert (tmp_path / "channel.dat.h5").is_file()

    window.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_real_execution_controller_roundtrip(tmp_path: Path) -> None:
    from cfdx.gui import create_application

    solver = tmp_path / "solver.py"
    solver.write_text(
        "import sys\n"
        "print('Iteration 3 Time = 0.5 CFL: 0.4', flush=True)\n"
        "assert sys.argv[1].endswith('.cfdx.h5')\n",
        encoding="utf-8",
    )
    app = create_application(["cfdx-e2e-test"])
    session = CFDXSession()
    session.case.name = "e2e"
    session.case.execution.solver = sys.executable
    window = CFDXMainWindow(session)
    window._case_path = tmp_path / "e2e.cfdx.h5"
    from cfdx.case_io import save_case
    save_case(session, window._case_path)
    window._dirty = False

    window.run_button.click()
    assert window.controller is not None
    assert window.controller.runner._thread is not None
    window.controller.runner._thread.join(timeout=5)

    assert session.state is SimulationState.CONVERGED
    assert session.iteration == 3
    assert session.time == pytest.approx(0.5)
    window.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_restart_command_uses_validated_dat(tmp_path: Path) -> None:
    from cfdx.case_io import save_case_with_dat
    from cfdx.gui import create_application

    app = create_application(["cfdx-restart-test"])
    session = CFDXSession()
    session.case.name = "restart"
    session.case.execution.solver = sys.executable
    source = tmp_path / "solver.dat"
    source.write_text("restart-state\n", encoding="utf-8")
    case_path, dat_path = save_case_with_dat(session, tmp_path / "restart.cfdx.h5", source)

    window = CFDXMainWindow(session)
    window._case_path = case_path
    window._restart_dat = dat_path
    command = window._solver_command()

    assert command[-2:] == ["--restart", str(dat_path)]
    window._dirty = False
    window.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_loads_hdf5_dat_checkpoint(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from PySide6.QtWidgets import QFileDialog
    from cfdx.case_io import save_case
    from cfdx.dat_io import DatField, DatRestart, write_dat_hdf5
    from cfdx.gui import create_application

    app = create_application(["cfdx-hdf5-dat-test"])
    session = CFDXSession()
    session.case.name = "hdf5-dat"
    case_path = tmp_path / "case.cfdx.h5"
    save_case(session, case_path)

    dat_path = tmp_path / "solution.dat.h5"
    write_dat_hdf5(
        dat_path,
        DatRestart(
            version=2,
            cells=1,
            iteration=7,
            time=1.25,
            fields={"T": DatField("T", 1, [350.0])},
        ),
    )

    class FakeView:
        def load_cfdx_dat(self, case: str, dat: str) -> list[str]:
            assert Path(case) == case_path
            assert Path(dat) == dat_path
            return ["T"]

    window = CFDXMainWindow(session)
    window._case_path = case_path
    window.view3d = FakeView()
    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *args, **kwargs: (str(dat_path), "CFDX DAT HDF5 (*.h5)")),
    )

    assert window._open_dat_checkpoint()
    assert window._restart_dat == dat_path
    assert "iteration=7" in window.result_status.text()
    assert "time=1.25" in window.result_status.text()

    window._dirty = False
    window.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_loads_dat_as_result_without_enabling_restart(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    from PySide6.QtWidgets import QFileDialog
    from cfdx.case_io import save_case
    from cfdx.dat_io import DatField, DatRestart, write_dat_hdf5
    from cfdx.gui import create_application

    app = create_application(["cfdx-dat-result-test"])
    session = CFDXSession()
    case_path = tmp_path / "case.cfdx.h5"
    save_case(session, case_path)
    dat_path = tmp_path / "solution.dat.h5"
    write_dat_hdf5(
        dat_path,
        DatRestart(
            version=2,
            cells=1,
            iteration=12,
            time=2.5,
            fields={"p": DatField("p", 1, [101325.0])},
        ),
    )

    class FakeView:
        def load_cfdx_dat(self, case: str, dat: str) -> list[str]:
            assert Path(case) == case_path
            assert Path(dat) == dat_path
            return ["p"]

    window = CFDXMainWindow(session)
    window._case_path = case_path
    window.view3d = FakeView()
    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileName",
        staticmethod(lambda *args, **kwargs: (str(dat_path), "CFDX DAT HDF5 (*.h5)")),
    )

    assert window._open_dat_result()
    assert window._restart_dat is None
    assert window._result_source == dat_path
    assert "DAT result:" in window.result_status.text()
    assert "iteration=12" in window.result_status.text()
    assert "time=2.5" in window.result_status.text()

    window._dirty = False
    window.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_gui_refuses_to_ignore_dat_when_restart_disabled(tmp_path: Path) -> None:
    from cfdx.case_io import save_case_with_dat
    from cfdx.gui import create_application

    app = create_application(["cfdx-restart-disabled-test"])
    session = CFDXSession()
    session.case.execution.solver = sys.executable
    session.case.execution.restart_option = None
    source = tmp_path / "solver.dat"
    source.write_text("restart-state\n", encoding="utf-8")
    case_path, dat_path = save_case_with_dat(session, tmp_path / "case.cfdx.h5", source)

    window = CFDXMainWindow(session)
    window._case_path = case_path
    window._restart_dat = dat_path

    with pytest.raises(ValueError, match="no solver restart option"):
        window._solver_command()

    window._dirty = False
    window.close()
    app.quit()
