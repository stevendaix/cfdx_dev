import importlib.util
from pathlib import Path

import pytest


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 optional")
def test_workbench_save_pairs_selected_restart_dat(tmp_path: Path) -> None:
    from cfdx import CFDXSession
    from cfdx.gui import create_application
    from cfdx.workbench import CFDXWorkbenchWindow

    app = create_application(["cfdx-workbench-save-bundle"])
    session = CFDXSession()
    source = tmp_path / "solver.dat"
    source.write_text(
        "CFDX-DAT 1 cells 1 iteration 4 time 0.25 field p 1 101325\n",
        encoding="utf-8",
    )
    case_path = tmp_path / "bundle.cfdx.h5"
    window = CFDXWorkbenchWindow(session)
    window.application.project_path = case_path
    window._restart_dat = source

    window._save_case()

    paired = tmp_path / "bundle.dat.h5"
    assert paired.is_file()
    assert window._restart_dat == paired
    assert window.application.dirty is False
    window.close()
    app.quit()
