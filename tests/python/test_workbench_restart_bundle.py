import importlib.util
from pathlib import Path

import pytest


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)
def test_workbench_save_pairs_selected_restart_dat(tmp_path: Path) -> None:
    from cfdx.gui import create_application
    from cfdx.workbench import CFDXWorkbenchWindow

    from cfdx import CFDXSession

    app = create_application(["cfdx-workbench-save-bundle"])
    session = CFDXSession()
    source = tmp_path / "solver.dat"
    source.write_text(
        "CFDX-DAT 3\n"
        "cells 1\n"
        "iteration 4\n"
        "time 0.25\n"
        "identity 8 6 1 42\n"
        "field U 3\n"
        "1 2 3\n"
        "field p 1\n"
        "101325\n"
        "optional_fields\n",
        encoding="utf-8",
    )
    case_path = tmp_path / "bundle.cfdx.h5"
    window = CFDXWorkbenchWindow(session)
    window.application.project_path = case_path
    window.application.load_dat(source, use_for_restart=True)

    window._save_case()

    paired = tmp_path / "bundle.dat.h5"
    assert paired.is_file()
    assert window.application.restart_dat == paired
    assert window.application.dirty is False
    window.close()
    app.quit()
