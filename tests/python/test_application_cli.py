from __future__ import annotations

import json
from pathlib import Path

from cfdx.cli.application_cli import main


def _case_file(tmp_path: Path) -> Path:
    from cfdx.case_io import save_case
    from cfdx.session import CFDXSession

    case = tmp_path / "cli.cfdx.h5"
    save_case(CFDXSession(), case)
    return case


def test_cli_status_uses_headless_application(tmp_path: Path, capsys) -> None:
    case = _case_file(tmp_path)
    assert main(["status", str(case)]) == 0
    output = capsys.readouterr().out
    assert "CFDX | Case: untitled" in output
    assert "State: CREATED" in output


def test_cli_info_exposes_application_state(tmp_path: Path, capsys) -> None:
    case = _case_file(tmp_path)
    assert main(["info", str(case)]) == 0
    data = json.loads(capsys.readouterr().out)
    assert data["simulation_state"] == "CREATED"
    assert data["project_path"].endswith("cli.cfdx.h5")


def test_cli_results_uses_application_results_contract(tmp_path: Path, capsys) -> None:
    case = _case_file(tmp_path)
    results = tmp_path / "results"
    results.mkdir()
    (results / "result_0001.vtu").write_text(
        """<?xml version="1.0"?><VTKFile type="UnstructuredGrid"><UnstructuredGrid><Piece NumberOfPoints="0" NumberOfCells="0"/></UnstructuredGrid></VTKFile>""",
        encoding="utf-8",
    )
    assert main(["results", str(case), str(results)]) == 0
    data = json.loads(capsys.readouterr().out)
    assert len(data["frames"]) == 1
