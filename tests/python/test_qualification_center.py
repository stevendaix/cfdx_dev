"""Strict registry parsing, provenance and parallel qualification browsing."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

from cfdx.qualification_center import (
    MAX_REGISTRY_BYTES,
    STATUSES,
    discover_qualification_registry,
    read_qualification_registry,
)


def _fixture_registry() -> dict:
    return {
        "schema_version": "1.0",
        "campaign": "Fixture engineering campaign",
        "status_policy": {status: f"Meaning of {status}" for status in STATUSES},
        "cases": [
            {
                "id": "EXT-NACA0012", "domain": "External flow", "case": "NACA0012",
                "status": "DIAGNOSTIC", "qoi": "Cd, Cl, Cp", "runtime": "unverified",
            },
            {
                "id": "TUR-CHANNEL", "domain": "Turbulence", "case": "Channel",
                "status": "BLOCKED", "qoi": "Cf", "runtime": "not run",
            },
        ],
    }


def _write(tmp_path: Path, data: dict) -> Path:
    path = tmp_path / "CFDX_QUALIFICATION_REGISTRY.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    return path


def test_actual_repository_registry_is_browsable_without_claiming_head_evidence() -> None:
    repo = Path(__file__).resolve().parents[2]
    source = discover_qualification_registry(repo)
    assert source == repo / "docs/validation/CFDX_QUALIFICATION_REGISTRY.json"
    assert source is not None
    registry = read_qualification_registry(source)
    assert len(registry.cases) == 25
    assert registry.counts() == {
        "PLANNED": 2, "READY": 10, "RUNNING": 0,
        "DIAGNOSTIC": 6, "PASS": 3, "BLOCKED": 4,
    }
    assert {case.id for case in registry.filter("NACA")} == {
        "EXT-NACA0012", "TUR-NACA0012", "N9-S8-NACA0012",
    }
    assert len(registry.filter(status="DIAGNOSTIC")) == 6
    assert all(case.status != "PASS" for case in registry.filter("NACA"))
    assert registry.source_signature[2] > 0


def test_filter_is_case_insensitive_and_status_domain_exact(tmp_path: Path) -> None:
    registry = read_qualification_registry(_write(tmp_path, _fixture_registry()))
    assert [case.id for case in registry.filter("naca")] == ["EXT-NACA0012"]
    assert [case.id for case in registry.filter(domain="Turbulence")] == [
        "TUR-CHANNEL"
    ]
    assert registry.filter("naca", status="BLOCKED") == ()
    assert registry.policy_for("DIAGNOSTIC") == "Meaning of DIAGNOSTIC"
    assert registry.cases[0].details[0] == ("qoi", "Cd, Cl, Cp")


@pytest.mark.parametrize(
    "edit,expected",
    [
        (lambda x: x.update(schema_version="2.0"), "schema version"),
        (lambda x: x["cases"][0].update(status="PASSED"), "unsupported qualification status"),
        (lambda x: x["cases"][1].update(id="EXT-NACA0012"), "duplicate qualification case id"),
        (lambda x: x["cases"][0].update(id="../../tmp/foo"), "invalid qualification case id"),
        (lambda x: x["cases"][0].update(qoi={"error": "vector"}), "bounded string"),
        (lambda x: x["status_policy"].pop("DIAGNOSTIC"), "incomplete status policy"),
        (lambda x: x.update(cases="not a list"), "invalid qualification registry case list"),
        (lambda x: x["cases"].append("not an object"), "must be an object"),
    ],
)
def test_reject_malformed_or_ambiguous_registry(tmp_path: Path, edit, expected):
    data = _fixture_registry()
    edit(data)
    with pytest.raises(ValueError, match=expected):
        read_qualification_registry(_write(tmp_path, data))


def test_reject_duplicate_json_keys_and_nan(tmp_path: Path) -> None:
    source = tmp_path / "bad.json"
    source.write_text('{"schema_version":"1.0","schema_version":"1.0"}')
    with pytest.raises(ValueError, match="duplicate JSON key"):
        read_qualification_registry(source)
    source.write_text('{"schema_version":NaN}')
    with pytest.raises(ValueError, match="non-finite JSON constant"):
        read_qualification_registry(source)


def test_registry_file_budget_is_fail_closed(tmp_path: Path) -> None:
    path = tmp_path / "oversize.json"
    with path.open("wb") as handle:
        handle.truncate(MAX_REGISTRY_BYTES + 1)
    with pytest.raises(ValueError, match="4 MiB"):
        read_qualification_registry(path)


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="Qt optional")
def test_qualification_panel_parallel_case_filters_and_failed_reload(tmp_path: Path) -> None:
    from cfdx.gui import create_application
    from cfdx.qualification_panel import QualificationCenterPanel

    app = create_application(["cfdx-qualification-test"])
    panel = QualificationCenterPanel()
    path = _write(tmp_path, _fixture_registry())
    panel.load_registry(path)
    assert panel.table.rowCount() == 2
    assert "not certify" in panel.caution.text()
    panel.search.setText("naca")
    assert panel.table.rowCount() == 1
    assert "EXT-NACA0012" in panel.details.toPlainText()
    assert "NOT CURRENT-HEAD CERTIFIED" in panel.details.toPlainText()
    panel.search.clear()
    panel.status_filter.setCurrentText("BLOCKED")
    assert panel.table.rowCount() == 1
    assert "TUR-CHANNEL" in panel.details.toPlainText()
    old_registry = panel.registry
    path.write_text("{not-json", encoding="utf-8")
    panel.reload()
    assert panel.registry is old_registry
    assert "previous snapshot retained" in panel.message.text()
    panel.close()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="Qt optional")
def test_workbench_validate_workspace_and_registry_dock() -> None:
    from PySide6.QtWidgets import QDockWidget
    from cfdx.gui import create_application
    from cfdx.session import CFDXSession
    from cfdx.workbench import CFDXWorkbenchWindow

    app = create_application(["cfdx-qualification-workbench-test"])
    window = CFDXWorkbenchWindow(CFDXSession())
    assert window.findChild(QDockWidget, "workbench.dock.qualification") is not None
    window.workspace_actions["validate"].trigger()
    assert window._workspace == "validate"
    assert not window.findChild(QDockWidget, "workbench.dock.qualification").isHidden()
    window.workspace_actions["setup"].trigger()
    window.close()
    app.quit()


def test_headless_cli_does_not_promote_registry_status_to_solver_evidence(
    tmp_path: Path,
) -> None:
    import subprocess
    import sys

    repo = Path(__file__).resolve().parents[2]
    source = _write(tmp_path, _fixture_registry())
    run = subprocess.run(
        [sys.executable, str(repo / "scripts/cfdx_qualification_status.py"),
         "--registry", str(source), "--format", "json", "--query", "naca"],
        capture_output=True, text=True, check=True,
    )
    payload = json.loads(run.stdout)
    assert payload["status_semantics"] == (
        "declarations_only_not_certified_on_current_head"
    )
    assert payload["total_cases"] == 2
    assert payload["filtered_cases"] == [
        {
            "id": "EXT-NACA0012", "case": "NACA0012",
            "domain": "External flow", "declared_status": "DIAGNOSTIC",
            "head_evidence": "unverified",
        }
    ]
