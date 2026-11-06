"""Regression tests for frozen turbulent thermohydraulic reference contracts."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
RUNNER = ROOT / "scripts" / "validate_thermohydraulic_reference.py"


def test_all_thermohydraulic_reference_contracts_are_consistent(tmp_path: Path) -> None:
    report_path = tmp_path / "thermohydraulic-reference-contracts.json"
    result = subprocess.run(
        [
            sys.executable,
            str(RUNNER),
            "--output",
            str(report_path),
        ],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    report = json.loads(report_path.read_text(encoding="utf-8"))
    assert report["not_a_solver_validation"] is True
    assert [case["case_id"] for case in report["cases"]] == [
        "TH-01",
        "TH-03",
        "TH-07",
        "TH-08",
    ]
    assert all(case["contract_status"] == "PASS" for case in report["cases"])
    assert all(
        case["physical_validation_status"] == "BLOCKED_NOT_A_SOLVER_RUN"
        for case in report["cases"]
    )
