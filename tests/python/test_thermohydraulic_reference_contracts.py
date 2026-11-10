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


def test_vertical_pipe_1998_must_never_be_mislabeled_as_th07_dns() -> None:
    contracts = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
    dns = json.loads((contracts / "TH-07.json").read_text(encoding="utf-8"))
    pipe = json.loads((contracts / "TH-PIPE-UPFLOW-1998.json").read_text(encoding="utf-8"))
    assert dns["case_id"] == "TH-07"
    assert "plane channel" in dns["physics"]["model"]
    assert dns["reference"]["type"] == "DNS"
    assert pipe["case_id"] != dns["case_id"]
    assert pipe["reference"]["type"] == "EXPERIMENT"
    assert "circular pipe" in pipe["physics"]["geometry"]
    assert pipe["reference"]["source_data"] == "NOT_INGESTED"
    assert pipe["campaign_implementation"] == "UNREGISTERED_REFERENCE_PREPARATION_ONLY"


def test_reference_status_tracks_experimental_sst_production_without_certification() -> None:
    contracts = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
    for case in ("TH-01", "TH-03", "TH-07"):
        data = json.loads((contracts / f"{case}.json").read_text(encoding="utf-8"))
        status = data["status"]
        assert status["physical_validation"].startswith("BLOCKED")
        assert "SST" in status["reason"]
        assert "experimental" in status["reason"].lower()


def test_vertical_pipe_1998_must_never_be_mislabeled_as_th07_dns() -> None:
    contracts = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
    dns = json.loads((contracts / "TH-07.json").read_text(encoding="utf-8"))
    pipe = json.loads((contracts / "TH-PIPE-UPFLOW-1998.json").read_text(encoding="utf-8"))
    assert dns["case_id"] == "TH-07"
    assert "plane channel" in dns["physics"]["model"]
    assert dns["reference"]["type"] == "DNS"
    assert pipe["case_id"] != dns["case_id"]
    assert pipe["reference"]["type"] == "EXPERIMENT"
    assert "circular pipe" in pipe["physics"]["geometry"]
    assert pipe["reference"]["source_data"] == "NOT_INGESTED"
    assert pipe["campaign_implementation"] == "UNREGISTERED_REFERENCE_PREPARATION_ONLY"


def test_reference_status_tracks_experimental_sst_production_without_certification() -> None:
    contracts = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
    for case in ("TH-01", "TH-03", "TH-07"):
        data = json.loads((contracts / f"{case}.json").read_text(encoding="utf-8"))
        status = data["status"]
        assert status["physical_validation"].startswith("BLOCKED")
        assert "SST" in status["reason"]
        assert "experimental" in status["reason"].lower()


def test_vertical_pipe_1998_must_never_be_mislabeled_as_th07_dns() -> None:
    contracts = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
    dns = json.loads((contracts / "TH-07.json").read_text(encoding="utf-8"))
    pipe = json.loads(
        (contracts / "TH-PIPE-UPFLOW-1998.json").read_text(encoding="utf-8")
    )
    assert dns["case_id"] == "TH-07"
    assert "plane channel" in dns["physics"]["model"]
    assert dns["reference"]["type"] == "DNS"
    assert pipe["case_id"] != dns["case_id"]
    assert pipe["reference"]["type"] == "EXPERIMENT"
    assert "circular pipe" in pipe["physics"]["geometry"]
    assert pipe["reference"]["source_data"] == "NOT_INGESTED"
    assert pipe["campaign_implementation"] == "UNREGISTERED_REFERENCE_PREPARATION_ONLY"


def test_reference_status_tracks_experimental_sst_production_without_certification() -> (
    None
):
    contracts = ROOT / "docs/validation/14-benchmarks/09-thermohydraulics"
    for case in ("TH-01", "TH-03", "TH-07"):
        data = json.loads((contracts / f"{case}.json").read_text(encoding="utf-8"))
        status = data["status"]
        assert status["physical_validation"].startswith("BLOCKED")
        assert "SST" in status["reason"]
        assert "experimental" in status["reason"].lower()
