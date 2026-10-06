from __future__ import annotations

import csv
import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CMP = ROOT / "scripts" / "compare_validation_results.py"
AGG = ROOT / "scripts" / "aggregate_validation_runs.py"


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["case", "Cd", "Cl"])
        writer.writeheader()
        writer.writerows(rows)


def run(script: Path, *args: str) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *args],
        capture_output=True,
        text=True,
        check=False,
    )


def test_comparator_success_zero_reference_and_json(tmp_path: Path) -> None:
    computed = tmp_path / "computed.csv"
    reference = tmp_path / "reference.csv"
    output = tmp_path / "report.json"
    write_csv(computed, [{"case": "naca", "Cd": "0.12", "Cl": "1e-6"}])
    write_csv(reference, [{"case": "naca", "Cd": "0.12", "Cl": "0"}])

    result = run(
        CMP, "--computed", str(computed), "--reference", str(reference),
        "--qoi", "Cd", "Cl", "--output", str(output),
    )
    assert result.returncode == 0, result.stderr
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["results"][0]["absolute_error"] == 0.0
    assert report["results"][1]["relative_error"] is None


def test_comparator_rejects_missing_case_and_qoi(tmp_path: Path) -> None:
    computed = tmp_path / "computed.csv"
    reference = tmp_path / "reference.csv"
    write_csv(computed, [{"case": "only-computed", "Cd": "0.1", "Cl": "0"}])
    write_csv(reference, [{"case": "only-reference", "Cd": "0.1", "Cl": "0"}])
    result = run(
        CMP, "--computed", str(computed), "--reference", str(reference),
        "--qoi", "Cd",
    )
    assert result.returncode != 0
    assert "case-set mismatch" in result.stderr


def test_comparator_rejects_nonfinite_and_duplicate_cases(tmp_path: Path) -> None:
    computed = tmp_path / "computed.csv"
    reference = tmp_path / "reference.csv"
    write_csv(computed, [
        {"case": "naca", "Cd": "nan", "Cl": "0"},
    ])
    write_csv(reference, [{"case": "naca", "Cd": "0.1", "Cl": "0"}])
    result = run(
        CMP, "--computed", str(computed), "--reference", str(reference),
        "--qoi", "Cd",
    )
    assert result.returncode != 0
    assert "not finite" in result.stderr

    write_csv(computed, [
        {"case": "naca", "Cd": "0.1", "Cl": "0"},
        {"case": "naca", "Cd": "0.2", "Cl": "0"},
    ])
    result = run(
        CMP, "--computed", str(computed), "--reference", str(reference),
        "--qoi", "Cd",
    )
    assert result.returncode != 0
    assert "duplicate case" in result.stderr


def test_aggregator_is_deterministic_and_preserves_provenance(tmp_path: Path) -> None:
    input_dir = tmp_path / "runs"
    input_dir.mkdir()
    write_csv(input_dir / "b.csv", [{"case": "b", "Cd": "0.2", "Cl": "0"}])
    write_csv(input_dir / "a.csv", [{"case": "a", "Cd": "0.1", "Cl": "0"}])
    output = tmp_path / "out" / "all.csv"

    result = run(
        AGG, "--input-dir", str(input_dir), "--output", str(output),
    )
    assert result.returncode == 0, result.stderr
    rows = list(csv.DictReader(output.open(newline="", encoding="utf-8")))
    assert [row["case"] for row in rows] == ["a", "b"]
    assert rows[0]["source_file"].endswith("a.csv")

    result = run(
        AGG, "--input-dir", str(input_dir), "--output", str(output),
    )
    assert result.returncode == 0, result.stderr


def test_aggregator_rejects_duplicate_keys(tmp_path: Path) -> None:
    input_dir = tmp_path / "runs"
    input_dir.mkdir()
    write_csv(input_dir / "a.csv", [{"case": "same", "Cd": "0.1", "Cl": "0"}])
    write_csv(input_dir / "b.csv", [{"case": "same", "Cd": "0.2", "Cl": "0"}])
    result = run(
        AGG, "--input-dir", str(input_dir), "--output", str(tmp_path / "out.csv"),
    )
    assert result.returncode != 0
    assert "duplicate validation key" in result.stderr
