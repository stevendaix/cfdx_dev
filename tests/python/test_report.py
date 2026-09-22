from pathlib import Path

from cfdx.report import Report


def test_report_is_deterministic(tmp_path: Path) -> None:
    report = Report("CFDX validation", {"Run": "Iteration: 10\nResidual: 1e-6"})
    expected = "# CFDX validation\n\n## Run\n\nIteration: 10\nResidual: 1e-6\n"
    assert report.to_markdown() == expected
    path = tmp_path / "report.md"
    report.write_markdown(path)
    assert path.read_text(encoding="utf-8") == expected
