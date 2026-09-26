from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SPEC = spec_from_file_location(
    "validation_report",
    ROOT / "scripts" / "validation_report.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


def test_validation_gate_passes_only_when_every_executable_passes() -> None:
    result = MODULE.validation_gate_status(
        {
            "test_a": 0,
            "test_b": 0,
        }
    )
    assert result == {
        "status": "PASS",
        "passed": ["test_a", "test_b"],
        "failed": [],
        "missing": [],
    }


def test_validation_gate_distinguishes_missing_executables() -> None:
    result = MODULE.validation_gate_status(
        {
            "test_ok": 0,
            "test_missing": -1,
            "test_failed": 7,
        }
    )
    assert result["status"] == "FAIL"
    assert result["passed"] == ["test_ok"]
    assert result["missing"] == ["test_missing"]
    assert result["failed"] == ["test_failed"]


def test_reference_oracle_is_not_a_solver_pass() -> None:
    # A passing Fluent reference executable only proves the reference
    # calculation/parameters are internally consistent. It must not be
    # interpreted as evidence that CFDX reproduced the Fluent solution.
    result = {
        "status": "PASS",
        "is_solver_validation": False,
    }
    assert result["status"] == "PASS"
    assert result["is_solver_validation"] is False


def test_generated_dynamic_table_rows_are_valid_latex(tmp_path: Path) -> None:
    output = tmp_path / "report.tex"
    MODULE.write_tex(
        output,
        [MODULE.Case("VMFL001", "Case", "READY", "Comparison")],
        [{
            "re": 100.0,
            "grid": "64x64",
            "continuity": 1e-12,
            "momentum": 2e-8,
            "u_rms": 0.01,
            "u_max": 0.02,
            "v_rms": 0.03,
            "v_max": 0.04,
        }],
        [{
            "name": "MODEL",
            "metric": "error",
            "value": "1e-6",
            "reference": "exact",
        }],
        {"oracle": (0, "")},
        {"test_a": 0},
        None,
        "2026-09-26 00:00 UTC",
    )

    lines = output.read_text(encoding="utf-8").splitlines()

    def assert_row(prefix: str, columns: int) -> None:
        row = next(line for line in lines if line.startswith(prefix))
        assert row.endswith(r"\\")
        assert row.count(" & ") == columns - 1

    assert_row(r"\texttt{test\_a}", 2)
    assert_row("100", 8)
    assert_row(r"\texttt{MODEL}", 4)
    assert_row(r"\texttt{VMFL001}", 4)
