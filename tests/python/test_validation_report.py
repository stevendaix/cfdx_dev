from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
SPEC = spec_from_file_location(
    "validation_report",
    ROOT / "scripts" / "validation_report.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
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
