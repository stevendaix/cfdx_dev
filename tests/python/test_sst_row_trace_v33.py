import importlib.util
from pathlib import Path
import pytest

S = Path(__file__).resolve().parents[2] / "scripts" / "analyze_sst_rows_v33.py"
spec = importlib.util.spec_from_file_location("sst_v33_row_audit", S)
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


def log(iterations=4, cells=(341, 329)):
    rows = []
    for i in range(1, iterations + 1):
        for cell in cells:
            for field in ("k", "omega"):
                rows.append(
                    f"SST_ROW_TRACE_V33 iteration={i} cell={cell} field={field} "
                    "full_residual=0.1 closure_error=0 source_identity_error=0 "
                    "local_relative=0.5 production_limited=0 "
                    "omega_viscosity_branch=1 branch_changed=0 "
                    "winner_accepted=0 frozen_source_derivative_valid=1 "
                    "dRk_dlogk_source=0.1 dRk_dlogomega_source=-0.2 "
                    "source_units=per_unit_volume cell_volume=0.01 "
                    "dRk_dlogk_fv_source=0.001 "
                    "dRk_dlogomega_fv_source=-0.002 qualified=0"
                )
    return "\n".join(rows)


def test_complete_diagnostic_does_not_qualify():
    result = mod.parse_log(log(), [341, 329])
    assert result["complete_inner_iterations"] == 4
    assert result["status"] == "DIAGNOSTIC_ONLY_NOT_QUALIFIED"
    assert result["outer_commits"] == "NOT_ESTABLISHED_FROM_ROW_TRACE"


def test_missing_iteration_is_rejected():
    with pytest.raises(ValueError, match="Missing inner"):
        mod.parse_log(log().replace("iteration=2 ", "iteration=5 "), [341, 329])


def test_truncated_omega_row_is_rejected():
    x = log()
    x = "\n".join(row for row in x.splitlines() if not
                  ("iteration=3 " in row and "cell=329 field=omega" in row))
    with pytest.raises(ValueError, match="Missing SST"):
        mod.parse_log(x, [341, 329])


def test_duplicate_row_is_rejected():
    x = log()
    with pytest.raises(ValueError, match="Duplicate row"):
        mod.parse_log(x + "\n" + x.splitlines()[0], [341, 329])


def test_nonfinite_value_is_rejected():
    with pytest.raises(ValueError, match="Non-finite"):
        mod.parse_log(log().replace("full_residual=0.1", "full_residual=nan", 1), [341, 329])


def test_false_qualification_is_rejected():
    with pytest.raises(ValueError, match="cannot certify"):
        mod.parse_log(log().replace("qualified=0", "qualified=1", 1), [341, 329])


def test_empty_log_is_rejected():
    with pytest.raises(ValueError, match="No V33"):
        mod.parse_log("BFS converged", [341, 329])


def test_duplicate_requested_cells_are_rejected():
    with pytest.raises(ValueError, match="distinct"):
        mod.parse_log(log(), [341, 341])


def test_first_attempt_must_be_explicit_when_a_retry_exists():
    combined = log(3) + "\n" + log(2)
    with pytest.raises(ValueError, match="Duplicate row"):
        mod.parse_log(combined, [341, 329])
    result = mod.parse_log(combined, [341, 329], first_attempt=True)
    assert result["complete_inner_iterations"] == 3


def test_fv_source_scaling_is_fail_closed():
    invalid = log().replace("dRk_dlogomega_fv_source=-0.002", "dRk_dlogomega_fv_source=-0.003", 1)
    with pytest.raises(ValueError, match="volume scaling"):
        mod.parse_log(invalid, [341, 329])


def test_fv_source_units_are_explicit():
    with pytest.raises(ValueError, match="Source derivative units"):
        mod.parse_log(log().replace("source_units=per_unit_volume", "source_units=unknown", 1), [341, 329])
