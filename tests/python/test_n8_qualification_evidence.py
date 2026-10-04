from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
SPEC = spec_from_file_location(
    "n8_physical_qualification",
    ROOT / "scripts" / "n8_physical_qualification.py",
)
assert SPEC is not None and SPEC.loader is not None
MODULE = module_from_spec(SPEC)
sys.modules[SPEC.name] = MODULE
SPEC.loader.exec_module(MODULE)


CTEST_OUTPUT = """    Start 21: test_n8_pressure_velocity_matrix
21: Test command: /build/test_phase9_acceptance
21: MODEL_RESULT COUPLED/BlockSchur/upwind/bounded solver_converged=true iterations=17 profile_L2=1.3586e-09 Umax=0.96875 |Uy|max=5.95589e-10 |Uz|max=0 gates_failed=0
21: MODEL_FAILURES PISO/upwind/bounded continuity_linf Umax
22: MODEL_RESULT Couette solver_converged=false iterations=3 gates_failed=2
22: MODEL_FAILURES Couette solver_not_converged execution_exception
23: N8_SCHUR case=0 cond_inf_Auu=1.51309 exact_solve_backward_error=5.64312e-17
23: n8_schur_benchmark cells=64 nnz=766 factorization=full
23: n8_schur_benchmark_lifecycle cells=64 coefficient_changed=true
1/3 Test #21: test_n8_pressure_velocity_matrix ....***Passed    1.23 sec
"""


def test_model_results_keep_distinct_pipe_prefixed_keys() -> None:
    results = MODULE.parse_key_value_records(CTEST_OUTPUT, "MODEL_RESULT")
    coupled = results[0]
    assert coupled["model"] == "COUPLED/BlockSchur/upwind/bounded"
    assert coupled["solver_converged"] is True
    assert coupled["iterations"] == 17
    assert coupled["profile_L2"] == 1.3586e-09
    assert coupled["|Uy|max"] == 5.95589e-10
    assert coupled["|Uz|max"] == 0
    assert coupled["gates_failed"] == 0


def test_model_failures_retain_reported_gate_names() -> None:
    failures = MODULE.parse_key_value_records(CTEST_OUTPUT, "MODEL_FAILURES")
    assert failures == [
        {"model": "PISO/upwind/bounded", "gates": ["continuity_linf", "Umax"]},
        {"model": "Couette", "gates": ["solver_not_converged", "execution_exception"]},
    ]


def test_records_without_a_payload_are_ignored() -> None:
    assert MODULE.parse_key_value_records("MODEL_RESULT\nMODEL_FAILURES\n", "MODEL_RESULT") == []


def test_evidence_routes_each_emitter_to_its_own_list() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    assert [record["model"] for record in evidence["physical_model_results"]] == [
        "COUPLED/BlockSchur/upwind/bounded",
        "Couette",
    ]
    assert len(evidence["physical_model_failures"]) == 2
    assert [record["case"] for record in evidence["schur_quantitative"]] == [0]
    assert evidence["schur_quantitative"][0]["cond_inf_Auu"] == 1.51309
    assert evidence["schur_quantitative"][0]["exact_solve_backward_error"] == 5.64312e-17
    # The lifecycle prefix extends the benchmark prefix, so it must not be counted
    # as a benchmark record as well.
    assert [record["cells"] for record in evidence["schur_production"]] == [64, 64]
    assert evidence["schur_production"][0]["factorization"] == "full"
    assert evidence["schur_production"][1]["coefficient_changed"] is True


def test_plain_output_without_the_ctest_prefix_is_parsed() -> None:
    line = "MODEL_RESULT Couette solver_converged=true |Uy|max=1.5 |Uz|max=0.5"
    assert MODULE.parse_key_value_records(line, "MODEL_RESULT") == [
        {
            "model": "Couette",
            "solver_converged": True,
            "|Uy|max": 1.5,
            "|Uz|max": 0.5,
        }
    ]


def test_free_form_values_absorb_the_rest_of_the_line() -> None:
    line = "MODEL_RESULT Ghia execution_exception=solver did not converge"
    assert MODULE.parse_key_value_records(line, "MODEL_RESULT") == [
        {"model": "Ghia", "execution_exception": "solver did not converge"}
    ]


def test_non_numeric_values_stay_strings_when_records_follow() -> None:
    line = "n8_schur_benchmark factorization=ilut cells=10 krylov=FGMRES"
    assert MODULE.parse_key_value_records(line, "n8_schur_benchmark") == [
        {"factorization": "ilut", "cells": 10, "krylov": "FGMRES"}
    ]


def test_evidence_coverage_accepts_complete_structured_records() -> None:
    evidence = {
        "physical_model_results": [{"model": "Couette", "solver_converged": True, "iterations": 10, "gates_failed": 0}],
        "schur_quantitative": [{"case": 0, "cond_inf_Auu": 1.5, "exact_solve_backward_error": 1e-16, "machine_epsilon": 2.22e-16}],
        "schur_production": [{"cells": 64, "unknowns": 256, "nnz": 1000, "schur_nnz": 200, "true_residual": 1e-12, "iterations": 4, "setup_us": 10, "solve_us": 20, "pressure_coarse_size": 8, "hierarchy_builds": 1, "numeric_updates": 0}],
    }
    coverage = MODULE.audit_evidence_coverage(evidence)
    assert coverage["status"] == "COMPLETE"
    assert coverage["records_checked"] == 3
    assert coverage["records_complete"] == 3
    assert coverage["missing_fields"] == {}


def test_evidence_coverage_reports_missing_fields_without_reinterpreting_results() -> None:
    evidence = {"physical_model_results": [{"model": "Couette", "solver_converged": True}], "schur_quantitative": [], "schur_production": []}
    coverage = MODULE.audit_evidence_coverage(evidence)
    assert coverage["status"] == "INCOMPLETE"
    assert coverage["records_checked"] == 1
    assert coverage["records_complete"] == 0
    assert coverage["missing_fields"]["physical_model_results"] == [{"index": 0, "fields": ["iterations", "gates_failed"]}]
    assert coverage["policy"] == "diagnostic_only"
