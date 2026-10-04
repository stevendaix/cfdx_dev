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
23: N8_SCHUR case=0 cond_inf_Auu=1.51309 exact_oracle_discrepancy=0.658375 exact_solve_backward_error=5.64312e-17 machine_epsilon=2.22045e-16
23: N8_SCHUR case=0 method=LSC algebra_error=3.08719e-16 exact_schur_error=8.5631
23: n8_schur_benchmark cells=64 unknowns=256 nnz=766 schur_nnz=190 true_residual=1.47002e-16 iterations=1 setup_us=44 solve_us=259 pressure_coarse_size=64 hierarchy_builds=1 numeric_updates=0 factorization=ilut
23: n8_schur_benchmark_lifecycle cells=64 coefficient_changed=true hierarchy_builds_before=1 hierarchy_builds_after=2 numeric_updates=1 updated_iterations=1 updated_true_residual=1.41362e-16 graph_change_rebuild=explicit_setup_required
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
    # The quantitative test emits one oracle record plus one record per method.
    assert [record["case"] for record in evidence["schur_quantitative"]] == [0, 0]
    assert evidence["schur_quantitative"][0]["cond_inf_Auu"] == 1.51309
    assert evidence["schur_quantitative"][0]["exact_solve_backward_error"] == 5.64312e-17
    assert evidence["schur_quantitative"][1]["method"] == "LSC"
    assert evidence["schur_quantitative"][1]["exact_schur_error"] == 8.5631
    # The lifecycle prefix extends the benchmark prefix, so it must not be counted
    # as a benchmark record as well.
    assert [record["cells"] for record in evidence["schur_production"]] == [64, 64]
    assert evidence["schur_production"][0]["factorization"] == "ilut"
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


ORACLE_RECORD = {
    "case": 0,
    "cond_inf_Auu": 1.51309,
    "exact_oracle_discrepancy": 0.658375,
    "exact_solve_backward_error": 5.64312e-17,
    "machine_epsilon": 2.22045e-16,
}
METHOD_RECORD = {
    "case": 0,
    "method": "LSC",
    "algebra_error": 3.08719e-16,
    "exact_schur_error": 8.5631,
}
BENCHMARK_RECORD = {
    "cells": 64,
    "unknowns": 256,
    "nnz": 766,
    "schur_nnz": 190,
    "true_residual": 1.47002e-16,
    "iterations": 1,
    "setup_us": 44,
    "solve_us": 259,
    "pressure_coarse_size": 64,
    "hierarchy_builds": 1,
    "numeric_updates": 0,
}
LIFECYCLE_RECORD = {
    "cells": 64,
    "coefficient_changed": True,
    "hierarchy_builds_before": 1,
    "hierarchy_builds_after": 2,
    "numeric_updates": 1,
    "updated_iterations": 1,
    "updated_true_residual": 1.41362e-16,
    "graph_change_rebuild": "explicit_setup_required",
}


def complete_evidence() -> dict[str, object]:
    return {
        "physical_model_results": [
            {
                "model": "Couette",
                "solver_converged": True,
                "iterations": 10,
                "gates_failed": 0,
            }
        ],
        "schur_quantitative": [ORACLE_RECORD, METHOD_RECORD],
        "schur_production": [BENCHMARK_RECORD, LIFECYCLE_RECORD],
    }


def test_evidence_coverage_accepts_complete_structured_records() -> None:
    coverage = MODULE.audit_evidence_coverage(complete_evidence())
    assert coverage["status"] == "COMPLETE"
    assert coverage["records_checked"] == 5
    assert coverage["records_complete"] == 5
    assert coverage["missing_fields"] == {}
    assert coverage["categories_without_records"] == []
    assert coverage["malformed_records"] == {}
    assert coverage["policy"] == "diagnostic_only"


def test_evidence_coverage_separates_the_two_schur_quantitative_shapes() -> None:
    # The quantitative test prints an oracle record and one record per method.
    # Requiring the oracle fields of both shapes would report every healthy
    # campaign as incomplete, so each shape is audited against its own layout.
    coverage = MODULE.audit_evidence_coverage(complete_evidence())
    assert "schur_quantitative_oracle" not in coverage["missing_fields"]
    assert "schur_quantitative_method" not in coverage["missing_fields"]


def test_evidence_coverage_reports_missing_fields_without_reinterpreting_results() -> None:
    evidence = complete_evidence()
    evidence["physical_model_results"] = [
        {"model": "Couette", "solver_converged": True, "iterations": 10, "gates_failed": 0},
        {"model": "Ghia", "solver_converged": False, "execution_exception": "boom"},
    ]
    coverage = MODULE.audit_evidence_coverage(evidence)
    assert coverage["status"] == "INCOMPLETE"
    assert coverage["records_checked"] == 6
    assert coverage["records_complete"] == 5
    assert coverage["missing_fields"]["physical_model_results"] == [
        {"index": 1, "fields": ["iterations", "gates_failed"]}
    ]
    assert coverage["policy"] == "diagnostic_only"


def test_evidence_coverage_flags_categories_without_records() -> None:
    coverage = MODULE.audit_evidence_coverage({})
    assert coverage["status"] == "INCOMPLETE"
    assert coverage["records_checked"] == 0
    assert coverage["categories_without_records"] == [
        "physical_model_results",
        "schur_quantitative_oracle",
        "schur_quantitative_method",
        "schur_production",
        "schur_lifecycle",
    ]


def test_evidence_coverage_reports_a_truncated_campaign() -> None:
    # main() stops at the first failing gate, so a partial campaign legitimately
    # has no production evidence; that absence has to stay visible.
    evidence = {"schur_quantitative": [ORACLE_RECORD, METHOD_RECORD]}
    coverage = MODULE.audit_evidence_coverage(evidence)
    assert coverage["status"] == "INCOMPLETE"
    assert coverage["records_checked"] == 2
    assert coverage["categories_without_records"] == [
        "physical_model_results",
        "schur_production",
        "schur_lifecycle",
    ]


def test_evidence_coverage_accounts_for_malformed_records() -> None:
    coverage = MODULE.audit_evidence_coverage({"schur_production": ["not a record"]})
    assert coverage["status"] == "INCOMPLETE"
    assert coverage["malformed_records"] == {"schur_production": [0]}
    assert coverage["records_checked"] == 0
    assert coverage["records_complete"] == 0


def test_evidence_coverage_tolerates_a_non_list_evidence_section() -> None:
    coverage = MODULE.audit_evidence_coverage({"physical_model_results": "oops"})
    assert coverage["status"] == "INCOMPLETE"
    assert coverage["categories_without_records"] == [
        "physical_model_results",
        "schur_quantitative_oracle",
        "schur_quantitative_method",
        "schur_production",
        "schur_lifecycle",
    ]


def test_evidence_coverage_is_complete_on_a_captured_campaign() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    assert MODULE.audit_evidence_coverage(evidence)["status"] == "COMPLETE"

