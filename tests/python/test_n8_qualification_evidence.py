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
21: MODEL_CONFIG algorithm=COUPLED/BlockSchur/upwind/bounded nx=8 ny=16 bounded=true preconditioner=coupled_block_schur preconditioner_id=11 alpha_u=0.7 alpha_p=0.3 pressure_correctors=1 fractional_steps=1
21: MODEL_CONFIG algorithm=COUPLED/MGR/upwind/bounded nx=8 ny=16 bounded=true preconditioner=mgr preconditioner_id=13 alpha_u=0.7 alpha_p=0.3 pressure_correctors=1 fractional_steps=1
21: MODEL_PLAN algorithm=COUPLED/BlockSchur/upwind/bounded pressure_requested_krylov=auto pressure_requested_preconditioner=auto coupled_resolved=true coupled_krylov=fgmres coupled_preconditioner=coupled_block_schur pressure_resolved=false pressure_krylov=none pressure_preconditioner=none pressure_null_space=none
21: MODEL_PLAN algorithm=COUPLED/MGR/upwind/bounded pressure_requested_krylov=auto pressure_requested_preconditioner=native_amg coupled_resolved=true coupled_krylov=fgmres coupled_preconditioner=mgr pressure_resolved=false pressure_krylov=none pressure_preconditioner=none pressure_null_space=none
21: MODEL_RESULT COUPLED/BlockSchur/upwind/bounded solver_converged=true iterations=17 profile_L2=1.3586e-09 Umax=0.96875 |Uy|max=5.95589e-10 |Uz|max=0 gates_failed=0
21: MODEL_FAILURES COUPLED/MGR/upwind/bounded continuity_linf Umax
21: MODEL_RESULT Couette solver_converged=false iterations=3 gates_failed=2
21: MODEL_FAILURES Couette solver_not_converged execution_exception
21: MODEL_SUMMARY successful=1 failed=2
22: N8_SCHUR case=0 cond_inf_Auu=1.51309 exact_oracle_discrepancy=0.658375 exact_solve_backward_error=5.64312e-17 machine_epsilon=2.22045e-16
22: N8_SCHUR case=0 method=LSC algebra_error=3.08719e-16 exact_schur_error=8.5631
22: n8_schur_benchmark cells=64 unknowns=256 nnz=766 schur_nnz=190 true_residual=1.47002e-16 iterations=1 setup_us=44 solve_us=259 pressure_coarse_size=64 hierarchy_builds=1 numeric_updates=0 factorization=ilut
22: n8_schur_benchmark_lifecycle cells=64 coefficient_changed=true hierarchy_builds_before=1 hierarchy_builds_after=2 numeric_updates=1 updated_iterations=1 updated_true_residual=1.41362e-16 graph_change_rebuild=explicit_setup_required
1/3 Test #21: test_n8_pressure_velocity_matrix ....***Passed    1.23 sec
"""

MODEL_CONFIGURATION = {
    "algorithm": "COUPLED/BlockSchur/upwind/bounded",
    "nx": 8,
    "ny": 16,
    "bounded": True,
    "preconditioner": "coupled_block_schur",
    "preconditioner_id": 11,
    "alpha_u": 0.7,
    "alpha_p": 0.3,
    "pressure_correctors": 1,
    "fractional_steps": 1,
}


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
        {"model": "COUPLED/MGR/upwind/bounded", "gates": ["continuity_linf", "Umax"]},
        {"model": "Couette", "gates": ["solver_not_converged", "execution_exception"]},
    ]


def test_records_without_a_payload_are_ignored() -> None:
    assert MODULE.parse_key_value_records("MODEL_RESULT\nMODEL_FAILURES\n", "MODEL_RESULT") == []


def test_run_configuration_names_the_requested_preconditioner() -> None:
    records = MODULE.parse_key_value_records(CTEST_OUTPUT, "MODEL_CONFIG")
    assert records[0] == MODEL_CONFIGURATION
    assert records[1]["algorithm"] == "COUPLED/MGR/upwind/bounded"
    assert records[1]["preconditioner"] == "mgr"
    assert MODULE.parse_key_value_records(CTEST_OUTPUT, "MODEL_SUMMARY") == [
        {"successful": 1, "failed": 2}
    ]


def test_evidence_routes_each_emitter_to_its_own_list() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    assert evidence["run_configuration"] == [
        MODEL_CONFIGURATION,
        {**MODEL_CONFIGURATION, "algorithm": "COUPLED/MGR/upwind/bounded",
         "preconditioner": "mgr", "preconditioner_id": 13},
    ]
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


def test_model_resolution_accepts_a_fully_accounted_campaign() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    resolution = MODULE.audit_model_resolution(evidence)
    # Couette is resolved but never configured, so the audit stays incomplete.
    assert resolution["status"] == "INCOMPLETE"
    assert resolution["undeclared_models"] == ["Couette"]
    assert resolution["unresolved_models"] == []
    assert resolution["tallies_agree"] is True
    assert resolution["policy"] == "diagnostic_only"


def test_model_resolution_flags_a_configured_model_without_a_verdict() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    del evidence["physical_model_results"][:]
    del evidence["physical_model_failures"][:]
    resolution = MODULE.audit_model_resolution(evidence)
    assert resolution["status"] == "INCOMPLETE"
    assert resolution["unresolved_models"] == [
        "COUPLED/BlockSchur/upwind/bounded",
        "COUPLED/MGR/upwind/bounded",
    ]
    assert resolution["tallies_agree"] is False


def test_model_resolution_detects_a_truncated_capture() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    # One failure record disappears from the capture: the model still has its
    # result record, so the observed failure count drops below the printed tally.
    del evidence["physical_model_failures"][1]
    resolution = MODULE.audit_model_resolution(evidence)
    assert resolution["status"] == "INCOMPLETE"
    assert resolution["observed_tally"] == {"resolved": 3, "failed": 1, "total": 3}
    assert resolution["reported_tally"] == {"successful": 1, "failed": 2, "total": 3}
    assert resolution["tallies_agree"] is False


def test_model_resolution_is_complete_without_a_reported_tally() -> None:
    evidence = {
        "run_configuration": [MODEL_CONFIGURATION],
        "physical_model_results": [{"model": MODEL_CONFIGURATION["algorithm"]}],
        "physical_model_failures": [],
    }
    resolution = MODULE.audit_model_resolution(evidence)
    assert resolution["status"] == "COMPLETE"
    assert resolution["models_configured"] == 1
    assert resolution["models_resolved"] == 1
    assert resolution["reported_tally"] is None
    assert resolution["tallies_agree"] is None


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


RESOLVED_PLAN = {
    "algorithm": "COUPLED/BlockSchur/upwind/bounded",
    "pressure_requested_krylov": "auto",
    "pressure_requested_preconditioner": "auto",
    "coupled_resolved": True,
    "coupled_krylov": "fgmres",
    "coupled_preconditioner": "coupled_block_schur",
    "pressure_resolved": False,
    "pressure_krylov": "none",
    "pressure_preconditioner": "none",
    "pressure_null_space": "none",
}


def complete_evidence() -> dict[str, object]:
    return {
        "run_configuration": [MODEL_CONFIGURATION],
        "resolved_plans": [RESOLVED_PLAN],
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
    assert coverage["records_checked"] == 7
    assert coverage["records_complete"] == 7
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
    assert coverage["records_checked"] == 8
    assert coverage["records_complete"] == 7
    assert coverage["missing_fields"]["physical_model_results"] == [
        {"index": 1, "fields": ["iterations", "gates_failed"]}
    ]
    assert coverage["policy"] == "diagnostic_only"


def test_evidence_coverage_flags_categories_without_records() -> None:
    coverage = MODULE.audit_evidence_coverage({})
    assert coverage["status"] == "INCOMPLETE"
    assert coverage["records_checked"] == 0
    assert coverage["categories_without_records"] == [
        "run_configuration",
        "resolved_plans",
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
        "run_configuration",
        "resolved_plans",
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
        "run_configuration",
        "resolved_plans",
        "physical_model_results",
        "schur_quantitative_oracle",
        "schur_quantitative_method",
        "schur_production",
        "schur_lifecycle",
    ]


def test_evidence_coverage_is_complete_on_a_captured_campaign() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    assert MODULE.audit_evidence_coverage(evidence)["status"] == "COMPLETE"



def test_resolved_plans_are_captured_with_both_requests() -> None:
    records = MODULE.parse_key_value_records(CTEST_OUTPUT, "MODEL_PLAN")
    assert records[0] == RESOLVED_PLAN
    # The pressure sub-problem request is announced separately from the coupled
    # one, so an explicit pressure request can be compared with its resolution.
    assert records[1]["pressure_requested_preconditioner"] == "native_amg"
    assert records[1]["coupled_preconditioner"] == "mgr"
    assert records[1]["pressure_resolved"] is False
    assert records[1]["pressure_preconditioner"] == "none"


def test_linear_plan_accepts_explicit_requests_that_were_honored() -> None:
    evidence = MODULE.extract_n8_evidence([{"name": "t", "output": CTEST_OUTPUT}])
    plan = MODULE.audit_linear_plan(evidence)
    assert plan["status"] == "COMPLETE"
    assert plan["models_compared"] == 2
    assert plan["substitutions"] == []
    assert plan["models_without_plan"] == []
    assert plan["automatic_resolutions"] == []
    assert plan["policy"] == "explicit_request_must_be_honored"


def test_linear_plan_records_an_automatic_resolution_instead_of_hiding_it() -> None:
    # The segmented algorithms request `auto` for the pressure sub-problem and
    # the dispatcher resolves it to native AMG. That resolution is legitimate
    # and has to be reported rather than inferred from the request.
    segmented = {
        "run_configuration": [
            {**MODEL_CONFIGURATION, "algorithm": "SIMPLE/upwind/bounded",
             "preconditioner": "auto"},
        ],
        "resolved_plans": [
            {
                **RESOLVED_PLAN,
                "algorithm": "SIMPLE/upwind/bounded",
                "coupled_resolved": False,
                "coupled_preconditioner": "none",
                "pressure_resolved": True,
                "pressure_krylov": "cg",
                "pressure_preconditioner": "native_amg",
            }
        ],
    }
    plan = MODULE.audit_linear_plan(segmented)
    assert plan["status"] == "COMPLETE"
    assert plan["models_compared"] == 1
    assert plan["substitutions"] == []
    assert plan["automatic_resolutions"] == [
        {
            "algorithm": "SIMPLE/upwind/bounded",
            "subproblem": "pressure",
            "requested": "auto",
            "resolved": "native_amg",
        }
    ]


def test_linear_plan_flags_an_explicit_request_that_was_substituted() -> None:
    evidence = {
        "run_configuration": [
            {**MODEL_CONFIGURATION, "preconditioner": "mgr"},
        ],
        "resolved_plans": [RESOLVED_PLAN],
    }
    plan = MODULE.audit_linear_plan(evidence)
    assert plan["status"] == "VIOLATION"
    assert plan["substitutions"] == [
        {
            "algorithm": "COUPLED/BlockSchur/upwind/bounded",
            "subproblem": "coupled",
            "requested": "mgr",
            "resolved": "coupled_block_schur",
        }
    ]


def test_linear_plan_flags_a_configured_model_without_a_resolved_plan() -> None:
    plan = MODULE.audit_linear_plan(
        {"run_configuration": [MODEL_CONFIGURATION], "resolved_plans": []}
    )
    assert plan["status"] == "VIOLATION"
    assert plan["models_without_plan"] == ["COUPLED/BlockSchur/upwind/bounded"]
    assert plan["substitutions"] == []


def test_linear_plan_flags_a_resolution_without_an_announced_request() -> None:
    unresolved_request = {
        key: value
        for key, value in RESOLVED_PLAN.items()
        if key != "pressure_requested_preconditioner"
    }
    plan = MODULE.audit_linear_plan(
        {
            "run_configuration": [MODEL_CONFIGURATION],
            "resolved_plans": [
                {
                    **unresolved_request,
                    "pressure_resolved": True,
                    "pressure_krylov": "cg",
                    "pressure_preconditioner": "native_amg",
                }
            ],
        }
    )
    assert plan["status"] == "VIOLATION"
    assert plan["substitutions"] == [
        {
            "algorithm": "COUPLED/BlockSchur/upwind/bounded",
            "subproblem": "pressure",
            "requested": "undeclared",
            "resolved": "native_amg",
        }
    ]


def test_linear_plan_is_complete_when_an_explicit_request_is_honored() -> None:
    plan = MODULE.audit_linear_plan(
        {
            "run_configuration": [MODEL_CONFIGURATION],
            "resolved_plans": [RESOLVED_PLAN],
        }
    )
    assert plan["status"] == "COMPLETE"
    assert plan["models_compared"] == 1
    assert plan["automatic_resolutions"] == []
