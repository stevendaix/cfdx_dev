"""V11 evidence/physics-preservation tests; they do not pretend CFD convergence."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import h5py
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src/cfdx/python"))
from generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
from qualify_bfs_v11_numerical_matrix import (  # noqa: E402
    PROFILES,
    case_physics_signature,
    joint_returned_residual_summary,
    nested_pcd_summary,
    prepare_case,
    run_one,
)


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    root = tmp_path_factory.mktemp("bfs_v11_case")
    src = root / "reference.cfdx.h5"
    write_case(
        src,
        make_structured_mesh("smoke"),
        mesh_source="v11 fixture",
        convection_scheme="second_order_upwind",
        sst_inner_iterations=800,
    )
    return src


@pytest.mark.parametrize("profile", list(PROFILES))
def test_profiles_preserve_physical_case_and_hdf5(reference, tmp_path, profile):
    generated = tmp_path / (profile + ".cfdx.h5")
    detail = prepare_case(reference, generated, profile, 1, sst_inner_budget=800)
    assert detail["physics_signature_sha256"] == case_physics_signature(reference)
    assert detail["sst_inner_budget"] == 800
    with h5py.File(reference) as before, h5py.File(generated) as after:
        for key in before:
            assert np.array_equal(before[key][()], after[key][()]), key
        setup = json.loads(after.attrs["case_setup_json"])
        raw = setup["numerics"]["raw_settings"]
        assert raw["v8_global_nonlinear_residual"] == "true"
        assert raw["v7_post_sst_momentum_gate"] == "true"
        assert raw["sst_transport_max_iterations"] == "800"
        if profile.startswith("simple_full_face_hilo"):
            assert raw["v18_sparse_flux_pressure_jacobian"] == "true"
            assert raw["v19_full_pressure_flux_step"] == "true"
            assert raw["pressure_high_low_production"] == "true"
        if profile.endswith("sst_newton"):
            assert raw["sst_omega_newton_sink"] == "true"
            assert raw["sst_sequential_omega_first"] == "true"


def test_missing_second_order_is_refused(reference, tmp_path):
    wrong = tmp_path / "wrong.cfdx.h5"
    import shutil

    shutil.copy2(reference, wrong)
    with h5py.File(wrong, "r+") as f:
        s = json.loads(f.attrs["case_setup_json"])
        for item in s["numerics"]["selection"]["entries"]:
            if item.get("family") == "convection":
                item["configuration_key"] = "numerics.convection.upwind"
        f.attrs["case_setup_json"] = json.dumps(s)
    with pytest.raises(ValueError, match="bounded order-two"):
        prepare_case(
            wrong, tmp_path / "out.cfdx.h5", "simple_baseline", 1, sst_inner_budget=800
        )


def test_invalid_budget_is_rejected(reference, tmp_path):
    with pytest.raises(ValueError, match="positive SST"):
        prepare_case(
            reference,
            tmp_path / "out.cfdx.h5",
            "simple_baseline",
            1,
            sst_inner_budget=0,
        )


def test_run_one_timeout_reaps_process_group(reference, tmp_path):
    script = tmp_path / "simulate.py"
    script.write_text("#!/usr/bin/env python3\nimport time\ntime.sleep(5)\n")
    script.chmod(0o755)
    record = run_one(script, reference, tmp_path / "out", 1, 0.2, inner_trace=False)
    assert record["execution"] == "TIMEOUT_INCOMPLETE"
    assert record["exit_code"] is None
    assert record["phase_observation"]["proof_status"] == "NOT_QUALIFIED"


def test_guarded_profiles_do_not_qualify_unfinished_sst(reference, tmp_path):
    # Data creation is not a numerical/experimental certificate.
    detail = prepare_case(
        reference,
        tmp_path / "g.cfdx.h5",
        "block_local_sst_newton",
        1,
        sst_inner_budget=800,
    )
    assert detail["physical_qualification"] == "NOT_QUALIFIED"


def test_nested_pcd_original_operator_trace_summary():
    content = (
        "COUPLED_NESTED_PCD_FP status=1 original_true_relative_rhs=0.325 accepted=1\n"
        "COUPLED_NESTED_PCD_FP status=1 original_true_relative_rhs=0.0001 accepted=0\n"
        "COUPLED_NESTED_PCD_KP status=0 independent_true_relative_rhs=1e-9 accepted=1\n"
    )
    metric = nested_pcd_summary(content)
    assert metric["COUPLED_NESTED_PCD_FP"]["count"] == 2
    assert metric["COUPLED_NESTED_PCD_FP"]["max_true_relative_rhs"] == 0.325
    assert metric["COUPLED_NESTED_PCD_FP"]["above_one_e_minus_three"] == 1
    assert metric["COUPLED_NESTED_PCD_FP"]["rejected_actions"] == 1
    assert metric["COUPLED_NESTED_PCD_KP"]["count"] == 1
    with pytest.raises(ValueError, match="nonfinite"):
        nested_pcd_summary(
            "COUPLED_NESTED_PCD_FP original_true_relative_rhs=nan accepted=1"
        )
    with pytest.raises(ValueError, match="incomplete"):
        nested_pcd_summary("COUPLED_NESTED_PCD_KP accepted=1")


def test_anchor_profile_preserves_all_physical_hdf5(reference, tmp_path):
    source = tmp_path / "anchor.cfdx.h5"
    prepare_case(
        reference, source, "pcd_fp_recirculation_anchor", 1, sst_inner_budget=800
    )
    with h5py.File(source) as f:
        raw = json.loads(f.attrs["case_setup_json"])["numerics"]["raw_settings"]
        assert raw["pcd_fp_recirculation_anchor_experimental"] == "true"
        assert raw["sst_transport_max_iterations"] == "800"


def test_auxiliary_fp_does_not_masquerade_as_original_inverse():
    trace = (
        "COUPLED_NESTED_PCD_FP status=0 original_true_relative_rhs=0.11 "
        "auxiliary_anchored=1 auxiliary_true_relative_rhs=6e-9 "
        "auxiliary_anchors=2 accepted=1\n"
    )
    metrics = nested_pcd_summary(trace)
    assert metrics["COUPLED_NESTED_PCD_FP"]["max_true_relative_rhs"] == 0.11
    assert metrics["fp_auxiliary_anchored"]["anchored_actions"] == 1
    assert metrics["fp_auxiliary_anchored"]["max_auxiliary_true_relative_rhs"] == 6e-9
    with pytest.raises(ValueError, match="violated"):
        nested_pcd_summary(trace.replace("0.11", "1.0"))
    with pytest.raises(ValueError, match="violated"):
        nested_pcd_summary(trace.replace("6e-9", "1e-3"))
    with pytest.raises(ValueError, match="incomplete"):
        nested_pcd_summary(trace.replace(" auxiliary_anchors=2", ""))


def test_joint_returned_residual_gate_does_not_promote_local_sst_or_pressure():
    lines = (
        "PRESSURE_HILO_PRODUCTION iteration=1 final_face_relative=1e-16 "
        "requested_tolerance=1e-11 true_face_converged=1\n"
        "BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_END committed=1 converged=1\n"
        "GLOBAL_NONLINEAR_V8 iteration=1 U=0.5275 p_continuity=0.00176 "
        "k=1e-13 omega=4e-14 nut_closure=0 merit=0.5275 accepted=0 qualified=0\n"
    )
    result = joint_returned_residual_summary(lines)
    assert result["sst_inner_completed_observations"] == 1
    assert result["true_face_pressure_converged_observations"] == 1
    assert result["global_qualified_observations"] == 0
    assert result["global_accepted_observations"] == 0
    assert result["first_returned_state"]["worst_field"] == "U"
    assert result["best_observed_merit"] == 0.5275


def test_joint_gate_fails_closed_on_incomplete_or_forged_evidence():
    assert joint_returned_residual_summary("")["independent_global_observations"] == 0
    line = (
        "GLOBAL_NONLINEAR_V8 iteration=1 U=0.2 p_continuity=0.01 "
        "k=1e-13 omega=4e-14 nut_closure=0 merit=0.2 accepted=0 qualified=0"
    )
    for bad in (
        line.replace("U=0.2", "U=nan"),
        line.replace("merit=0.2", "merit=0.01"),
        line.replace("accepted=0 qualified=0", "accepted=0 qualified=1"),
        line.replace("omega=4e-14 ", ""),
    ):
        with pytest.raises(ValueError):
            joint_returned_residual_summary(bad)
    assert (
        joint_returned_residual_summary(
            "BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_BEGIN budget=800"
        )["global_qualified_observations"]
        == 0
    )
