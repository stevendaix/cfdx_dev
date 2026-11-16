"""Fail-closed tests for BFS native SST wall-model comparison."""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = (
    Path(__file__).resolve().parents[2] / "scripts/qualify_bfs_sst_wall_contract_v6.py"
)
spec = importlib.util.spec_from_file_location(
    "qualify_bfs_sst_wall_contract_v6", SCRIPT
)
assert spec and spec.loader
mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(mod)


@pytest.fixture
def native_case():
    return {
        "turbulence_model": "k_omega_sst",
        "boundary_conditions": [
            {
                "patch_name": "lowerWall",
                "type": "wall",
                "sst_scalar_conditions": [
                    {"field": "k", "type": "wall_k_zero"},
                    {"field": "omega", "type": "wall_omega_sst"},
                ],
            },
        ],
        "numerics": {
            "momentum_scheme": "first_order",
            "raw_settings": {"sst_wall_omega_factor": "6"},
        },
    }


def test_reference_physics_constants_are_exact():
    assert mod.REFERENCE["reattachment_x_over_h"] == 6.26
    assert mod.REFERENCE["reattachment_uncertainty_absolute"] == 0.10
    assert mod.REFERENCE["h_m"] == 0.0127
    assert mod.REFERENCE["u_inlet_m_s"] == 44.2


def test_native_sst_mismatch_blocks_reattachment(native_case):
    result = mod.inspect(native_case)
    assert not result["physical_qualification"]
    assert not result["openfoam_wall_functions_equivalent"]
    assert result["cfdx_wall_contract"]["wall_omega_factor"] == "6"
    assert not result["wall_yplus_qualified"]
    assert not result["wall_cf_qualified"]


def test_pressure_pass_does_not_become_sst_pass(native_case):
    log = (
        "PRESSURE_HILO_PHYSICAL_FLUX iteration=1 corrector=0 "
        "exact_relative=8.4e-14 stored_fp64_relative=8.5e-14 "
        "alpha_p=0.1 accepted_iteration=0\n"
    )
    result = mod.inspect(native_case, log, {"converged": False})
    assert result["pressure_substeps_pass_original_1e_11"]
    assert not result["solver_nonlinear_converged"]
    assert not result["physical_qualification"]


def test_solver_claim_converged_without_wall_evidence_still_blocked(native_case):
    result = mod.inspect(native_case, convergence={"converged": True})
    assert result["solver_nonlinear_converged"]
    assert not result["physical_qualification"]


def test_bad_flux_or_nonfinite_is_not_promoted(native_case):
    bad = "PRESSURE_HILO_PHYSICAL_FLUX exact_relative=3e-10 stored_fp64_relative=2e-10"
    assert not mod.inspect(native_case, bad)["pressure_substeps_pass_original_1e_11"]
    with pytest.raises(mod.BFSWallContractError):
        mod.inspect(native_case, bad.replace("3e-10", "nan"))


def test_missing_walls_is_rejected(native_case):
    native_case["boundary_conditions"] = []
    with pytest.raises(mod.BFSWallContractError, match="solid wall"):
        mod.inspect(native_case)


def test_native_hdf5_case_loading(tmp_path, native_case):
    h5py = pytest.importorskip("h5py")
    path = tmp_path / "mesh.cfdx.h5"
    with h5py.File(path, "w") as stream:
        stream.attrs["case_setup_json"] = json.dumps(native_case)
    assert mod.load_case(path) == native_case
