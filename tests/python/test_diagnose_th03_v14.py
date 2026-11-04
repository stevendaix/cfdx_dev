"""V14 diagnostic gates stay observational and fail closed."""
from pathlib import Path

import pytest

from scripts.diagnose_th03_v14 import analyse
from scripts.generate_th03_openfoam_blockmesh import make_case_setup, make_structured_mesh, write_case


def test_v14_rejects_logs_without_th03_evidence(tmp_path: Path):
    log = tmp_path / 'noise.log'
    log.write_text('Converged YES Iterations 99\n', encoding='utf-8')
    with pytest.raises(ValueError, match='no CFDX TH-03'):
        analyse(log)


def test_v14_separates_krylov_from_returned_state_and_qualification(tmp_path: Path):
    log = tmp_path / 'coupled.log'
    log.write_text('''INCOMPRESSIBLE_ITER iter=1 pressure_relative_residual=3.8e-11 continuity_normalized=0.01
SST_TRANSPORT_DIAGNOSTIC iteration=30 k_returned_nonlinear_raw=0.012 omega_returned_nonlinear_raw=15 k_local_relative=0.22 omega_local_relative=0.11 k_local_worst_cell=3 omega_local_worst_cell=8 floor_k=0 floor_omega=0 finite_events=0 k_linear_status=converged omega_linear_status=converged k_linear_relative_before_relaxation=1e-12 k_linear_relative_after_relaxation=0.1 omega_linear_relative_before_relaxation=2e-12 omega_linear_relative_after_relaxation=0.2
SST_PROVISIONAL_OUTER iteration=1 qualified=0
Converged NO Iterations 1
''', encoding='utf-8')
    result = analyse(log)
    assert result['status'] == 'DIAGNOSTIC_ONLY_NOT_QUALIFIED'
    assert result['pressure']['last_pressure_relative_residual'] == 3.8e-11
    assert result['sst']['all_provisional_steps_unqualified']
    assert result['sst']['last_linear_contract']['k_linear_relative_after_relaxation'] == 0.1
    assert result['sst']['last_nonlinear_sample']['omega_returned_nonlinear_raw'] == 15.


def test_v14_rejects_invalid_picard_damping():
    for val in (0, -0.1, float('nan'), float('inf'), 1.01):
        with pytest.raises(ValueError, match='Picard damping'):
            make_case_setup(2, 0.25, 0.1, sst_picard_damping=val)


def test_v14_controls_unchanged_by_default_and_exported_only_on_request(tmp_path: Path):
    original = make_case_setup(2, 0.25, 0.1).numerics.raw_settings
    assert original['sst_scale_aware_convergence'] == 'false'
    assert 'sst_outer_provisional' not in original
    assert 'sst_picard_damping' not in original
    modified = make_case_setup(2, 0.25, 0.1,
        sst_provisional=True, sst_picard_damping=0.25,
        sst_relative_gate=True).numerics.raw_settings
    assert modified['sst_outer_provisional'] == 'enabled'
    assert modified['sst_picard_damping'] == '0.25'
    assert modified['sst_scale_aware_convergence'] == 'true'
    outfile = tmp_path / 'case.h5'
    write_case(outfile, make_structured_mesh('smoke'), mesh_source='V14 test',
        max_iterations=2, sst_provisional=True, sst_picard_damping=0.25,
        sst_relative_gate=True)
    assert outfile.exists()