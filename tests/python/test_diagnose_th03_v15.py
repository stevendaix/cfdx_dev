"""V15 diagnostics prove strict non-qualification and monitor arithmetic gaps."""
from pathlib import Path

import pytest

from scripts.diagnose_th03_v15 import analyse


def test_v15_parser_detects_residual_gap_insufficient_to_explain_failure(tmp_path: Path):
    log = tmp_path / 'v15.log'
    log.write_text('''INCOMPRESSIBLE_ITER iter=1 pressure_relative_residual=1e-11 continuity_normalized=0.1
SST_TRANSPORT_DIAGNOSTIC iteration=800 k_returned_nonlinear_raw=1e-9 omega_returned_nonlinear_raw=1e-9 k_local_relative=1.57e-8 omega_local_relative=1e-9 k_local_compensated=1.569e-8 omega_local_compensated=1e-9 k_local_roundoff_gap=1e-11 omega_local_roundoff_gap=1e-12 compensated_gate=0 k_local_worst_cell=210 omega_local_worst_cell=222
SST_WALL_DISCRETE_AUDIT iteration=1 faces=12 yplus_viscous_proxy_max=2.3 qualified=0
SST_NONLINEAR_LINESEARCH iteration=1 alpha=0.5 old_merit=0.1 new_merit=0.08 accepted=1 qualified=0
SST_NONLINEAR_LINESEARCH iteration=2 alpha=0 old_merit=0.08 new_merit=0.08 accepted=0 qualified=0
Converged NO Iterations 1
''', encoding='utf-8')
    result = analyse(log)
    assert result['v15']['precision_samples'] == 1
    assert result['v15']['last_precision_sample']['roundoff_alone_cannot_explain_k_gate']
    assert result['v15']['wall_audit_samples'] == 1
    assert result['v15']['line_search_accepted_samples'] == 1
    assert result['v15']['line_search_rejected_samples'] == 1
    assert not result['v15']['line_search_contract_failures']
    assert result['v15']['physical_qualification'] == 'NOT_QUALIFIED'


def test_v15_parser_flag_false_linesearch_claim(tmp_path: Path):
    log = tmp_path / 'bad.log'
    log.write_text('''SST_NONLINEAR_LINESEARCH iteration=1 alpha=1 old_merit=0.1 new_merit=0.2 accepted=1
''', encoding='utf-8')
    with pytest.raises(ValueError, match='no CFDX TH-03'):
        analyse(log)
    log.write_text('''INCOMPRESSIBLE_ITER iter=1 pressure_relative_residual=1e-10
SST_NONLINEAR_LINESEARCH iteration=1 alpha=1 old_merit=0.1 new_merit=0.2 accepted=1
''', encoding='utf-8')
    assert analyse(log)['v15']['line_search_contract_failures'] == [
        'accepted-without-strict-reduction']


def test_v15_parser_rejects_invalid_tolerance(tmp_path: Path):
    with pytest.raises(ValueError):
        analyse(tmp_path / 'missing.log', sst_tolerance=0.0)