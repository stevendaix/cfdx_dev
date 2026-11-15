import importlib.util
from pathlib import Path
import pytest

P=Path(__file__).resolve().parents[2]/'scripts'/'audit_bfs_v88_native.py'
spec=importlib.util.spec_from_file_location('audit_bfs_v88_native',P)
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)

def one(ratio,allowed=1, iteration=1,retry=0):
    reason=('delta_velocity_preflight_only_not_certified' if allowed
            else 'local_velocity_increment_envelope')
    return (f'BFS_LOCAL_STEP_TRIAL_V88 iteration={iteration} retry={retry} '
            f'continue_sst={allowed} max_delta_u_over_inlet={ratio} '
            'delta_u_limit=0.5 worst_cell=7406 delta_pressure_span_over_qin=12.1 '
            f'alpha_u=0.125 alpha_p=0.05 reason={reason} qualified=0 committed=0\n')

def test_rejected_trial_stays_private():
    log=one(1.9,0)+'CFDX nonlinear iteration 1 rejected; retry=1\n'+one(0.49,1)+\
        'BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_BEGIN\n'+\
        'GLOBAL_PROGRESS_COMMITTED_V85 iteration=1 committed=1 converged=0 qualified=0\n'
    rep=module.audit(log)
    assert rep['pre_sst_rejected']==1
    assert rep['accepted_outer_count']==1
    assert not rep['finally_converged']

@pytest.mark.parametrize('corrupt',[ 'committed=1','qualified=1','continue_sst=1',
                                     'max_delta_u_over_inlet=nan','delta_u_limit=nan'])
def test_corrupt_evidence_rejected(corrupt):
    original=one(1.9,0)
    old = ('committed=0' if corrupt.startswith('committed=') else
           'qualified=0' if corrupt.startswith('qualified=') else
           'continue_sst=0' if corrupt.startswith('continue_sst=') else
           'max_delta_u_over_inlet=1.9' if corrupt.startswith('max_delta_') else
           'delta_u_limit=0.5')
    with pytest.raises(ValueError):module.audit(original.replace(old,corrupt))

def test_missing_required_evidence_rejected():
    with pytest.raises(ValueError):module.audit(one(1.9,0).replace('worst_cell=7406 ',''))

def test_sst_on_rejected_trial_rejected():
    with pytest.raises(ValueError):module.audit(one(1.9,0)+'BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_BEGIN\n')

def test_false_commit_sequence_rejected():
    with pytest.raises(ValueError):module.audit('GLOBAL_PROGRESS_COMMITTED_V85 iteration=2 converged=0 qualified=0\n')