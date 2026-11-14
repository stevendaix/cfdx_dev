import sys
from pathlib import Path
import pytest
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from diagnose_bfs_sst_v29 import first_attempt


def trace(active='1', available='1', selected='jfnk',
          candidate='0.895', l2='0.17', old='0.926', new='0.890',
          qualified='0', cfl_after='0.5', target='1e-8', stall='0.82'):
    return (f'SST_JFNK_V22 iteration=10 v29_active={active} '
            f'v29_probe={active} v29_available={available} '
            'v29_reason=dual_descending v29_cfl_before=0.5 '
            f'v29_cfl_after={cfl_after} v29_shifted=1355 v29_wall_skip=69 '
            f'v29_linf={candidate} v29_l2={l2} '
            f'v26_l2_before=0.178 old_merit={old} new_merit={new} '
            f'selected={selected} accepted=1 qualified={qualified}\n'
            'SST_FP64_PRECISION_STALL_V14 iteration=27 identical_steps=3 '
            f'k_local_relative={stall} target={target} accepted=0 qualified=0\n'
            'SST_JFNK_V22 iteration=10 v29_active=0 selected=picard\n')


def test_valid_proposal_is_not_convergence():
    r = first_attempt(trace(), expected_active=True)
    assert r['v29_eligible'] and not r['v29_selected']
    assert r['first_stagnation_iteration'] == 27
    assert not r['sst_converged'] and not r['outer_state_accepted']


def test_inactive_path_is_observational():
    r = first_attempt(trace(active='0', available='0'), expected_active=False)
    assert not r['v29_eligible'] and not r['v29_active']


def test_rejects_candidate_failing_either_gate():
    with pytest.raises(ValueError):
        first_attempt(trace(candidate='0.94'), expected_active=True)
    with pytest.raises(ValueError):
        first_attempt(trace(l2='0.19'), expected_active=True)


def test_rejects_false_qualification_or_threshold_relaxation():
    with pytest.raises(ValueError):
        first_attempt(trace(qualified='1'), expected_active=True)
    with pytest.raises(ValueError):
        first_attempt(trace(target='1e-5'), expected_active=True)


def test_rejects_noncommitted_cfl_change_and_unknown_branch():
    with pytest.raises(ValueError):
        first_attempt(trace(cfl_after='2'), expected_active=True)
    with pytest.raises(ValueError):
        first_attempt(trace(selected='invented'), expected_active=True)


def test_rejects_fake_or_missing_stagnation():
    with pytest.raises(ValueError):
        first_attempt(trace(stall='1e-10'), expected_active=True)
    with pytest.raises(ValueError):
        first_attempt(trace().split('SST_FP64_PRECISION_STALL_V14')[0], expected_active=True)


def test_cannot_select_unavailable_or_inactive_candidate():
    with pytest.raises(ValueError):
        first_attempt(trace(available='0', selected='sst_pseudo_transient_v29'), expected_active=True)
    with pytest.raises(ValueError):
        first_attempt(trace(active='0', available='1'), expected_active=False)
