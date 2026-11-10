import sys
from pathlib import Path
import pytest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
from qualify_bfs_committed_v13 import parse

TRIAL="GLOBAL_NONLINEAR_V8 iteration=1 U=0.5 p_continuity=0.002 k=0.0001 omega=0.001 nut_closure=0 merit=0.5 accepted=0 qualified=0"
COMMIT="GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=0.5 p_continuity=0.002 k=0.0001 omega=0.001 nut_closure=0 merit=0.5 converged=0 qualified=0"

def test_tentative_without_commit_is_not_accepted():
    v=parse(TRIAL)
    assert v['committed_count']==0 and v['uncommitted_trial_count']==1
    assert not v['physical_qualification']

def test_commit_counted_and_not_a_physics_pass():
    v=parse(TRIAL+'\n'+COMMIT)
    assert v['committed_count']==1 and v['uncommitted_trial_count']==0
    assert not v['all_committed_converged']

def test_rejected_and_retried_trial_not_counted():
    v=parse(TRIAL+'\nCFDX nonlinear iteration 1 rejected; retry=1\n'+TRIAL+'\n'+COMMIT)
    assert v['committed_count']==1 and v['uncommitted_trial_count']==1
    assert v['rejected_retry_markers']==1

def test_forged_commit_rejected():
    with pytest.raises(ValueError,match='without a preceding'):
        parse(COMMIT)
    with pytest.raises(ValueError,match='differs'):
        parse(TRIAL+'\n'+COMMIT.replace('U=0.5','U=0.7').replace('merit=0.5','merit=0.7'))

def test_duplicate_commit_rejected():
    with pytest.raises(ValueError,match='duplicate/reordered'):
        parse(TRIAL+'\n'+COMMIT+'\n'+TRIAL+'\n'+COMMIT)

def test_nonfinite_tentative_rejected():
    with pytest.raises(ValueError,match='nonfinite'):
        parse(TRIAL.replace('U=0.5','U=nan'))
