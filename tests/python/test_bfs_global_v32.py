"""V32 regression: do not confuse inner SST progress with outer CFDX commits."""
import importlib.util
from pathlib import Path
import pytest
p=Path(__file__).resolve().parents[2]/'scripts/audit_bfs_global_v32.py'
spec=importlib.util.spec_from_file_location('bfs_global_v32',p)
x=importlib.util.module_from_spec(spec);spec.loader.exec_module(x)

def line(i,U,p,k,w,nut,conv=0):
    vals=[U,p,k,w,nut]
    return (f'GLOBAL_NONLINEAR_COMMITTED_V13 iteration={i} U={U} p_continuity={p} '
            f'k={k} omega={w} nut_closure={nut} merit={max(vals)} converged={conv} qualified=0\n')

def test_sst_inner_progress_without_outer_commit_is_not_qualified():
    r=x.audit('SST_JFNK_V22 iteration=10 selected=picard accepted=1\n'
              'CFDX nonlinear iteration 1 rejected; retry=1\n')
    assert r['committed_outer_states']==0
    assert r['sst_internal_accepts']==1
    assert r['outer_rejections']==1
    assert not r['qualified']

def test_three_commits_without_final_convergence_is_still_unqualified():
    s=''.join(line(i,v,.01*v,v,v,v) for i,v in ((1,.2),(2,.1),(3,.05)))
    r=x.audit(s)
    assert r['committed_outer_states']==3 and not r['qualified']

def test_three_strict_commits_and_all_gates_can_qualify():
    s=line(1,.2,.1,.1,.1,.1)+line(2,.1,.05,.05,.05,.05)+line(3,1e-10,1e-11,1e-10,1e-10,1e-10,1)
    assert x.audit(s)['qualified']

def test_false_converged_flag_fails_closed():
    with pytest.raises(ValueError,match='false physical'):
        x.audit(line(1,.1,.01,.1,.1,.1,1))

def test_unordered_or_noncontracting_commits_fail_closed():
    with pytest.raises(ValueError,match='did not contract'):
        x.audit(line(1,.1,.01,.1,.1,.1)+line(2,.11,.01,.11,.11,.11))
    with pytest.raises(ValueError,match='non-increasing'):
        x.audit(line(2,.1,.01,.1,.1,.1)+line(1,.05,.01,.05,.05,.05))
