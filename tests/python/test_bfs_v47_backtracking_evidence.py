"""Mutation tests: no false converged CFDX V47 evidence on partial logs."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bfs_composite_v47 import audit

def need_raises(log, **kwargs):
    try: audit(log, **kwargs)
    except ValueError: return
    raise AssertionError('invalid V47 evidence was accepted')

trial = ('COUPLED_COMPOSITE_TRIAL_V47 iteration={i} retry={r} beta={b} '
         'search_kind={kind} carrier_linear_candidate_certified=1 damped_linear_solution_certified=0 '
         'flux_rederived=1 sst_transport_pending=1 physical_qualified=0 committed=0\n')
commit = ('GLOBAL_NONLINEAR_COMMITTED_V13 iteration={i} U={u} p_continuity=0.002 '
          'k=1e-11 omega=2e-11 nut_closure=0 merit={u} converged=0 qualified=0\n')
v8 = ('GLOBAL_NONLINEAR_V8 iteration={i} U={u} p_continuity=0.002 k=1e-11 omega=2e-11 '
      'nut_closure=0 merit={u} accepted=0 qualified=0\n')
v14 = ('GLOBAL_CONTRACTION_TRIAL_V14 iteration={i} old={old} new={new} '
       'ratio=0.5 allowed=1 qualified=0\n')
valid = (trial.format(i=1,r=0,b=1.0,kind='geometric') +
         v14.format(i=1,old=0,new=0.4) + v8.format(i=1,u=0.4) + commit.format(i=1,u=0.4) +
         trial.format(i=2,r=0,b=1.0,kind='geometric') + trial.format(i=2,r=1,b=0.5,kind='geometric') +
         trial.format(i=2,r=2,b=0.25,kind='geometric') + trial.format(i=2,r=3,b=0.125,kind='geometric') +
         trial.format(i=2,r=4,b=0.0625,kind='geometric') +
         'GLOBAL_JOINT_V13 iteration=2 growth=0.5 accepted=1 qualified=0\n' +
         v14.format(i=2,old=0.4,new=0.1) + v8.format(i=2,u=0.1) + commit.format(i=2,u=0.1) +
         'Converged NO Iterations 2\n')
report = audit(valid,require_deep_retry=True)
assert report['committed_outer_steps']==2
assert report['deepest_retry']==4
assert report['physical_converged'] is False
need_raises(valid.replace('beta=0.0625','beta=0.125'))
need_raises(valid.replace('damped_linear_solution_certified=0','damped_linear_solution_certified=1'))
need_raises(valid.replace('merit=0.1','merit=0.3'))
need_raises(valid.replace('U=0.1','U=0.5'))
need_raises(valid.replace('committed=0','committed=1'))
need_raises(valid.replace('Converged NO','Converged YES'))
need_raises(valid.replace('GLOBAL_NONLINEAR_V8 iteration=2 U=0.1','GLOBAL_NONLINEAR_V8 iteration=2 U=0.2'))
need_raises(valid.replace('old=0.4 new=0.1','old=0.4 new=0.3'))
need_raises(valid.replace('growth=0.5 accepted=1','growth=0.5 accepted=0'))
need_raises(valid.replace('Converged NO Iterations 2','Converged NO Iterations 3'),require_complete_log=True)
need_raises(valid.replace('retry=3 beta=0.125','retry=3 beta=0.0625'))
print('BFS_V47_BACKTRACKING_EVIDENCE_PASS mutations=11 no_physical_certificate=1')
