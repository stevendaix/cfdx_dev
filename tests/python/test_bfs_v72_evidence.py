@@ -0,0 +1,38 @@
#!/usr/bin/env python3
import pathlib
import sys
ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from audit_bfs_v72_fixed_point import audit

HISTORY=(
    'GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=23 committed=0',
    'GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=22 accepted=1',
    'BFS_SST_RC_PICARD_V69 iteration=23 private_transport=1',
    'INCOMPRESSIBLE_ITER iter=23 momentum=0.1',
)
A='BFS_JOINT_FIXED_POINT_V72 iteration=23 retry=8 cycle=1 gamma=0.03125 pressure_fraction=0.25 mass_before=0.0008 mass_after=0.0006 omega_after=0.0007 merit_ratio=0.93 worst_growth=0.93 rAU_change=4.0 nut_change=0.001 phi_change=0.0001 evaluated=1 v65_guards=1 mass_improved=1 selected=1 rejected_reason=none private_trial=1 committed=0 qualified=0'
C='BFS_JOINT_CERTIFICATE_V72 iteration=23 cycle=1 mass=0.0006 boundary_net=0.0 face_identity=0.0 continuity_face_gap=0.0 rc_fixed_point=0.0 momentum_scaled=0.1 k_scaled=0.1 omega_scaled=0.01 replay=1 mass_closed=0 six_blocks_closed=0 linear_momentum_reissued=0 linear_sst_reissued=0 commit_authorized=0 private_trial=1 committed=0 qualified=0'
S='BFS_JOINT_SUMMARY_V72 iteration=23 attempted_cycles=1 accepted_private_cycles=1 initial_mass=0.0008 final_mass=0.0006 mass_closed=0 all_physical_blocks_closed=0 status=cycle_budget_exhausted_without_full_certificate commit_authorized=0 committed=0 qualified=0'
BASE='\n'.join(HISTORY)+'\n'
LOG=BASE+'\n'.join((A,C,S))+'\n'
result=audit(BASE,LOG)
assert result['accepted_private_cycles']==1
bad=[
    LOG.replace('committed=0 qualified=0','committed=1 qualified=0',1),
    LOG.replace('selected=1','selected=0'),
    LOG.replace('face_identity=0.0','face_identity=0.001'),
    LOG.replace('continuity_face_gap=0.0','continuity_face_gap=0.001'),
    LOG.replace('rc_fixed_point=0.0','rc_fixed_point=0.001'),
    LOG.replace('linear_momentum_reissued=0','linear_momentum_reissued=1'),
    LOG.replace('linear_sst_reissued=0','linear_sst_reissued=1'),
    LOG.replace('mass_closed=0','mass_closed=1'),
    LOG.replace('mass=0.0006','mass=0.0007'),
    LOG.replace('merit_ratio=0.93','merit_ratio=nan'),
    LOG.replace('GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=22','GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=21'),
]
for j, falsified in enumerate(bad):
    try:audit(BASE,falsified)
    except (ValueError,KeyError):pass
    else:raise AssertionError('falsification '+str(j)+' accepted')
print('V72_EVIDENCE_FALSIFICATION_PASS injected_failures=',len(bad))
