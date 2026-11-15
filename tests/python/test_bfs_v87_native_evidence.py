@@ -0,0 +1,41 @@
"""Falsification tests: an early screen cannot become a certification."""
import unittest
from scripts.audit_bfs_v87_native import audit
OK='BFS_EARLY_FIELD_TRIAL_V87 iteration=3 continue_sst=1 max_speed_over_inlet=1.1 max_speed_cell=10 reason=velocity_preflight_only_not_certified qualified=0'
BAD='BFS_EARLY_FIELD_TRIAL_V87 iteration=3 continue_sst=0 max_speed_over_inlet=3.356 max_speed_cell=7406 reason=velocity_envelope qualified=0'
SST='BFS_PHASE_V10 iteration=3 phase=SST_TRANSPORT_BEGIN budget=10 target=1e-8'
ROLL='CFDX nonlinear iteration 3 rejected; retry=1 alpha_u=1.25e-1 alpha_p=5e-2 reason=V87 pre-SST screen: velocity_envelope'
COMMIT='GLOBAL_PROGRESS_COMMITTED_V85 iteration=3 committed=3 k_debt=1 omega_debt=2 converged=0 qualified=0'
POST='BFS_FIELD_TRUST_TRIAL_V86 iteration=3 accepted=1 reason=provisional_step_within_field_trust_region max_speed_over_inlet=1.2 k_max_over_inlet_squared=0.001 qualified=0'
class TestV87Evidence(unittest.TestCase):
 def test_valid_rejection(self):
  r=audit(BAD+'\n'+ROLL)
  self.assertEqual(r['early_rejections_with_rollback'],1)
  self.assertEqual(r['early_passes_to_sst'],0)
 def test_reject_running_sst(self):
  with self.assertRaisesRegex(ValueError,'SST ran after rejected'):
   audit(BAD+'\n'+SST)
 def test_valid_retry_and_commit(self):
  r=audit('\n'.join((BAD,ROLL,OK,SST,POST,COMMIT)))
  self.assertEqual(r['committed_iterations'],[3]);self.assertEqual(r['early_passes_to_sst'],1)
 def test_no_false_certification(self):
  with self.assertRaisesRegex(ValueError,'falsely certified'):
   audit(BAD.replace('qualified=0','qualified=1')+'\n'+ROLL)
 def test_above_limit(self):
  with self.assertRaisesRegex(ValueError,'Unsafe'):
   audit(OK.replace('max_speed_over_inlet=1.1','max_speed_over_inlet=3.2'))
 def test_below_limit_rejected(self):
  with self.assertRaisesRegex(ValueError,'Incoherent'):
   audit(BAD.replace('3.356','1.5')+'\n'+ROLL)
 def test_fake_certification_post(self):
  with self.assertRaisesRegex(ValueError,'False V86'):
   audit(POST.replace('qualified=0','qualified=1'))
 def test_fake_commit(self):
  with self.assertRaisesRegex(ValueError,'False final'):
   audit(COMMIT.replace('converged=0','converged=1'))
 def test_missing_rollback(self):
  with self.assertRaisesRegex(ValueError,'Previous early trial'):
   audit(BAD+'\n'+OK)
 def test_incomplete_timeout(self):
  r=audit(OK);self.assertTrue(r['incomplete_last_preflight']);self.assertFalse(r['physical_convergence_certified'])
if __name__=='__main__':unittest.main()
