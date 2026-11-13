from pathlib import Path
import sys
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from analyze_bfs_homotopy_v48 import inspect

BASE=(
 'COUPLED_COMPOSITE_TRIAL_V47 iteration=1 retry=0 beta=1 carrier_linear_candidate_certified=1\n'
 'GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=0.5 p_continuity=0.1 k=0.01 omega=0.02 nut_closure=0 merit=0.5 converged=0 qualified=0\n'
 'COUPLED_COMPOSITE_TRIAL_V47 iteration=2 retry=1 beta=0.5 carrier_linear_candidate_certified=1\n'
 'COUPLED_RAU_HOMOTOPY_V48 iteration=2 retry=1 beta=0.5 previous_identity_l1=0 previous_identity_linf=0 effective_rAU_rederived=1 exact_physical_residual_oracle=0 committed=0\n'
 'COUPLED_NATIVE_FLUX_GAP_V48 iteration=2 retry=1 beta=0.5 relative_l1=0.03 absolute_linf=0.001 worst_face=7 native_operator_qualified=0 committed=0\n'
 'GLOBAL_NONLINEAR_COMMITTED_V13 iteration=2 U=0.4 p_continuity=0.11 k=0.01 omega=0.02 nut_closure=0 merit=0.4 converged=0 qualified=0\n'
 'CFDX production solver error: V14 independently rejected later state\n'
)
class TestV48Forensics(unittest.TestCase):
 def test_valid(self):
  o=inspect(BASE,2)
  self.assertEqual(o['status'],'ABORTED_UNQUALIFIED')
  self.assertEqual(o['native_gap_samples'],1)
  self.assertFalse(o['physical_qualification'])
 def test_reject_missing_gap(self):
  with self.assertRaises(ValueError):
   inspect('\n'.join(s for s in BASE.splitlines() if not s.startswith('COUPLED_NATIVE_FLUX_GAP_V48')))
 def test_reject_identity_drift(self):
  with self.assertRaises(ValueError):
   inspect(BASE.replace('previous_identity_l1=0 previous_identity_linf=0',
                        'previous_identity_l1=0.01 previous_identity_linf=0',1))
 def test_reject_false_certification(self):
  with self.assertRaises(ValueError):inspect(BASE.replace('converged=0 qualified=0','converged=1 qualified=1'))
 def test_reject_wrong_merit(self):
  with self.assertRaises(ValueError):inspect(BASE.replace('merit=0.4','merit=0.3'))
 def test_reject_noncontraction(self):
  with self.assertRaises(ValueError):inspect(BASE.replace('merit=0.4','merit=0.7').replace('U=0.4','U=0.7'))
 def test_reject_bad_attempt(self):
  with self.assertRaises(ValueError):inspect(BASE.replace('COUPLED_RAU_HOMOTOPY_V48 iteration=2 retry=1','COUPLED_RAU_HOMOTOPY_V48 iteration=2 retry=0'))
 def test_reject_missing_commit(self):
  with self.assertRaises(ValueError):inspect(BASE.replace('GLOBAL_NONLINEAR_COMMITTED_V13 iteration=2','GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1'))
 def test_reject_truncated_requirement(self):
  with self.assertRaises(ValueError):inspect(BASE,3)
if __name__=='__main__':unittest.main()
