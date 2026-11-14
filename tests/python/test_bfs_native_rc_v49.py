"""No external dependencies: mutation tests of fail-closed V49 evidence."""
import pathlib
import sys
import unittest

ROOT=pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from analyze_bfs_native_rc_v49 import validate

LINE=('NATIVE_RC_ORACLE_V49 iteration=1 retry=0 closed=1 '
      'fixed_point_iterations=24 fixed_point_defect=1.8e-12 '
      'repeated_same_state=0 native_flux_scale_linf=0.0075 '
      'U=1.25 continuity=0.008 k=0.063 omega=0.052 '
      'carrier_flux_gap_l1=0.0048 carrier_flux_gap_linf=0.0016 '
      'worst_face=3 full_six_field_oracle_qualified=0 '
      'physical_converged=0 committed=0 qualified=0')
COMMIT=('GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=1.89 '
        'p_continuity=0.012 k=1e-11 omega=1e-13 '
        'nut_closure=0 merit=1.89 converged=0 qualified=0')
LOG=LINE+'\n'+COMMIT+'\n'

class V49Evidence(unittest.TestCase):
    def test_valid_diagnostic_never_qualified(self):
        r=validate(LOG,LOG)
        self.assertEqual(r['status'],'DIAGNOSTIC_ONLY_PHYSICAL_UNQUALIFIED')
        self.assertFalse(r['qualified'])
        self.assertEqual(r['native_flux_fixed_point_closed'],1)
    def test_missing_oracle_rejected(self):
        with self.assertRaises(ValueError):validate(COMMIT)
    def test_no_commit_rejected(self):
        with self.assertRaises(ValueError):validate(LINE)
    def test_fake_physical_convergence_rejected(self):
        with self.assertRaises(ValueError):validate(LOG.replace('physical_converged=0','physical_converged=1'))
    def test_wrong_operator_certification_rejected(self):
        with self.assertRaises(ValueError):validate(LOG.replace('full_six_field_oracle_qualified=0','full_six_field_oracle_qualified=1'))
    def test_bad_fixed_point_rejected(self):
        with self.assertRaises(ValueError):validate(LOG.replace('1.8e-12','1e-3'))
    def test_nan_rejected(self):
        with self.assertRaises(ValueError):validate(LOG.replace('k=0.063','k=nan'))
    def test_duplicate_trial_rejected(self):
        with self.assertRaises(ValueError):validate(LINE+'\n'+LOG)
    def test_reference_mutation_rejected(self):
        with self.assertRaises(ValueError):validate(LOG,LOG.replace('U=1.89','U=0.89'))
    def test_invalid_replay_flag_rejected(self):
        with self.assertRaises(ValueError):
            validate(LOG.replace('repeated_same_state=0','repeated_same_state=2'))
    def test_missing_flux_scale_rejected(self):
        with self.assertRaises(ValueError):
            validate(LOG.replace('native_flux_scale_linf=0.0075 ',''))
    def test_missing_independent_residual_rejected(self):
        with self.assertRaises(ValueError):validate(LOG.replace('omega=0.052 ',''))
    def test_failure_explicit_not_qualified(self):
        failure=('NATIVE_RC_ORACLE_V49 iteration=1 retry=0 closed=0 '
                 'fixed_point_iterations=32 fixed_point_defect=0.4 '
                 'full_six_field_oracle_qualified=0 physical_converged=0 '
                 'committed=0 qualified=0')
        r=validate(failure+'\n'+COMMIT, min_closed=0)
        self.assertEqual(r['native_flux_fixed_point_closed'],0)

if __name__=='__main__':unittest.main()
