import json
import sys
import unittest
from pathlib import Path

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'scripts'))
from diagnose_bfs_sst_v34 import compare,first_outer_attempt,sample


def fake_log(v34=False,merit=0.889):
    common=(
        'SST_JFNK_V22 iteration=10 selected=jfnk qualified=0\n'
        f'SST_TRANSPORT_DIAGNOSTIC iteration=10 k_local_relative={merit} '
        'omega_local_relative=0.2 qualified=0\n'
        'SST_ACCEPTANCE_GATES iteration=10 qualified_bfs=0 outer_update_ok=0 '
        'tolerance=1e-8 qualified=0\n'
    )
    if v34:
        common += ('SST_BRANCH_PRECOND_V34 iteration=10 evaluated=10 '
                   'corrected=0 skipped_nonsmooth=0 max_difference=0.002 '
                   'max_unilateral_gap=0.001 cells=341,329 qualified=0\n')
    return common+'CFDX nonlinear iteration 1 rejected; retry=1\n'


class V34ProtocolTests(unittest.TestCase):
    def test_failed_outer_state_is_not_promoted(self):
        report=compare(fake_log(),fake_log(True))
        self.assertEqual(report['status'],'NOT_CONVERGED_FIRST_OUTER_REJECTED')
        self.assertEqual(report['corrected_blocks'],0)
        self.assertEqual(report['merit_difference'],0.0)

    def test_missing_rejection_is_not_a_valid_baseline(self):
        with self.assertRaises(ValueError):
            compare(fake_log().split('CFDX nonlinear')[0],fake_log(True))

    def test_forged_physical_certificate_is_refused(self):
        with self.assertRaises(ValueError):
            compare(fake_log(),fake_log(True).replace('qualified_bfs=0','qualified_bfs=1'))

    def test_missing_branch_instrumentation_is_refused(self):
        with self.assertRaises(ValueError):
            compare(fake_log(),fake_log())

    def test_unbounded_residual_reassembly_is_refused(self):
        with self.assertRaises(ValueError):
            compare(fake_log(),fake_log(True).replace('evaluated=10','evaluated=500'))

    def test_threshold_tampering_is_refused(self):
        with self.assertRaises(ValueError):
            compare(fake_log(),fake_log(True).replace('tolerance=1e-8','tolerance=1e-6'))


if __name__=='__main__':
    unittest.main()
