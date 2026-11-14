"""Fail-closed tests for non-qualifying global-merit periodicity forensics."""
import pathlib
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT/'scripts'))
from bfs_outer_period_v42 import audit


def make_log(path: pathlib.Path, values):
    path.write_text(''.join(
        f'GLOBAL_NONLINEAR_COMMITTED_V13 iteration={i} U={m} '
        f'p_continuity=0.0001 k=1e-12 omega=1e-12 nut_closure=0 '
        f'merit={m} converged=0 qualified=0\n'
        for i, m in enumerate(values, 1)), encoding='utf-8')


class TestGlobalMeritForensics(unittest.TestCase):
    def test_realistic_four_phase(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/'a.log'
            vals=[0.001533, 0.003081, 0.004324, 0.002098]*5
            make_log(p, vals)
            result=audit(p)
            self.assertEqual(result['period'], 4)
            self.assertEqual(result['status'], 'REPEATING_UNCONVERGED_GLOBAL_MERIT')
            self.assertFalse(result['full_field_cycle_proven'])
            self.assertFalse(result['global_physical_convergence_proven'])

    def test_no_cycle_in_contraction(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/'a.log'
            make_log(p,[1.0/(2**i) for i in range(25)])
            self.assertEqual(audit(p)['status'],'NO_PERIOD_DETECTED')

    def test_reject_incomplete_history(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/'a.log'
            make_log(p,[0.002,0.004]*8)
            p.write_text(p.read_text().replace('iteration=7 ','iteration=8 ',1))
            with self.assertRaisesRegex(ValueError,'non-contiguous'):
                audit(p)

    def test_reject_fabricated_nonfinite(self):
        with tempfile.TemporaryDirectory() as d:
            p=pathlib.Path(d)/'a.log'
            make_log(p,[0.002,0.004]*8)
            p.write_text(p.read_text().replace('U=0.004','U=nan',1))
            with self.assertRaisesRegex(ValueError,'invalid global'):
                audit(p)

    def test_missing_file(self):
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'missing or empty'):
                audit(pathlib.Path(d)/'missing.log')


if __name__ == '__main__':
    unittest.main()
