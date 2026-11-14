"""Fail-closed mutation checks for the real V51 forensic evidence format."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'scripts'))
from analyze_bfs_fphys_v51 import audit


class Evidence(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.baseline = '\n'.join(
            f'GLOBAL_NONLINEAR_COMMITTED_V13 iteration={i} merit={1/i}'
            for i in range(1, 4))
        cls.candidate = ('FPHYS_CANDIDATE_V51 iteration=1 retry=0 cells=2 entries=12 '
                         'closed=1 replay_identical=1 flux_defect=1e-12 merit=0.2 '
                         'continuity=0.1 k=0.02 omega=0.03 raw_u=1 raw_k=2 raw_omega=3 '
                         'qualified=0 committed=0')
        cls.parity = ('FPHYS_PARITY_V51 iteration=1 max_normalized_gap=1e-14 '
                      'qualified=0')
        cls.cross = [
            ('FPHYS_CROSS_V51 iteration=1 column=%d Ux=1 Uy=1 Uz=0 '
             'continuity=1 k=1 omega=1 qualified=0') % i
            for i in range(6)]
        cls.observed = cls.baseline + '\n' + '\n'.join(
            [cls.candidate] * 4 + [cls.parity] + cls.cross)

    def test_valid(self):
        report = audit(self.observed, self.baseline, 3, True)
        self.assertFalse(report['qualified'])
        self.assertEqual(report['committed'], 3)

    def test_mutation_closed(self):
        self.reject(self.observed.replace('closed=1', 'closed=0'))

    def test_mutation_nonrepeat(self):
        self.reject(self.observed.replace('replay_identical=1', 'replay_identical=0'))

    def test_mutation_flux(self):
        self.reject(self.observed.replace('flux_defect=1e-12', 'flux_defect=1e-3'))

    def test_mutation_shape(self):
        self.reject(self.observed.replace('entries=12', 'entries=10'))

    def test_mutation_parity(self):
        self.reject(self.observed.replace('max_normalized_gap=1e-14', 'max_normalized_gap=1e-2'))

    def test_mutation_qualified(self):
        self.reject(self.observed.replace('qualified=0', 'qualified=1'))

    def test_mutation_committed(self):
        self.reject(self.observed.replace('committed=0', 'committed=1'))

    def test_mutation_missing_cross(self):
        self.reject(self.observed.replace(self.cross[5], ''))

    def test_mutation_sst(self):
        self.reject(self.observed.replace(' k=1 omega=1 qualified=0',
                                         ' k=0 omega=0 qualified=0'))

    def test_mutation_baseline(self):
        self.reject(self.observed.replace('iteration=3 merit=', 'iteration=3 merit=9'))

    def reject(self, trial):
        with self.assertRaises(ValueError):
            audit(trial, self.baseline, 3, True)


if __name__ == '__main__':
    unittest.main()
