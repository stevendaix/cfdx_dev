"""Mutation tests: do not let a V45 damped trial masquerade as physical convergence."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from compare_bfs_cycle_escape_v45 import compare

FP='COUPLED_MATRIX_FINGERPRINT_V40 A=aa b=bb x0=cc rows=608 nnz=7080\n'

def commit(i,merit):
    return (f'GLOBAL_NONLINEAR_COMMITTED_V13 iteration={i} U={merit} '
            f'p_continuity=0.00001 k=1e-12 omega=1e-13 nut_closure=0 '
            f'merit={merit} converged=0 qualified=0\n')

BASE=FP+commit(1,0.4)+commit(2,0.3)+commit(3,0.2)+commit(4,0.21)+(
    'Converged NO Iterations 4\n')
EVENT=('COUPLED_PHASE_ESCAPE_V45 iteration=3 retry=0 triggered=1 '
       'remaining=2 beta=0.5 linear_candidate_certified=1 '
       'updated_linear_solution_certified=0 physical_qualified=0\n')
TRIAL=FP+commit(1,0.4)+commit(2,0.3)+(
    'OUTER_FIELD_PERIOD_V43 iteration=2 period_four_candidate=1\n')+(
    'COUPLED_PERIOD4_FULL_V44 iteration=3 k_omega_averaged=1 '
    'nut_recomputed=1 strain_recomputed=1 committed=0 physical_qualified=0\n')+EVENT+commit(3,0.19)+commit(4,0.17)+(
    'COUPLED_PHASE_ESCAPE_V45 iteration=4 retry=0 triggered=0 remaining=1 '
    'beta=0.5 linear_candidate_certified=1 updated_linear_solution_certified=0 '
    'physical_qualified=0\nConverged NO Iterations 4\n')

class TestV45Evidence(unittest.TestCase):
    def test_good(self):
        out=compare(BASE,TRIAL,requested_beta=.5,steps=2)
        self.assertEqual(out['first_trigger'],3)
        self.assertFalse(out['physical_converged'])
        self.assertEqual(out['post_trigger_min_merit'],.17)
    def test_pre_history_mutation(self):
        with self.assertRaisesRegex(ValueError,'BEFORE'):
            compare(BASE,TRIAL.replace('iteration=2 U=0.3 p_continuity=0.00001',
                                      'iteration=2 U=0.3 p_continuity=0.00002'),
                    requested_beta=.5,steps=2)
    def test_stale_initial_matrix(self):
        with self.assertRaisesRegex(ValueError,'fingerprint'):
            compare(BASE,TRIAL.replace('A=aa','A=dd'),requested_beta=.5,steps=2)
    def test_false_certificate(self):
        with self.assertRaises(ValueError):
            compare(BASE,TRIAL.replace('updated_linear_solution_certified=0',
                                        'updated_linear_solution_certified=1'),
                    requested_beta=.5,steps=2)
    def test_no_complete_final(self):
        with self.assertRaises(ValueError):
            compare(BASE,TRIAL.replace('Converged NO Iterations 4',''),
                    requested_beta=.5,steps=2)
    def test_window_leak(self):
        with self.assertRaisesRegex(ValueError,'window'):
            compare(BASE,TRIAL.replace('iteration=4 retry=0',
                                       'iteration=5 retry=0'),requested_beta=.5,steps=2)
    def test_no_committed_orbit(self):
        with self.assertRaises(ValueError):
            compare(BASE,TRIAL.replace('period_four_candidate=1',
                                       'period_four_candidate=0'),requested_beta=.5,steps=2)
    def test_unattested_trigger(self):
        with self.assertRaises(ValueError):
            compare(BASE,TRIAL.replace('COUPLED_PERIOD4_FULL_V44',
                                       'COUPLED_PERIOD4_FAKE'),requested_beta=.5,steps=2)

if __name__=='__main__':
    unittest.main()
