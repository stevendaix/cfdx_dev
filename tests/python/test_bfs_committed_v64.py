@@ -0,0 +1,32 @@
import unittest
from scripts.audit_bfs_committed_v64 import audit

def stream(n=3):
    out=[]
    for i in range(1,n+1):
        out.extend([
            f'BFS_PHASE_V10 iteration={i} phase=COUPLED_LINEAR_END status=0 true_relative=1e-12',
            f'BFS_PHASE_V10 iteration={i} phase=SST_TRANSPORT_END converged=1 k_residual=1e-11 omega_residual=2e-11',
            f'POST_SST_MOMENTUM iteration={i} relative_nut_change=0.02 old_nut_operator=0.2 refreshed_nut_operator=0.25',
            f'GLOBAL_NONLINEAR_COMMITTED_V13 iteration={i} U={1 if i==1 else 0.5 if i==2 else 0.6} '
            f'p_continuity=0.1 k=1e-11 omega=1e-11 merit=0.6 converged=0 qualified=0',
        ])
    return '\n'.join(out)

class TestBfsCommittedV64(unittest.TestCase):
    def test_real_rows_only(self):
        a=audit(stream());self.assertEqual(a['momentum_increase_iterations'],[3])
        self.assertFalse(a['solver_qualified_last']);self.assertEqual(a['count'],3)
    def test_missing(self):
        with self.assertRaises(ValueError):audit('')
    def test_truncated(self):
        with self.assertRaises(ValueError):audit(stream().replace('p_continuity=0.1',''))
    def test_missing_sst(self):
        with self.assertRaises(ValueError):audit(stream().replace('phase=SST_TRANSPORT_END','phase=SST_ABSENT'))
    def test_broken_sequence(self):
        with self.assertRaises(ValueError):audit(stream().replace('iteration=2','iteration=4'))
    def test_forged_physics(self):
        with self.assertRaises(ValueError):audit(stream().replace('converged=0 qualified=0','converged=0 qualified=1'))
    def test_nonfinite(self):
        with self.assertRaises(ValueError):audit(stream().replace('U=1 ','U=nan '))
if __name__=='__main__':unittest.main()
