@@ -0,0 +1,62 @@
"""Fail-closed perturbation tests for the V61 actual-capture Schur oracle."""
import sys
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csc_matrix
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from pcd_v61_exact_schur_oracle import analyze,METHODS


class OracleTest(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.p=Path(self.tmp.name)
        # Genuine saddle point with cross-field gradients, no synthetic RHS
        # success certificate can be extrapolated to a real turbulent BFS.
        self.A=np.eye(8)*2.
        self.A[0,6]=-.7;self.A[2,7]=.1
        self.A[6,0]=.2;self.A[7,2]=-.4
        self.A[6,6]=.05;self.A[7,7]=.08
        self.b=np.arange(1,9,dtype=float).reshape(-1,1)
        self.x=np.zeros((8,1))
        self.write()
    def write(self):
        for m in METHODS:
            mmwrite(self.p/f'{m}_A.mtx',csc_matrix(self.A))
            mmwrite(self.p/f'{m}_b.mtx',self.b)
            mmwrite(self.p/f'{m}_x0.mtx',self.x)
    def test_solvable_same_input_without_pcd_qualification(self):
        result=analyze(self.p)
        self.assertTrue(result['same_A_b_x0_sha256'])
        self.assertLess(result['direct_LU_original_true_relative'],1e-12)
        self.assertEqual(result['cells'],2)
        self.assertFalse(result['pcd_converged'])
        self.assertFalse(result['physical_BFS_converged'])
    def test_modified_one_preconditioner_capture_rejected(self):
        f=self.p/'simplec_A.mtx'
        changed=self.A.copy();changed[0,0]+=0.1
        mmwrite(f,csc_matrix(changed))
        with self.assertRaises(ValueError):analyze(self.p)
    def test_missing_one_capture_rejected(self):
        (self.p/'pcd_b.mtx').unlink()
        with self.assertRaises(ValueError):analyze(self.p)
    def test_singular_coupled_matrix_rejected(self):
        self.A[-1]=0
        self.write()
        with self.assertRaises(ValueError):analyze(self.p)
    def test_nonfinite_original_rejected(self):
        self.b[0,0]=float('nan')
        self.write()
        with self.assertRaises(ValueError):analyze(self.p)
    def test_wrong_warm_start_rejected(self):
        self.x=np.zeros((7,1))
        self.write()
        with self.assertRaises(ValueError):analyze(self.p)
    def test_large_dense_schur_must_be_opt_in(self):
        with self.assertRaises(ValueError):analyze(self.p,maximum_cells=1)

if __name__=='__main__':unittest.main()
