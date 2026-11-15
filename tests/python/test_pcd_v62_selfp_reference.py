"""V62 fail-closed falsifications of six-input selfp/ILU reference oracle."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csc_matrix, eye, diags
from scripts.pcd_v62_selfp_fgmres_reference import METHODS, analyze

class SelfpV62Tests(unittest.TestCase):
    def setUp(self):
        self.dir=tempfile.TemporaryDirectory(prefix='cfdx-v62-selfp-')
        self.p=Path(self.dir.name)
        n=4
        U=diags([3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.]).tocsc()
        G=csc_matrix(np.tile(np.eye(n),(3,1))*0.4)
        D=G.T.tocsc()
        C=-eye(n,format='csc')
        from scipy.sparse import bmat
        A=bmat([[U,G],[D,C]],format='csc')
        rhs=np.array([np.linspace(0.1,2.,A.shape[0])]).T
        zeros=np.zeros_like(rhs)
        for model in METHODS:
            mmwrite(self.p/f'{model}_A.mtx',A)
            mmwrite(self.p/f'{model}_b.mtx',rhs)
            mmwrite(self.p/f'{model}_x0.mtx',zeros)
    def tearDown(self):self.dir.cleanup()
    def test_original_six_reference(self):
        out=analyze(self.p,max_cells=16)
        self.assertEqual(out['selfp_rank'],4)
        self.assertTrue(out['six_original_captures_identical'])
        self.assertTrue(out['approximate_ilu_blocks']['meets_reference_target'])
        self.assertLess(out['approximate_ilu_blocks']['true_relative_residual'],1e-9)
        self.assertFalse(out['PCD_production_converged'])
        self.assertFalse(out['BFS_physical_qualified'])
    def test_mismatch_rejected(self):
        mmwrite(self.p/'lsc_b.mtx',np.ones((16,1)))
        with self.assertRaisesRegex(ValueError,'not identical'):
            analyze(self.p)
    def test_missing_original_input_rejected(self):
        (self.p/'simplec_A.mtx').unlink()
        with self.assertRaisesRegex(ValueError,'missing'):
            analyze(self.p)
    def test_changed_pressure_gauge_rejected(self):
        # Destroy a pressure row/column; the approximate Schur is singular.
        for model in METHODS:
            matrix=csc_matrix(np.asarray(__import__('scipy').io.mmread(self.p/f'{model}_A.mtx').toarray()))
            matrix=matrix.tolil();matrix[12,:]=0.;matrix[:,12]=0.
            mmwrite(self.p/f'{model}_A.mtx',matrix.tocsc())
        with self.assertRaisesRegex(ValueError,'singular'):
            analyze(self.p)
    def test_unbounded_cases_rejected(self):
        for kw in (dict(max_cells=257),dict(ilu_drop_tol=.5),
                   dict(ilu_fill_factor=100),dict(relative_target=.1)):
            with self.assertRaises(ValueError):analyze(self.p,**kw)

if __name__=='__main__':unittest.main()