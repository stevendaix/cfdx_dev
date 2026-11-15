import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from scripts.bfs_pressure_precision_oracle_v79 import run, factor_and_split

class PressureOracleV79Tests(unittest.TestCase):
    def test_balanced_matrix_and_nonpromotion(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp)
            A=csr_matrix(np.array([[2.,-1.],[-1.,2.]]))
            mmwrite(d/'m.mtx', A);np.savetxt(d/'b.rhs', np.array([1.,0.]),fmt='%.17g')
            rep=run(d/'m.mtx',d/'b.rhs',d/'proof',1e-11,6)
            self.assertTrue(rep['algebraically_closed'])
            self.assertFalse(rep['physical_face_flux_certified'])
            self.assertFalse(rep['production_commit_authorized'])
            self.assertIn('CFDX_PRESSURE_SPLIT_V79', (d/'proof'/'oracle_v79.split').read_text())
    def test_fail_closed_wrong_rhs(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp); A=csr_matrix(np.eye(2));mmwrite(d/'m.mtx',A)
            np.savetxt(d/'b.rhs',np.array([1.,2.,3.]))
            with self.assertRaises(ValueError):factor_and_split(d/'m.mtx',d/'b.rhs')
    def test_fail_closed_singular(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp);A=csr_matrix(np.array([[1.,-1.],[-1.,1.]]))
            mmwrite(d/'m.mtx',A); np.savetxt(d/'b.rhs',np.array([1.,-1.]))
            with self.assertRaises(RuntimeError):factor_and_split(d/'m.mtx',d/'b.rhs')
    def test_fail_closed_threshold(self):
        with tempfile.TemporaryDirectory() as temp:
            d=Path(temp);mmwrite(d/'m.mtx',csr_matrix(np.eye(2)))
            np.savetxt(d/'b.rhs',np.array([1.,2.]))
            with self.assertRaises(ValueError):factor_and_split(d/'m.mtx',d/'b.rhs',0.0)

if __name__=='__main__': unittest.main()