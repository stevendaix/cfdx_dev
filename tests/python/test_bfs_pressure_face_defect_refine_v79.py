import tempfile
import unittest
from pathlib import Path
import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix
from scripts.bfs_pressure_face_defect_refine_v79 import refine, read_defect

class FaceDefectV79Tests(unittest.TestCase):
    def test_one_newton_step_on_an_inexact_face_rhs(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t)
            a=csr_matrix(np.array([[2.,-1.],[-1.,2.]]))
            mmwrite(d/'a.mtx',a);np.savetxt(d/'b.rhs',np.array([1.,0.]),fmt='%.17g')
            (d/'start.split').write_text('CFDX_PRESSURE_SPLIT_V79\n2 4\n0 0\n0 0\n')
            (d/'defect.vec').write_text('CFDX_PRESSURE_FACE_DEFECT_V79\n2 4\n1\n0\n')
            rep=refine(d/'a.mtx',d/'b.rhs',d/'start.split',d/'defect.vec',d/'out')
            self.assertLess(rep['csr_after'],1e-14)
            self.assertFalse(rep['production_commit_authorized'])
            self.assertEqual(rep['physical_face_after'],'REQUIRES_NATIVE_REEVALUATION')
    def test_reject_truncated_wrong_geometry(self):
        with tempfile.TemporaryDirectory() as t:
            d=Path(t);(d/'defect.vec').write_text('CFDX_PRESSURE_FACE_DEFECT_V79\n3 4\n1\n0\n')
            with self.assertRaises(ValueError):read_defect(d/'defect.vec',2,4)
            (d/'defect.vec').write_text('CFDX_PRESSURE_FACE_DEFECT_V79\n2 4\n1\n')
            with self.assertRaises(ValueError):read_defect(d/'defect.vec',2,4)

if __name__=='__main__':unittest.main()