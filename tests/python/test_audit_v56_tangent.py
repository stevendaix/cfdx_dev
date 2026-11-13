import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from audit_v56_tangent import audit

VALID = ("V56_BFS152_FACE_AUDIT: schemes=6 faces=634 verified=3804 nonsmooth=0 "
         "max_primal_error=0 max_tangent_error=3.7e-13 max_cell_parity=7.4e-16 "
         "max_cell_tangent_error=6.9e-10 max_additivity=8.2e-20 "
         "kink_face_id=3 kink_detected=1\n")

class TangentEvidenceGate(unittest.TestCase):
    def test_valid_is_only_local(self):
        a = audit(VALID)
        self.assertTrue(a["local_face_tangent_verified"])
        self.assertFalse(a["global_six_field_jv_qualified"])
        self.assertFalse(a["six_schur_global_newton_qualified"])
        self.assertFalse(a["bfs_physically_converged"])
    def test_missing(self):
        with self.assertRaises(ValueError): audit("")
    def test_double_report(self):
        with self.assertRaises(ValueError): audit(VALID+VALID)
    def test_forged_success(self):
        with self.assertRaises(ValueError): audit(VALID+"six_schur_qualified=true\n")
    def test_nan(self):
        with self.assertRaises(ValueError): audit(VALID.replace("3.7e-13","nan"))
    def test_missing_kink(self):
        with self.assertRaises(ValueError): audit(VALID.replace("kink_detected=1","kink_detected=0"))
    def test_bad_assembly(self):
        with self.assertRaises(ValueError): audit(VALID.replace("6.9e-10","0.1"))
    def test_wrong_mesh(self):
        with self.assertRaises(ValueError): audit(VALID.replace("faces=634","faces=600"))

if __name__ == "__main__":
    unittest.main()
