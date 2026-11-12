from pathlib import Path
import importlib.util
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location(
    "v42_recovery", ROOT / "scripts" / "audit_bfs_pcd_recovery_v42.py")
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
FINGERPRINT = ("COUPLED_MATRIX_FINGERPRINT_V40 A=aa670a634d0fee97 "
               "b=cbcfa394bfde5d1d x0=bca7010f84b2cdef rows=608 nnz=7080\n")
BAD = "COUPLED_BLOCK_TRUE_RESIDUAL_V16 full_relative=6.27493e-3 krylov_status=1\n"
RECOVER = ("COUPLED_PCD_BLOCKLOCAL_RECOVERY_V42 attempted=1 adopted=1 "
           "primary_model=PCD recovery_model=BlockLocal primary_status=1 "
           "primary_iterations=300 primary_budget=300 recovery_status=0 "
           "recovery_iterations=40 true_relative=5.201229499689e-13 "
           "unchanged_tolerance=1e-12 original_A_unchanged=1 physical_qualified=0\n")


class TestV42Recovery(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.b = Path(self.temp.name) / "pcd.log"
        self.r = Path(self.temp.name) / "retry.log"
        self.b.write_text(FINGERPRINT + BAD)
        self.r.write_text(FINGERPRINT + RECOVER)

    def test_verified_linear_only(self):
        r = module.audit(self.b, self.r)
        self.assertEqual(r["qualification"], "LINEAR_RECOVERY_ONLY")
        self.assertTrue(r["frozen_linear_system_converged"])
        self.assertFalse(r["pure_pcd_convergence_proven"])
        self.assertFalse(r["global_physical_convergence_proven"])

    def test_not_the_same_frozen_system(self):
        self.r.write_text(self.r.read_text().replace("aa670a634d0fee97", "bb670a634d0fee97"))
        with self.assertRaisesRegex(ValueError, "different frozen"):
            module.audit(self.b, self.r)

    def test_unqualified_retry_rejected(self):
        self.r.write_text(self.r.read_text().replace("adopted=1", "adopted=0"))
        with self.assertRaisesRegex(ValueError, "honest"):
            module.audit(self.b, self.r)

    def test_false_convergence_rejected(self):
        self.r.write_text(self.r.read_text().replace("5.201229499689e-13", "2e-12"))
        with self.assertRaisesRegex(ValueError, "true-residual"):
            module.audit(self.b, self.r)

    def test_missing_evidence_rejected(self):
        self.r.write_text(FINGERPRINT)
        with self.assertRaisesRegex(ValueError, "never ran"):
            module.audit(self.b, self.r)


if __name__ == "__main__":
    unittest.main()
