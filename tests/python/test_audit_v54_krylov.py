"""Fail-closed mutation coverage for V54 physical/Krylov fidelity evidence."""
from __future__ import annotations
import sys
from pathlib import Path
import unittest
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"scripts"))
from audit_v54_krylov import audit

SAMPLE = ("FPHYS_JFNK_V52 iteration=1 step=0 v54_frozen_amat=1 "
          "v54_jacobian_evaluations=24 v54_basis_audits=1 "
          "v54_centered_columns=12 v54_one_sided_columns=0 "
          "v54_fidelity=7.9e-2 v54_halving=6e-5 "
          "v54_sided_nonsmooth=4.4e-1 v54_fidelity_ok=0 "
          "linear_true_relative=inf arnoldi_relative=inf "
          "private_accepted=0 production_committed=0 qualified=0 "
          "common_state_fingerprint=18386539465314058853 "
          "reason=V54_frozen_Amat_fails_independent_Fphys_Krylov_basis_fidelity\n")
class TestAuditV54(unittest.TestCase):
    def test_negative_proof_is_valid_but_not_qualified(self):
        data=audit(SAMPLE,expected_cells=2)
        self.assertEqual(data["status"],"BLOCKED_JV_FIDELITY")
        self.assertFalse(data["qualified"])
        self.assertFalse(data["six_schur_6N_requalified"])
    def test_mutations_must_fail(self):
        changes=(
            ("v54_frozen_amat=1","v54_frozen_amat=0"),
            ("v54_basis_audits=1","v54_basis_audits=0"),
            ("v54_jacobian_evaluations=24","v54_jacobian_evaluations=0"),
            ("v54_centered_columns=12","v54_centered_columns=11"),
            ("private_accepted=0","private_accepted=1"),
            ("production_committed=0","production_committed=1"),
            ("qualified=0","qualified=1"),
            ("v54_fidelity_ok=0","v54_fidelity_ok=1"),
            ("v54_halving=6e-5","v54_halving=nan"),
            ("v54_sided_nonsmooth=4.4e-1","v54_sided_nonsmooth=-1"),
            ("common_state_fingerprint=18386539465314058853","common_state_fingerprint=0"),
            ("reason=V54_frozen_Amat_fails_independent_Fphys_Krylov_basis_fidelity",
             "reason=Arnoldi_passed"),
        )
        for before,after in changes:
            with self.subTest(change=before),self.assertRaises(ValueError):
                audit(SAMPLE.replace(before,after),expected_cells=2)
        with self.assertRaises(ValueError):
            audit(SAMPLE+SAMPLE,expected_cells=2)
    def test_consistent_fidelity_proof_is_not_a_convergence_certificate(self):
        stable=(SAMPLE.replace("v54_fidelity=7.9e-2","v54_fidelity=1e-8")
               .replace("v54_sided_nonsmooth=4.4e-1","v54_sided_nonsmooth=1e-7")
               .replace("v54_fidelity_ok=0","v54_fidelity_ok=1"))
        data=audit(stable,expected_cells=2)
        self.assertEqual(data["status"],"JACOBIAN_FIDELITY_ONLY")
        self.assertFalse(data["qualified"])
        self.assertFalse(data["global_newton_committed"])
if __name__=="__main__":unittest.main()
