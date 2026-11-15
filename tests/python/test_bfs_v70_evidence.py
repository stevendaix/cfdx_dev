@@ -0,0 +1,83 @@
"""Falsification gates for V70 native replay, no fake CFD certification."""
import unittest
from scripts.audit_bfs_v70_conservative import check_v70_evidence


def baseline():
    rows = []
    rows += [f"GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration={i} committed=0"
             for i in range(149)]
    rows += [f"GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration={i}"
             for i in range(22)]
    rows += [f"BFS_SST_RC_PICARD_V69 cycle={i} committed=0"
             for i in range(6)]
    rows += [f"INCOMPRESSIBLE_ITER iter={i}" for i in range(164)]
    return "\n".join(rows) + "\n"


def sample():
    return (baseline() +
        "BFS_SST_RC_CONSERVATIVE_V70 iteration=23 retry=8 cycle=1 fraction=1 "
        "mass_before=0.01 mass_frozen_after=1e-13 mass_true_after_sst=0.005 "
        "pressure_delta=1 face_delta=0.1 merit_ratio=0.9 worst_growth=1.0 "
        "candidate_eligible=0 private_trial=1 committed=0 qualified=0\n" +
        "BFS_SST_RC_TRANSACTION_V70 iteration=23 retry=8 "
        "U_p_phi_k_omega_nut_unchanged=1 full_outer_rollback_available=1 "
        "committed=0 qualified=0\n")


class TestV70Audit(unittest.TestCase):
    def test_clean(self):
        result = check_v70_evidence(sample(), baseline())
        self.assertEqual(result["v70_projection_cycles"], 1)
        self.assertFalse(result["real_v70_commit"])

    def test_historical_change(self):
        with self.assertRaises(ValueError):
            check_v70_evidence(sample().replace("INCOMPRESSIBLE_ITER iter=3",
                                                  "INCOMPRESSIBLE_ITER iter=5", 1), baseline())

    def test_hidden_publication(self):
        with self.assertRaises(ValueError):
            check_v70_evidence(sample().replace("private_trial=1 committed=0",
                                                  "private_trial=1 committed=1"), baseline())

    def test_missing_rollback(self):
        with self.assertRaises(ValueError):
            check_v70_evidence(sample().replace("U_p_phi_k_omega_nut_unchanged=1",
                                                  "U_p_phi_k_omega_nut_unchanged=0"), baseline())

    def test_frozen_divergence_failure(self):
        with self.assertRaises(ValueError):
            check_v70_evidence(sample().replace("mass_frozen_after=1e-13",
                                                  "mass_frozen_after=1e-4"), baseline())

    def test_partial_improves_without_false_full_certificate(self):
        trial = sample().replace("fraction=1 ", "fraction=0.5 ").replace(
            "mass_frozen_after=1e-13", "mass_frozen_after=0.005")
        result = check_v70_evidence(trial, baseline())
        self.assertEqual(result["v70_projection_cycles"], 1)
        self.assertFalse(result["real_v70_commit"])

    def test_partial_cannot_increase_mass(self):
        trial = sample().replace("fraction=1 ", "fraction=0.5 ").replace(
            "mass_frozen_after=1e-13", "mass_frozen_after=0.02")
        with self.assertRaises(ValueError):
            check_v70_evidence(trial, baseline())

    def test_false_eligibility(self):
        with self.assertRaises(ValueError):
            check_v70_evidence(sample().replace("candidate_eligible=0", "candidate_eligible=1"), baseline())

    def test_nonnumeric(self):
        with self.assertRaises(ValueError):
            check_v70_evidence(sample().replace("pressure_delta=1",
                                                  "pressure_delta=nan"), baseline())

    def test_missing_probe(self):
        with self.assertRaises(ValueError):
            check_v70_evidence(baseline(), baseline())


if __name__ == "__main__":
    unittest.main()
