"""Mutation tests for same-state 6N JFNK Schur audit; no CFD solve."""
from __future__ import annotations
import copy
from pathlib import Path
import sys
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
from audit_v52_frozen_jfnk import NAMES,audit_logs

AUDIT={"status":"COMPARABLE_4N_ONLY","all_six_same_original_4N_A_b_x0_SHA256":True,"n_cells":152}
def fake(name):
    return ("COUPLED_MATRIX_FINGERPRINT_V40 A=123 b=456 x0=789 rows=608\n"
            f"FPHYS_JFNK_V52 iteration=1 step=0 selected_schur={name} "
            "frozen_prestep=1 common_state_fingerprint=123456 "
            "native_schur_ready=1 production_committed=0 qualified=0 "
            "linear_true_relative=0.47 merit_before=218.0 evaluations=33 private_accepted=0\n")

class TestFrozenJFNK(unittest.TestCase):
    def setUp(self):self.logs={name:fake(name) for name in NAMES}
    def reject(self,logs=None,audit=None):
        with self.assertRaises(ValueError):audit_logs(logs or self.logs,audit or AUDIT)
    def test_baseline(self):
        r=audit_logs(self.logs,AUDIT)
        self.assertEqual(r["status"],"SAME_FROZEN_X_5_METHODS_NO_NEWTON_CONVERGENCE")
        self.assertFalse(r["physically_qualified"])
        self.assertEqual(len(r["methods"]),6)
    def test_mismatched_initial_state(self):
        logs=copy.deepcopy(self.logs);logs["lsc"]=logs["lsc"].replace("123456","777777")
        self.reject(logs)
    def test_mismatched_four_n_a(self):
        logs=copy.deepcopy(self.logs);logs["bfbt"]=logs["bfbt"].replace("A=123","A=999")
        self.reject(logs)
    def test_method_substituted(self):
        logs=copy.deepcopy(self.logs);logs["simple"]=logs["simple"].replace("selected_schur=simple","selected_schur=pcd")
        self.reject(logs)
    def test_missing_schur(self):
        logs=copy.deepcopy(self.logs);logs.pop("lsc");self.reject(logs)
    def test_invalid_four_n(self):
        prior=copy.deepcopy(AUDIT);prior["status"]="INCOMPARABLE_OR_INCOMPLETE"
        self.reject(audit=prior)
    def test_fake_qualification(self):
        logs=copy.deepcopy(self.logs);logs["simplec"]=logs["simplec"].replace("qualified=0","qualified=1")
        self.reject(logs)
    def test_fake_newton_accepted(self):
        logs=copy.deepcopy(self.logs);logs["simplec"]=logs["simplec"].replace("private_accepted=0","private_accepted=1")
        self.reject(logs)
    def test_nonfinite_residual(self):
        logs=copy.deepcopy(self.logs);logs["simplec"]=logs["simplec"].replace("relative=0.47","relative=nan")
        self.reject(logs)
    def test_duplicate_evaluation(self):
        logs=copy.deepcopy(self.logs);logs["lsc"]+=logs["lsc"]
        self.reject(logs)

if __name__=="__main__":unittest.main()
