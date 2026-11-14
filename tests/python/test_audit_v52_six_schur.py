import copy
import hashlib
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/"scripts"))
from audit_v52_six_schur import audit, file_sha256
from generate_th03_bfs_case import SCHUR_MODELS


class SixSchurV52Tests(unittest.TestCase):
    def setUp(self):
        self.folder=tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        path=Path(self.folder.name)
        self.data={"identical_physical_inputs":True,"cases":{},"results":[]}
        for name in SCHUR_MODELS:
            self.data["cases"][name]={"physical_sha256":"identical-test-input", "n_cells":2}
            captures={}
            for key in ("A","b","x0"):
                fp=path/f"{name}_{key}.mtx"
                fp.write_bytes(f"fixed-{key}".encode("ascii"))
                captures[key]={"path":str(fp),"sha256":file_sha256(fp)}
            log=path/f"{name}.log"
            log.write_text("COUPLED_BLOCK_TRUE_RESIDUAL_V16 full_relative=5e-13\n"+
                 f"BFS_PHASE_V10 iteration=1 phase=COUPLED_LINEAR_END status=0 iterations=10 true_relative=5e-13 actual_krylov=fgmres actual_preconditioner=coupled_block_schur resolved_schur={name}\n")
            self.data["results"].append({"schur":name,
                "operator_capture":captures,"log":str(log),
                "log_sha256":file_sha256(log),"execution":"FAILED",
                "first_blocker":{"first_blocker":"outer_unconverged"},
                "nonlinear_convergence_proven":False,
                "operator_fingerprint_v40":{"A":"a","b":"b","x0":"c",
                    "rows":8,"nnz":100}})

    def test_equal_original_operator_does_not_qualify_six_field_newton(self):
        report=audit(self.data)
        self.assertEqual(report["status"],"COMPARABLE_4N_ONLY")
        self.assertTrue(report["all_six_same_original_4N_A_b_x0_SHA256"])
        self.assertFalse(report["physically_qualified"])
        self.assertFalse(report["globalized_6N_preconditioners_compared"])
        self.assertEqual(len(report["methods"]),6)

    def test_mutated_matrix_is_rejected(self):
        item=self.data["results"][0]["operator_capture"]["A"]
        Path(item["path"]).write_text("CHANGED")
        with self.assertRaises(ValueError): audit(self.data)

    def test_changed_initial_guess_is_incomparable(self):
        other=self.data["results"][0]["operator_capture"]["x0"]
        Path(other["path"]).write_text("different-initial-guess")
        other["sha256"]=file_sha256(other["path"])
        self.assertEqual(audit(self.data)["status"],"INCOMPARABLE_OR_INCOMPLETE")

    def test_duplicate_method_rejected(self):
        self.data["results"][0]["schur"]=self.data["results"][1]["schur"]
        with self.assertRaises(ValueError): audit(self.data)

    def test_missing_log_rejected(self):
        Path(self.data["results"][0]["log"]).unlink()
        with self.assertRaises(ValueError): audit(self.data)

    def test_fingerprint_mismatch_not_comparable(self):
        self.data["results"][0]["operator_fingerprint_v40"]["A"]="not-equal"
        self.assertEqual(audit(self.data)["status"],"INCOMPARABLE_OR_INCOMPLETE")

    def test_six_physical_cases_must_match(self):
        self.data["cases"][SCHUR_MODELS[0]]["physical_sha256"]="different"
        with self.assertRaises(ValueError): audit(self.data)

    def test_silent_schur_substitution_rejected(self):
        obj=self.data["results"][0]
        path=Path(obj["log"])
        path.write_text(path.read_text().replace(
            "resolved_schur="+obj["schur"],"resolved_schur=block_local"))
        obj["log_sha256"]=file_sha256(path)
        with self.assertRaises(ValueError): audit(self.data)

    def test_missing_capture_rejected(self):
        self.data["results"][0].pop("operator_capture")
        with self.assertRaises(ValueError): audit(self.data)

if __name__=="__main__": unittest.main()
