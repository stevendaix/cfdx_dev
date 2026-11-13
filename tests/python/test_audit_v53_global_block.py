import pathlib
import sys
import tempfile
import unittest
ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / 'scripts'))
from audit_v53_global_block import audit

class TestV53Audit(unittest.TestCase):
    def make_logs(self, folder):
        for method in ("block_local", "lsc", "bfbt", "simple", "simplec", "pcd"):
            text=(f"COUPLED_MATRIX_FINGERPRINT_V40 A=same b=same x0=same rows=608 nnz=7080\n"
                  f"Converged NO Iterations 1\n")
            if method!="pcd":
                text+=(f"FPHYS_JFNK_V52 selected_schur={method} v53_global_block=1 "
                    "native_schur_ready=1 frozen_prestep=1 production_committed=0 qualified=0 "
                    "linear_true_relative=1.0 arnoldi_relative=0.8 jv_linearity=0.2 "
                    "merit_before=2.0 v53_locality_error=0 v53_probes=67 "
                    "v53_transport_nnz=1720 common_state_fingerprint=12345 "
                    "private_accepted=0\n")
            (folder/f"v53_{method}.log").write_text(text)
    def test_six_and_mutations(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder=pathlib.Path(tmp);self.make_logs(folder)
            pcd=folder/"v53_pcd.log"
            pcd.write_text(pcd.read_text().replace("Converged NO Iterations 1",
                "COUPLED_BLOCK_TRUE_RESIDUAL_V16 full_relative=0.006\n"
                "CFDX production solver error: coupled momentum-continuity solve did not converge"))
            self.assertFalse(audit(folder)["physically_qualified"])
            p=folder/"v53_bfbt.log";original=p.read_text()
            for old,new in (("v53_global_block=1","v53_global_block=0"),
                ("private_accepted=0","private_accepted=1"),
                ("common_state_fingerprint=12345","common_state_fingerprint=other"),
                ("linear_true_relative=1.0","linear_true_relative=nan"),
                ("A=same","A=other"),
                ("v53_transport_nnz=1720","v53_transport_nnz=2"),
                ("qualified=0","qualified=1"),
                ("Converged NO","Converged YES")):
                p.write_text(original.replace(old,new))
                with self.assertRaises(ValueError):audit(folder)
            p.write_text(original)
            (folder/"v53_pcd.log").unlink()
            with self.assertRaises(ValueError):audit(folder)
if __name__=='__main__':unittest.main()
