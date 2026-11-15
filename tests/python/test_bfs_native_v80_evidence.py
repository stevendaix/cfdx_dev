@@ -0,0 +1,32 @@
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from plot_bfs_native_v80 import parse

CSR='PRESSURE_NATIVE_CSR_V80 iteration=1 corrector=0 initial_relative=2e-10 final_relative=1e-18 production_commit_authorized=0'
FACE='PRESSURE_NATIVE_FACE_V80 iteration=1 corrector=0 initial_relative=3e-10 final_relative=3e-14 true_face_converged=1 requested_tolerance=1e-11 production_commit_authorized=0'
FLUX='PRESSURE_HILO_PHYSICAL_FLUX iteration=1 corrector=0 exact_relative=5e-14 stored_fp64_relative=7e-14 accepted_iteration=0'
FULL='BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_BEGIN\nBFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_END committed=1 converged=1\nGLOBAL_NONLINEAR_V8 iteration=1 U=1 p_continuity=0.01 k=0.1 omega=0.1 merit=1 accepted=0 qualified=0\nGLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=1 p_continuity=0.01 k=0.1 omega=0.1 converged=0 qualified=0'

class NativeEvidenceTests(unittest.TestCase):
    def test_pressure_only(self):
        x=parse('\n'.join([CSR,FACE,FLUX]));self.assertTrue(x['pressure_native_verified']);self.assertEqual(x['accepted_outer_iterations'],[]);self.assertFalse(x['converged_cfd'])
    def test_outer_transaction_requires_sst(self):
        x=parse('\n'.join([CSR,FACE,FLUX,FULL]));self.assertEqual(x['accepted_outer_iterations'],[1]);self.assertFalse(x['converged_cfd'])
    def test_missing_face(self):
        with self.assertRaises(ValueError):parse('\n'.join([CSR,FLUX]))
    def test_forged_commit(self):
        with self.assertRaises(ValueError):parse('\n'.join([CSR.replace('authorized=0','authorized=1'),FACE,FLUX]))
    def test_missing_sst_certificate(self):
        with self.assertRaises(ValueError):parse('\n'.join([CSR,FACE,FLUX,FULL.replace('phase=SST_TRANSPORT_END committed=1 converged=1','phase=SST_TRANSPORT_END committed=0 converged=0')]))
    def test_missing_global_evaluation(self):
        with self.assertRaises(ValueError):parse('\n'.join([CSR,FACE,FLUX,FULL.replace('GLOBAL_NONLINEAR_V8 iteration=1 U=1 p_continuity=0.01 k=0.1 omega=0.1 merit=1 accepted=0 qualified=0\n','')]))
    def test_weak_gate(self):
        with self.assertRaises(ValueError):parse('\n'.join([CSR,FACE.replace('tolerance=1e-11','tolerance=1e-9'),FLUX]))
    def test_bad_fp64_flux(self):
        self.assertFalse(parse('\n'.join([CSR,FACE,FLUX.replace('stored_fp64_relative=7e-14','stored_fp64_relative=1e-8')]))['pressure_native_verified'])
    def test_nonfinite(self):
        with self.assertRaises(ValueError):parse('\n'.join([CSR,FACE.replace('final_relative=3e-14','final_relative=nan'),FLUX]))

if __name__=='__main__': unittest.main()
