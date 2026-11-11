"""No external CFD fixtures: fail-closed V39 pressure-orientation tests."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.io import mmwrite

SCRIPTS = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPTS))
from diagnose_pcd_orientations_v39 import analyze, main  # noqa: E402
from diagnose_pcd_schur_v38 import normalized_gap  # noqa: E402


def save_cfp(path: Path, array: np.ndarray) -> None:
    from scipy.sparse import csr_matrix
    import struct
    m = csr_matrix(array)
    n = m.shape[0]
    path.write_bytes(struct.pack("<8sQQ", b"CFDXFP01", n, m.nnz)
                     + m.indptr.astype("<u4").tobytes()
                     + m.indices.astype("<u4").tobytes()
                     + m.data.astype("<f8").tobytes()
                     + np.zeros(n*2, dtype="<f8").tobytes())


class V39Orientations(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cfdx_v39_")
        root = Path(self.tmp.name)
        n = 2
        u = np.diag(np.array([2., 4., 3., 5., 6., 7.]))
        g = np.zeros((6, n))
        g[0, 0], g[3, 0], g[1, 1] = 0.5, 0.2, -0.3
        d = g.T
        c = -np.eye(n) + d @ np.linalg.solve(u, g)
        self.A = root / "A.mtx"
        self.b = root / "b.mtx"
        self.mp = root / "mp.cfp"
        self.kp = root / "kp.cfp"
        self.fp = root / "fp.cfp"
        self.out = root / "result.json"
        mmwrite(self.A, np.block([[u,g],[d,c]]))
        mmwrite(self.b, np.arange(1.,9.).reshape((-1,1)))
        save_cfp(self.mp, np.diag([2.,3.]))
        save_cfp(self.kp, np.array([[4.,1.],[-0.5,3.]]))
        save_cfp(self.fp, np.array([[2.,1.5],[-0.7,1.]]))

    def tearDown(self):
        self.tmp.cleanup()

    def audit(self):
        return analyze(self.A, self.b, self.mp, self.kp, self.fp)

    def args(self):
        return ["--coupled-A",str(self.A),"--coupled-b",str(self.b),
                "--mp",str(self.mp),"--kp",str(self.kp),"--fp",str(self.fp),
                "--output",str(self.out)]

    def test_noncommuting_orders_produce_distinct_defects(self):
        report=self.audit()
        self.assertEqual(report["schema"],"cfdx-pcd-orientation-audit/v39")
        actions=report["actions"]
        self.assertNotAlmostEqual(actions["v39_negative_Kinv_Fp_Minv"]["frobenius_defect_over_sqrt_n"],
                                  actions["old_negative_Minv_Fp_Kinv"]["frobenius_defect_over_sqrt_n"])
        self.assertFalse(report["physical_converged"])
        self.assertFalse(report["selected_for_production"])
        self.assertFalse(report["full_nonlinear_U_p_k_omega_jacobian"])

    def test_oracle_exact_on_diagonal_velocity_test(self):
        report = self.audit()
        defect = report["actions"]["scalar_velocity_diagonal_exact_dense_inverse_oracle_only"]
        self.assertLess(defect["frobenius_defect_over_sqrt_n"],1e-12)
        block=report["actions"]["cell_local_3x3_exact_dense_inverse_oracle_only"]
        self.assertLess(block["frobenius_defect_over_sqrt_n"],1e-12)

    def test_cli_report_is_atomic_and_contains_gate(self):
        self.assertEqual(main(self.args()),0)
        report=json.loads(self.out.read_text())
        self.assertTrue(report["algebraic_gate"]["pass"])
        self.assertFalse(report["algebraic_gate"]["requested"])

    def test_strict_gate_fails_without_qualifying_physics(self):
        self.assertEqual(main(self.args()+["--maximum-v39-action-defect","1e-9"]),2)
        report=json.loads(self.out.read_text())
        self.assertFalse(report["algebraic_gate"]["pass"])
        self.assertFalse(report["physical_converged"])

    def test_optional_native_mesh_hotspots(self):
        try:
            import h5py
        except ImportError:
            self.skipTest("h5py not installed; optional geometry enrichment")
        mesh = self.out.with_suffix(".h5")
        with h5py.File(mesh, "w") as h5:
            h5.create_dataset("points", data=np.array([
                [0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
                [2.,0.,0.],[3.,0.,0.],[3.,1.,0.],[2.,1.,0.]]))
            h5.create_dataset("cell_offsets", data=np.array([0,1,2]))
            h5.create_dataset("cell_faces", data=np.array([0,1]))
            h5.create_dataset("face_offsets", data=np.array([0,4,8]))
            h5.create_dataset("face_vertices", data=np.arange(8))
        j = analyze(self.A, self.b, self.mp, self.kp, self.fp,
                    mesh_file=mesh)
        self.assertTrue(j["mesh_geometry_used"])
        for row in j["largest_S_vs_Fp_row_norm_ratios"]:
            self.assertEqual(len(row["center_m"]),3)
            self.assertAlmostEqual(row["center_m"][0],0.5+2.0*row["cell"])
        with h5py.File(mesh,"r+") as h5:
            del h5["cell_faces"]
        with self.assertRaisesRegex(ValueError,"lacks native CSR"):
            analyze(self.A,self.b,self.mp,self.kp,self.fp,mesh_file=mesh)

    def test_nonlumped_mass_rejected(self):
        save_cfp(self.mp,np.array([[2.,0.1],[0.1,3.]]))
        with self.assertRaisesRegex(ValueError,"exact diagonal pressure mass"):
            self.audit()

    def test_wrong_coupled_dimension_rejected(self):
        mmwrite(self.A,np.eye(7))
        with self.assertRaisesRegex(ValueError,"exact 3-component"):
            self.audit()

    def test_zero_fp_placeholder_fail_closed(self):
        # Captured matrix-only RHS must remain zero; test V38 protects this.
        raw=bytearray(self.fp.read_bytes())
        raw[-16:-8]=np.float64(4.0).tobytes()
        self.fp.write_bytes(raw)
        with self.assertRaisesRegex(ValueError,"ZERO placeholder"):
            self.audit()

    def test_nonzero_probe(self):
        r=np.array([1.,2.])
        s=np.array([[2.,0.],[0.,3.]])
        out=normalized_gap(r,np.linalg.solve(s,r),s)
        self.assertLess(out["schur_true_residual_relative"],1e-13)


if __name__=="__main__":
    unittest.main()
