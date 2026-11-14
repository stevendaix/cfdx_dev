"""Fail-closed, synthetic regression for the forensic exact-Schur audit."""
from __future__ import annotations

import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from scipy.io import mmwrite
from scipy.sparse import csr_matrix

SCRIPT = Path(__file__).resolve().parents[2] / "scripts"
sys.path.insert(0, str(SCRIPT))
from diagnose_pcd_schur_v38 import analyze, qualification_gates  # noqa: E402


def save_operator(path, a, dummy_rhs=None):
    matrix = csr_matrix(a)
    n = matrix.shape[0]
    import struct
    vec = np.zeros(n, dtype="<f8")
    if dummy_rhs is not None:
        vec[0] = dummy_rhs
    raw = (
        struct.pack("<8sQQ", b"CFDXFP01", n, matrix.nnz)
        + matrix.indptr.astype("<u4").tobytes()
        + matrix.indices.astype("<u4").tobytes()
        + matrix.data.astype("<f8").tobytes()
        + vec.tobytes() + np.zeros(n, dtype="<f8").tobytes()
    )
    path.write_bytes(raw)


class SchurV38Tests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="cfdx_pcd_v38_")
        self.dir = Path(self.tmp.name)
        n = 2
        u = np.diag(np.arange(1., 7.))
        g = np.zeros((6, n))
        g[0, 0] = 0.5
        g[1, 1] = -0.2
        d = np.zeros((n, 6))
        d[0, 0] = 0.7
        d[1, 1] = -0.9
        C = -np.eye(n) + d @ np.linalg.solve(u, g)
        A = np.block([[u, g], [d, C]])
        self.A = self.dir / "A.mtx"
        self.b = self.dir / "b.mtx"
        self.mp = self.dir / "mp.cfp"
        self.kp = self.dir / "kp.cfp"
        self.fp = self.dir / "fp.cfp"
        mmwrite(self.A, A)
        mmwrite(self.b, np.arange(1., 9.).reshape((-1,1)))
        for path in (self.mp, self.kp, self.fp):
            save_operator(path, np.eye(2))

    def tearDown(self):
        self.tmp.cleanup()

    def run_audit(self):
        return analyze(self.A, self.b, self.mp, self.kp, self.fp)

    def test_exact_schur_and_both_actions(self):
        result = self.run_audit()
        self.assertEqual(result["schema"], "cfdx-pcd-exact-schur-audit/v38")
        self.assertEqual(result["pressure_cells"], 2)
        self.assertLess(result["schur_action_error"]["legacy_frobenius_action_error_over_sqrt_n"], 1e-12)
        self.assertLess(result["schur_action_error"]["commutator_frobenius_action_error_over_sqrt_n"], 1e-12)
        self.assertFalse(result["verified_physical_convergence"])

    def test_captured_rhs_must_be_placeholder_not_physical(self):
        save_operator(self.fp, np.eye(2), dummy_rhs=0.1)
        with self.assertRaisesRegex(ValueError, "ZERO placeholder"):
            self.run_audit()

    def test_wrong_pressure_dimension_is_rejected(self):
        save_operator(self.fp, np.eye(3))
        with self.assertRaisesRegex(ValueError, "dimensions mismatch"):
            self.run_audit()

    def test_wrong_saddle_dimension_is_rejected(self):
        mmwrite(self.A, np.eye(7))
        with self.assertRaisesRegex(ValueError, "exact 3-component"):
            self.run_audit()

    def test_algebraic_gate_can_reject_without_claiming_physics(self):
        result = self.run_audit()
        self.assertEqual(qualification_gates(result, max_action_error=1e-5,
                                             require_fp_nullity_match=True), [])
        result["spectral"]["Fp"]["numerical_nullity"] = 1
        failures = qualification_gates(result, require_fp_nullity_match=True)
        self.assertEqual(len(failures), 1)
        self.assertIn("nullity", failures[0])
        result["schur_action_error"]["legacy_frobenius_action_error_over_sqrt_n"] = 2.0
        self.assertEqual(len(qualification_gates(result, max_action_error=0.1)), 1)
        with self.assertRaises(ValueError):
            qualification_gates(result, max_action_error=float("nan"))

    def test_negative_mass_is_rejected(self):
        save_operator(self.mp, np.diag([-1., 1.]))
        with self.assertRaisesRegex(ValueError, "mass must be positive"):
            self.run_audit()


if __name__ == "__main__":
    unittest.main()
