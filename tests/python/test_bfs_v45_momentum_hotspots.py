"""Mutation regressions for fail-closed row/COMMITTED association."""
from __future__ import annotations
import sys
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'scripts'))
from analyze_bfs_momentum_hotspots_v45 import analyze


def make_row(component, relative):
    value = relative * 2
    return (f'POST_SST_MOMENTUM_ROW_V45 iteration=1 component={component} cell=12 '
            f'original_signed={value} original_abs={value} rhs_scale=2 rhs_relative={relative} '
            f'diagonal_term={value+1} off_diagonal_terms=0 rhs=1 '
            f'compensated_signed={value} pressure_source_integral=0.25 '
            f'body_source_integral=0.25 pressure_body_rhs_remainder=0.5 '
            f'nu_t=1 k=1 omega=1 committed=0 qualified=0\n')

VALID = (make_row(0,0.03)+make_row(1,0.04)+make_row(2,0.9)+
         'GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 U=0.04 p_continuity=0.001 k=1e-12 '
         'omega=1e-13 nut_closure=0 merit=0.04 converged=0 qualified=0\n'
         'Converged NO Iterations 1\n')

class TestHotspotAudit(unittest.TestCase):
    def test_valid_extruded_2d_ignores_inactive_z(self):
        result = analyze(VALID,sample=(1,))
        self.assertEqual(result['samples']['1']['dominant_component'], 1)
        self.assertEqual(result['committed_rows'],1)
        self.assertFalse(result['physical_converged'])
    def test_catches_relative_mismatch(self):
        with self.assertRaisesRegex(ValueError,'relative scale mismatch'):
            analyze(VALID.replace('rhs_relative=0.04','rhs_relative=0.041'))
    def test_catches_inconsistent_raw_csr(self):
        with self.assertRaisesRegex(ValueError,'independent CSR residual'):
            analyze(VALID.replace('compensated_signed=0.08','compensated_signed=0.081'))
    def test_catches_misreported_gate(self):
        with self.assertRaisesRegex(ValueError,'COMMITTED momentum gate'):
            analyze(VALID.replace('U=0.04','U=0.05'))
    def test_catches_missing_component(self):
        with self.assertRaisesRegex(ValueError,'lacks reassembled'):
            analyze(VALID.replace(make_row(1,0.04),''))
    def test_catches_fake_certification(self):
        with self.assertRaisesRegex(ValueError,'incorrectly claims'):
            analyze(VALID.replace('committed=0 qualified=0','committed=1 qualified=0'))
    def test_catches_fake_global_converged(self):
        with self.assertRaisesRegex(ValueError,'final physical status'):
            analyze(VALID.replace('Converged NO','Converged YES'))
    def test_catches_truncation(self):
        with self.assertRaisesRegex(ValueError,'final physical status'):
            analyze(VALID.replace('Converged NO Iterations 1',''))
    def test_catches_forced_xyz_gate(self):
        with self.assertRaisesRegex(ValueError,'COMMITTED momentum gate'):
            analyze(VALID,active_components=(0,1,2))

if __name__=='__main__':
    unittest.main()
