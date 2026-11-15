import tempfile
from pathlib import Path
import unittest
from scripts.plot_bfs_sst_v81 import parse_line, get_number, read_native

END='BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_END committed=0 converged=0 iterations=40 k_residual=1 omega_residual=2\n'
START='BFS_PHASE_V10 iteration=1 phase=SST_TRANSPORT_BEGIN budget=40 target=1e-8\n'
DIAG=('SST_TRANSPORT_DIAGNOSTIC iteration=40 '
      'k_returned_row_scaled=0.02 omega_returned_row_scaled=100 '
      'k_local_relative=0.03 omega_local_relative=0.005 '
      'dk_inf=0.01 domega_inf=1000\n')
class TestV81(unittest.TestCase):
    def check(self, text, should_pass):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'native.log';p.write_text(text)
            if should_pass:
                x=read_native(p)
                self.assertEqual(x['end']['committed'],'0')
            else:
                with self.assertRaises(ValueError):read_native(p)
    def test_reference(self):self.check(START+DIAG+END,True)
    def test_no_end(self):self.check(START+DIAG,False)
    def test_false_commit(self):self.check(START+DIAG+END.replace('committed=0','committed=1'),False)
    def test_missing_physical(self):self.check(START+END,False)
    def test_loosened_target(self):self.check(START.replace('1e-8','1e-5')+DIAG+END,False)
    def test_nonfinite(self):self.check(START+DIAG.replace('omega_returned_row_scaled=100','omega_returned_row_scaled=nan')+END,False)
    def test_duplicate_metric(self):self.check(START+DIAG.replace('dk_inf=','dk_inf=0.4 dk_inf=')+END,False)
    def test_budget_falsification(self):
        budget=('SST_CELL_BUDGET iteration=40 cell=1 field=omega convection=0 diffusion=-2 '
                'physical_production=-3 physical_destruction=2 physical_cross=0 '
                'residual=4 decomposition_error=0\n')
        self.check(START+DIAG+budget+END,False)
    def test_invalid_iter(self):self.check(START+DIAG.replace('iteration=40','iteration=41')+END,False)
    def test_negative(self):self.check(START+DIAG.replace('dk_inf=0.01','dk_inf=-0.01')+END,False)
    def test_duplicate_field(self):
        with self.assertRaises(ValueError):parse_line('A k=2 k=3')
    def test_missing_number(self):
        with self.assertRaises(ValueError):get_number({},'omega')
if __name__=='__main__':unittest.main()