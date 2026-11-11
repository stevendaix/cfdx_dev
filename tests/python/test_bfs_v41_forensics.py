import importlib.util
import pathlib
import sys
import tempfile
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[2]
spec=importlib.util.spec_from_file_location('bfs_v41_forensics', ROOT/'scripts'/'bfs_v41_forensics.py')
module=importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)


class V41ForensicsTest(unittest.TestCase):
    def test_missing_log_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            with self.assertRaises(ValueError):
                module.parse_production_log(pathlib.Path(temp)/'missing.log')

    def test_only_inner_convergence_not_global(self):
        with tempfile.TemporaryDirectory() as temp:
            p=pathlib.Path(temp)/'valid.log'
            p.write_text('SST_CERTIFIED_INNER iteration=383 k_nonlinear_raw=1e-13 '
                         'omega_nonlinear_raw=2e-10 physical_qualification=NOT_QUALIFIED\n'
                         'GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 merit=0.527 '
                         'converged=0 qualified=0\n')
            d=module.parse_production_log(p)
            self.assertEqual(d['first_sst_inner_certificate']['iteration'],383)
            self.assertEqual(d['committed_outer_states'],1)
            self.assertFalse(d['global_physically_qualified'])

    def test_incomplete_global_certificate_fails_closed(self):
        with tempfile.TemporaryDirectory() as temp:
            p=pathlib.Path(temp)/'bad.log'
            p.write_text('GLOBAL_NONLINEAR_COMMITTED_V13 iteration=1 converged=1\n')
            with self.assertRaises(ValueError):
                module.parse_production_log(p)

    def test_rescue_is_diagnostic_not_certification(self):
        with tempfile.TemporaryDirectory() as temp:
            p=pathlib.Path(temp)/'rescue.log'
            p.write_text('SST_TWO_CYCLE_RESCUE_V41 iteration=62 cell=414 '
                         'next_picard_damping=0.5 physical_qualified=0\n'
                         'SST_FP64_PRECISION_STALL_V14 iteration=249 '
                         'omega_returned_row_scaled=1.00929563449e-8 accepted=0 qualified=0\n')
            d=module.parse_production_log(p)
            self.assertEqual(d['first_two_cycle_rescue']['cell'],414)
            self.assertAlmostEqual(d['first_precision_stall']['omega_returned_row_scaled'],
                                   1.00929563449e-8)
            self.assertFalse(d['global_physically_qualified'])

    def test_no_duplicate_source(self):
        with tempfile.TemporaryDirectory() as temp:
            p=pathlib.Path(temp)/'x.log'; p.write_text('hello\n')
            with self.assertRaises(ValueError):
                module.audit({'a':p,'b':p})


if __name__=='__main__':
    unittest.main()
