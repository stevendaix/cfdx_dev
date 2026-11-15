"""Fail-closed V84 native progress evidence/plot tests."""
from pathlib import Path
import importlib.util
import os
import tempfile
import unittest
import numpy as np

SCRIPT=Path(__file__).resolve().parents[2]/'scripts/plot_bfs_v84_progress.py'
SPEC=importlib.util.spec_from_file_location('cfdx_v84_progress',SCRIPT)
M=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

COMMIT=('GLOBAL_PROGRESS_COMMITTED_V84 iteration={iteration} committed=2 '
        'k_debt=1 omega_debt=0 merit={merit} U={U} U_global=0.4464 '
        'k=0.862 omega=0.008 p_continuity=1e-16 converged=0 qualified=0\n')

class V84Evidence(unittest.TestCase):
    def write(self,text):
        self.tmp=tempfile.TemporaryDirectory()
        p=Path(self.tmp.name)/'native.stderr';p.write_text(text)
        self.addCleanup(self.tmp.cleanup)
        return p

    def test_accepted_native_only(self):
        p=self.write(COMMIT.format(iteration=2,merit='0.995',U='0.995'))
        trials,records=M.scan_ledger(p)
        self.assertEqual(len(records),1)
        self.assertEqual(records[0]['iteration'],2)

    def test_rejected_trial_is_not_accepted(self):
        p=self.write('GLOBAL_PROGRESS_TRIAL_V84 iteration=2 accepted=0 reason=k_growth_budget\n')
        trials,records=M.scan_ledger(p)
        self.assertEqual(len(records),0)
        self.assertEqual(len(trials),1)

    def test_qualified_forbidden_in_progress_record(self):
        p=self.write(COMMIT.format(iteration=2,merit='0.995',U='0.995').replace('qualified=0','qualified=1'))
        with self.assertRaises(ValueError):M.scan_ledger(p)

    def test_negative_nonnfinite_or_unsorted_rejected(self):
        for merit in ('-0.5','nan'):
            with self.subTest(merit=merit):
                p=self.write(COMMIT.format(iteration=2,merit=merit,U='0.995'))
                with self.assertRaises(ValueError):M.scan_ledger(p)
        p=self.write(COMMIT.format(iteration=2,merit='0.995',U='0.995')*2)
        with self.assertRaises(ValueError):M.scan_ledger(p)

    def test_plot_source_is_cell_based_and_has_pyvista(self):
        content=SCRIPT.read_text()
        self.assertIn('render_pyvista',content)
        self.assertIn('mesh_centres',content)
        self.assertIn('negative_ux_cells',content)

    def test_native_geometry_and_cells(self):
        path=os.getenv('CFDX_V84_H5')
        if not path:self.skipTest('Set CFDX_V84_H5')
        c=M.mesh_centres(Path(path))
        self.assertEqual(c.shape,(20540,3))
        self.assertTrue(np.all(np.isfinite(c)))

    def test_no_native_commit_no_plots(self):
        p=self.write('GLOBAL_PROGRESS_TRIAL_V84 iteration=1 accepted=1\n')
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaises(RuntimeError):
                M.make_figures(p,Path(d),Path('unused.h5'),Path(d)/'figures','vtk')

if __name__=='__main__':unittest.main()