"""V78: refuse un diagnostic géométrique présenté comme validation."""
import sys
import unittest
import tempfile
import json
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
from bfs_geometry_comparability_v78 import gate_reference_geometry, pressure_failure_evidence

class TestBfsReferenceGateV78(unittest.TestCase):
    def setUp(self):
        self.ds=dict(h_m=.0127,h_upstream=8.,h_downstream=9.,x_min_h=-130.,x_max_h=50.)
    def test_openfoam_six_block_geometry_only(self):
        v=gate_reference_geometry(self.ds)
        self.assertTrue(v['pass_reference_domain'])
    def test_pilot_rejected_upstream_downstream(self):
        d=dict(self.ds,h_upstream=7.,h_downstream=8.,x_min_h=-6,x_max_h=40)
        v=gate_reference_geometry(d)
        self.assertFalse(v['pass_reference_domain'])
        self.assertFalse(v['upstream_height_8H']);self.assertFalse(v['downstream_height_9H'])
    def test_short_inlet_not_nasa_uniform_inflow(self):
        d=dict(self.ds,x_min_h=-4.)
        self.assertFalse(gate_reference_geometry(d)['long_inlet_for_uniform_profile'])
    def test_small_geometry_error_rejected(self):
        d=dict(self.ds,h_upstream=7.999)
        self.assertFalse(gate_reference_geometry(d)['pass_reference_domain'])
    def test_nonfinite_rejected(self):
        d=dict(self.ds,h_upstream=float('nan'))
        with self.assertRaises(ValueError):gate_reference_geometry(d)
    def test_real_pressure_rejections_not_summary_duplicate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root=Path(tmp); log=root/'native.log';e=root/'execution.json'
            line='CFDX nonlinear iteration 1 rejected; retry=1 reason=solve_steady_incompressible: pressure-correction solve did not converge (status=1, iterations=1000, residual=1e-12, relative=5e-10)'
            log.write_text(line+'\n'+line.replace('retry=1','retry=2')+'\nCFDX production solver error: pressure-correction solve did not converge (status=1, iterations=1000, residual=1e-12, relative=5e-10)\n')
            e.write_text(json.dumps({'converged':False,'process_exit_code':2}))
            rows=pressure_failure_evidence(log,e)['trials']
            self.assertEqual(len(rows),2)
            e.write_text(json.dumps({'converged':True,'process_exit_code':0}))
            with self.assertRaises(ValueError):pressure_failure_evidence(log,e)

    def test_incomplete_rejected(self):
        with self.assertRaises(ValueError):gate_reference_geometry({'h_m':.0127})

if __name__=='__main__':unittest.main()
