@@ -0,0 +1,65 @@
"""Diagnostic-only CFDX VTU postprocessing invariants for V77."""
import pathlib
import sys
import unittest
import numpy as np
import vtk

sys.path.insert(0,str(pathlib.Path(__file__).resolve().parents[2]/'scripts'))
from bfs_native_diagnostic_comparison_v77 import (wls_uv_proxy, interpolated_profiles,
    require, compare_metrics, audited_cell_centers)

class NativeBfsV77Tests(unittest.TestCase):
    def test_native_hex_centers_from_original_connectivity(self):
        points=vtk.vtkPoints()
        for p in ((0,0,0),(2,0,0),(2,1,0),(0,1,0),
                  (0,0,1),(2,0,1),(2,1,1),(0,1,1)):
            points.InsertNextPoint(*p)
        grid=vtk.vtkUnstructuredGrid();grid.SetPoints(points)
        hexa=vtk.vtkHexahedron()
        for i in range(8):hexa.GetPointIds().SetId(i,i)
        grid.InsertNextCell(hexa.GetCellType(),hexa.GetPointIds())
        mesh={'ncell':1,'points':np.array([points.GetPoint(i) for i in range(8)]),
              'faceoffset':np.array([0,4,8,12,16,20,24]),
              'vertices':np.array([0,1,2,3,4,5,6,7,0,1,5,4,1,2,6,5,
                                   2,3,7,6,3,0,4,7]),
              'cell_faces':np.arange(6), 'cell_offsets':np.array([0,6])}
        center=audited_cell_centers(grid,mesh)
        self.assertTrue(np.allclose(center,[[1.,.5,.5]]))
        mesh['vertices'][mesh['vertices']==0]=7
        with self.assertRaises(ValueError):audited_cell_centers(grid,mesh)

    def test_wls_affine_shear_exact(self):
        xx,yy=np.meshgrid(np.arange(5)*0.01,np.arange(5)*0.002)
        x=xx.ravel();y=yy.ravel()
        xyz=np.column_stack([x,y,np.zeros(len(x))]);nu=2e-3
        # Ux=2x+3y and Uy=-4x+5y -> Sxy=3-4=-1
        fields=dict(u_x=2*x+3*y,u_y=-4*x+5*y,nu_t=np.full(len(x),nu))
        shear,cond=wls_uv_proxy(xyz,fields,10)
        self.assertTrue(np.isfinite(shear).all())
        self.assertTrue(np.isfinite(cond).all())
        self.assertTrue(np.allclose(shear,1000*nu/(10**2),atol=1e-10))

    def test_interpolator_does_not_extrapolate_to_wall(self):
        xyz=np.array([[0.02,0.002,0],[0.02,0.012,0],
                      [0.04,0.002,0],[0.04,0.012,0]],float)
        fields=dict(u_x=np.array([1.,2.,1.,2.]))
        samples=[dict(xH='3',yH='0.6',uUref='1',uv_times1000='-1'),
                 dict(xH='3',yH='0.05',uUref='1',uv_times1000='-1')]
        data=interpolated_profiles(xyz,fields,0.01,10,samples,np.array([1.,2.,3.,4.]))
        self.assertIsNotNone(data[0]['cfdx_u_interp'])
        self.assertIsNone(data[1]['cfdx_u_interp'])
        self.assertIsNone(data[1]['cfdx_uv_wls_proxy'])

    def test_fails_closed_on_missing_geometry(self):
        with self.assertRaises(ValueError):require(False,'missing case')

    def test_pressure_metrics_diagnostic_not_qualification(self):
        row=[dict(patch='downstream_wall',x_over_h=float(x),cp_owner_proxy_shifted=0.0)
             for x in (3,5,7,9)]
        metrics=compare_metrics(row,[])
        self.assertEqual(metrics['cp_faces_evaluated'],4)
        self.assertTrue(metrics['cp_owner_proxy_rmse_vs_exp_on_CFDX_wall_face_x']>0)
        self.assertIsNone(metrics['u_station_4']['rmse_u_u_ref'])

if __name__=='__main__': unittest.main()
