"""Headless Qt/VTK smoke for the bounded, async CFDX surface LOD preview."""
from __future__ import annotations

import importlib.util
from time import monotonic, sleep

import h5py
import pytest

from test_workbench_lot9_industrial import make_cartesian_case


@pytest.mark.skipif(any(importlib.util.find_spec(name) is None for name in
    ("PySide6", "pyvista", "pyvistaqt")), reason="Qt/PyVista optional in local environment")
def test_real_qt_view_uses_sampled_surface_not_full_volume(tmp_path):
    from PySide6.QtWidgets import QApplication
    from cfdx.case_inspector import inspect_cfdx_file
    from cfdx.gui_3d import PyVistaQtView

    app = QApplication.instance() or QApplication(["cfdx-lot9-surface"])
    path=make_cartesian_case(tmp_path/"mesh.h5",[(0,.5),(.5,1)], [1,2])
    with h5py.File(path,"r+") as h:
        h.attrs["format"]="CFDX-HDF5-mesh-v1"
        h.attrs["n_cells"]=2
    summary=inspect_cfdx_file(path)
    view=PyVistaQtView()
    events=[]
    view.progressive_preview_ready.connect(lambda message: events.append(("ok",message)))
    view.progressive_preview_failed.connect(lambda message: events.append(("error",message)))
    view.load_progressive_surface(str(path),summary=summary)
    deadline=monotonic()+12
    while not events and monotonic()<deadline:
        app.processEvents(); sleep(.01)
    assert events and events[0][0]=="ok",events
    assert "NO full-volume fields" in events[0][1]
    assert view._source_dataset is None
    assert view._rendered_dataset.n_cells==12
    view.close()
    app.processEvents()


@pytest.mark.skipif(importlib.util.find_spec("vtk") is None, reason="VTK optional")
def test_sampled_surface_is_valid_native_vtk_polydata(tmp_path):
    """Headless VTK data-path test; no GPU render or OpenGL claim."""
    from vtk import vtkCellArray, vtkPoints, vtkPolyData
    from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray
    from cfdx.surface_streaming import read_surface_preview

    path=make_cartesian_case(tmp_path/"geometry.h5",[(0,1)],[2])
    surface=read_surface_preview(path)
    vtk_points=vtkPoints()
    vtk_points.SetData(numpy_to_vtk(surface.points,deep=True))
    vtk_faces=vtkCellArray()
    vtk_faces.SetCells(len(surface.original_face_ids),
                       numpy_to_vtkIdTypeArray(surface.faces,deep=True))
    poly=vtkPolyData()
    poly.SetPoints(vtk_points)
    poly.SetPolys(vtk_faces)
    assert poly.GetNumberOfPoints()==8
    assert poly.GetNumberOfPolys()==6
    assert poly.GetBounds()==(0.0,1.0,0.0,1.0,0.0,1.0)
