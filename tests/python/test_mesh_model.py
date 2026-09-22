from pathlib import Path
import h5py
import numpy as np
import pytest
from cfdx.mesh_model import read_mesh_catalog

def write_mesh(path: Path) -> None:
    with h5py.File(path, "w") as h:
        h.attrs["format"]="CFDX-HDF5-mesh-v1"
        h.attrs["n_cells"]="4"
        h.create_dataset("points", data=np.zeros((8,3)))
        h.create_dataset("owner", data=np.arange(6,dtype=np.uint64))
        h.create_dataset("neighbour", data=np.full(6,-1,dtype=np.int64))
        h.create_dataset("cell_offsets", data=np.arange(5,dtype=np.uint64))
        h.attrs["boundary_patches"]="inlet:0:2:1;wall:2:4:0"
        h.create_dataset("patch_face_ids",data=np.arange(6,dtype=np.uint64))
        h.create_dataset("patch_face_offsets",data=np.array([0,2,6],dtype=np.uint64))

def test_catalog_reads_real_hdf5_patch_metadata(tmp_path):
    path=tmp_path/"mesh.h5"; write_mesh(path)
    mesh=read_mesh_catalog(path)
    assert (mesh.n_points,mesh.n_faces,mesh.n_cells)==(8,6,4)
    assert mesh.patch_names==("inlet","wall")
    assert mesh.patch(0).type=="inlet"
    assert mesh.patch(0).face_count==2
    assert mesh.patch(0).selection.index==0

def test_catalog_rejects_missing_topology(tmp_path):
    path=tmp_path/"mesh.h5"
    with h5py.File(path,"w") as h: h.create_dataset("points",data=np.zeros((1,3)))
    with pytest.raises(ValueError,match="mandatory topology"): read_mesh_catalog(path)

def test_patch_selection_has_stable_identity_independent_of_index(tmp_path):
    path=tmp_path/"mesh.h5"; write_mesh(path)
    mesh=read_mesh_catalog(path)
    assert mesh.patch(0).stable_id=="patch:inlet"
    assert mesh.patch(0).selection.stable_id=="patch:inlet"


def test_catalog_rejects_duplicate_patch_names(tmp_path):
    path=tmp_path/"mesh.h5"; write_mesh(path)
    with h5py.File(path, "a") as h:
        h.attrs["boundary_patches"]="inlet:0:2:1;inlet:2:4:0"
    with pytest.raises(ValueError, match="unique"):
        read_mesh_catalog(path)


def test_cell_selection_has_stable_identity(tmp_path):
    path=tmp_path/"mesh.h5"; write_mesh(path)
    mesh=read_mesh_catalog(path)
    assert mesh.cell(2).stable_id=="cell:2"
    assert mesh.cell(2).selection.stable_id=="cell:2"
    assert len(set(mesh.cell_ids)) == mesh.n_cells


def test_catalog_reads_multi_region_hdf5(tmp_path):
    path=tmp_path/"regions.h5"
    with h5py.File(path, "w") as h:
        regions = h.create_group("regions")
        for region_name, n_cells, patch_name in (("fluid", 3, "inlet"), ("solid", 2, "wall")):
            g = regions.create_group(region_name)
            g.create_dataset("points", data=np.zeros((4,3)))
            g.create_dataset("owner", data=np.arange(2,dtype=np.uint64))
            g.create_dataset("neighbour", data=np.full(2,-1,dtype=np.int64))
            g.create_dataset("cell_offsets", data=np.arange(n_cells+1,dtype=np.uint64))
            g.attrs["n_cells"] = n_cells
            g.attrs["boundary_patches"] = f"{patch_name}:0:2:1"
            g.create_dataset("patch_face_ids", data=np.arange(2,dtype=np.uint64))
            g.create_dataset("patch_face_offsets", data=np.array([0,2],dtype=np.uint64))
            g.create_dataset("cell_ids", data=np.asarray([f"{region_name}-{i}" for i in range(n_cells)], dtype="S16"))
    mesh=read_mesh_catalog(path)
    assert tuple(r.name for r in mesh.regions)==("fluid","solid")
    assert mesh.n_cells==5
    assert mesh.patch(0).stable_id=="region:fluid/patch:inlet"
    assert mesh.patch(1).stable_id=="region:solid/patch:wall"
    assert mesh.cell_ids[0]=="region:fluid/cell:fluid-0"
    assert mesh.cell_ids[-1]=="region:solid/cell:solid-1"
