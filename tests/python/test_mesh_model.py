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
    with pytest.raises(ValueError,match="mandatory topology"): read_mesh_catalog(path)\n\ndef test_patch_selection_has_stable_identity_independent_of_index(tmp_path):\n    path=tmp_path/"mesh.h5"; write_mesh(path)\n    mesh=read_mesh_catalog(path)\n    assert mesh.patch(0).stable_id=="patch:inlet"\n    assert mesh.patch(0).selection.stable_id=="patch:inlet"\n\n\ndef test_catalog_rejects_duplicate_patch_names(tmp_path):\n    path=tmp_path/"mesh.h5"; write_mesh(path)\n    with h5py.File(path, "a") as h:\n        h.attrs["boundary_patches"]="inlet:0:2:1;inlet:2:4:0"\n    with pytest.raises(ValueError, match="unique"):\n        read_mesh_catalog(path)\n