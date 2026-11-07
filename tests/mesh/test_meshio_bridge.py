#!/usr/bin/env python3
import pathlib, subprocess, sys, tempfile
import h5py
import meshio
import numpy as np

bridge = pathlib.Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as d:
    root = pathlib.Path(d)
    src = root / "two_tetra.msh"
    dst = root / "mesh.h5"
    points = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]], dtype=float)
    cells = [("tetra", np.array([[0,1,2,3]], dtype=int))]
    mesh = meshio.Mesh(points, cells)
    meshio.write(src, mesh, file_format="gmsh", binary=False)
    subprocess.check_call([sys.executable, str(bridge), str(src), str(dst)])
    with h5py.File(dst, "r") as h5:
        assert h5.attrs["format"] == "CFDX-HDF5-mesh-v1"
        assert tuple(h5["points"].shape) == (4, 3)
        assert h5["cell_offsets"].shape[0] == 2
        assert h5["face_offsets"].shape[0] == 5
        assert h5["neighbour"].shape[0] == 4
    # A 2D Gmsh surface mesh stores its boundary conditions as physical
    # line elements. The importer must map those lines onto reconstructed
    # edge faces instead of reporting them as skipped cell types.
    src_2d = root / "two_triangles.msh"
    dst_2d = root / "mesh_2d.h5"
    points_2d = np.array(
        [[0, 0], [1, 0], [1, 1], [0, 1]],
        dtype=float,
    )
    cells_2d = [
        ("line", np.array([[0, 1], [1, 2], [2, 3], [3, 0]], dtype=int)),
        ("triangle", np.array([[0, 1, 2], [0, 2, 3]], dtype=int)),
    ]
    physical_data = [
        np.array([1, 2, 3, 4], dtype=int),
        np.array([5, 5], dtype=int),
    ]
    field_data = {
        "inlet": np.array([1, 1]),
        "outlet": np.array([2, 1]),
        "top_wall": np.array([3, 1]),
        "bottom_wall": np.array([4, 1]),
        "fluid": np.array([5, 2]),
    }
    mesh_2d = meshio.Mesh(
        points_2d,
        cells_2d,
        cell_data={"gmsh:physical": physical_data},
        field_data=field_data,
    )
    meshio.write(src_2d, mesh_2d, file_format="gmsh", binary=False)
    subprocess.check_call([sys.executable, str(bridge), str(src_2d), str(dst_2d)])
    with h5py.File(dst_2d, "r") as h5:
        patch_metadata = h5.attrs["boundary_patches"]
        if isinstance(patch_metadata, bytes):
            patch_metadata = patch_metadata.decode("ascii")
        for name in ("inlet", "outlet", "top_wall", "bottom_wall"):
            assert name in patch_metadata, (name, patch_metadata)
        patch_ids = np.asarray(h5["patch_face_ids"])
        patch_offsets = np.asarray(h5["patch_face_offsets"])
        assert patch_ids.size == 4, patch_ids
        assert patch_offsets.tolist() == [0, 1, 2, 3, 4], patch_offsets
        assert np.count_nonzero(np.asarray(h5["neighbour"]) < 0) == 4

print("meshio bridge test passed")
