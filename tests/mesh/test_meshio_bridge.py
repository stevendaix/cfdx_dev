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
    cells = [
        ("triangle", np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]], dtype=int)),
        ("tetra", np.array([[0,1,2,3]], dtype=int)),
    ]
    mesh = meshio.Mesh(
        points, cells,
        cell_data={"gmsh:physical": [
            np.array([1,1,1,1], dtype=int),
            np.array([2], dtype=int),
        ]},
        field_data={
            "wall": np.array([1,2], dtype=int),
            "fluid": np.array([2,3], dtype=int),
        },
    )
    meshio.write(src, mesh, file_format="gmsh", binary=False)
    subprocess.check_call([sys.executable, str(bridge), str(src), str(dst)])
    with h5py.File(dst, "r") as h5:
        assert h5.attrs["format"] == "CFDX-HDF5-mesh-v1"
        assert tuple(h5["points"].shape) == (4, 3)
        assert h5["cell_offsets"].shape[0] == 2
        assert h5["face_offsets"].shape[0] == 5
        assert h5["neighbour"].shape[0] == 4
print("meshio bridge test passed")
