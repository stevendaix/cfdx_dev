import pathlib
import subprocess
import sys
import tempfile

import h5py
import meshio
import numpy as np

bridge = pathlib.Path(sys.argv[1]).resolve()
with tempfile.TemporaryDirectory() as d:
    root = pathlib.Path(d)
    src = root / "two_tetra.msh"
    dst = root / "mesh.h5"
    points = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], dtype=float)
    cells = [("tetra", np.array([[0, 1, 2, 3]], dtype=int))]
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
    # Use an explicit Gmsh 2.2 ASCII fixture. Meshio's Gmsh 4.1 writer
    # requires gmsh:dim_tags for mixed line/triangle blocks; that writer
    # constraint is unrelated to the importer contract tested below.
    src_2d.write_text(
        "$MeshFormat\n"
        "2.2 0 8\n"
        "$EndMeshFormat\n"
        "$PhysicalNames\n"
        "5\n"
        '1 1 "inlet"\n'
        '1 2 "outlet"\n'
        '1 3 "top_wall"\n'
        '1 4 "bottom_wall"\n'
        '2 5 "fluid"\n'
        "$EndPhysicalNames\n"
        "$Nodes\n"
        "4\n"
        "1 0 0 0\n"
        "2 1 0 0\n"
        "3 1 1 0\n"
        "4 0 1 0\n"
        "$EndNodes\n"
        "$Elements\n"
        "6\n"
        "1 1 2 1 1 1 2\n"
        "2 1 2 2 1 2 3\n"
        "3 1 2 3 1 3 4\n"
        "4 1 2 4 1 4 1\n"
        "5 2 2 5 1 1 2 3\n"
        "6 2 2 5 1 1 3 4\n"
        "$EndElements\n",
        encoding="utf-8",
    )
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

    # SU2 marker names are strings. The importer normalizes them for
    # meshio, then restores the names as CFDX physical boundary patches.
    su2_source = root / "named_markers.su2"
    su2_output = root / "named_markers.h5"
    su2_source.write_text(
        "NDIME= 2\n"
        "NPOIN= 4\n"
        "0.0 0.0\n"
        "1.0 0.0\n"
        "1.0 1.0\n"
        "0.0 1.0\n"
        "NELEM= 2\n"
        "5 0 1 2\n"
        "5 0 2 3\n"
        "NMARK= 4\n"
        "MARKER_TAG= inlet\n"
        "MARKER_ELEMS= 1\n"
        "3 0 1\n"
        "MARKER_TAG= outlet\n"
        "MARKER_ELEMS= 1\n"
        "3 1 2\n"
        "MARKER_TAG= top_wall\n"
        "MARKER_ELEMS= 1\n"
        "3 2 3\n"
        "MARKER_TAG= bottom_wall\n"
        "MARKER_ELEMS= 1\n"
        "3 3 0\n",
        encoding="utf-8",
    )
    result = subprocess.run(
        [sys.executable, str(bridge), str(su2_source), str(su2_output)],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "tags of string type" not in result.stderr, result.stderr
    with h5py.File(su2_output, "r") as h5:
        patch_metadata = h5.attrs["boundary_patches"]
        if isinstance(patch_metadata, bytes):
            patch_metadata = patch_metadata.decode("ascii")
        for name in ("inlet", "outlet", "top_wall", "bottom_wall"):
            assert name in patch_metadata, (name, patch_metadata)
        assert np.asarray(h5["patch_face_ids"]).size == 4
        assert np.asarray(h5["patch_face_offsets"]).tolist() == [0, 1, 2, 3, 4]

print("meshio bridge test passed")
