#!/usr/bin/env python3
"""Generate deterministic CFD mesh geometries and validate meshio imports."""
from __future__ import annotations
import pathlib, subprocess, sys, tempfile
import meshio
import numpy as np

BRIDGE = pathlib.Path(sys.argv[1]).resolve()

def cube():
    p=np.array([[x,y,z] for z in (0,1) for y in (0,1) for x in (0,1)],float)
    return p, np.array([[0,1,3,2,4,5,7,6]],int)

def two_tetra():
    p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,1]],float)
    return p, np.array([[0,1,2,3],[1,2,3,4]],int)

def pyramid():
    p=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[.5,.5,1]],float)
    return p, np.array([[0,1,2,3,4]],int)

def wedge():
    p=np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,0,1],[0,1,1]],float)
    return p, np.array([[0,1,2,3,4,5]],int)

cases={
 "cube":(cube,{"cells":1,"boundary":6,"internal":0}),
 "two_tetra":(two_tetra,{"cells":2,"boundary":6,"internal":1}),
 "pyramid":(pyramid,{"cells":1,"boundary":5,"internal":0}),
 "wedge":(wedge,{"cells":1,"boundary":5,"internal":0}),
}

for name,(factory,expected) in cases.items():
    points, data=factory()
    with tempfile.TemporaryDirectory() as d:
        root=pathlib.Path(d)
        for fmt, suffix, kwargs in [
            ("gmsh_ascii","msh",{"file_format":"gmsh","binary":False}),
            ("gmsh_binary","msh",{"file_format":"gmsh","binary":True}),
            ("vtu_ascii","vtu",{"file_format":"vtu","binary":False}),
        ]:
            src=root/f"{name}.{suffix}"
            out=root/f"{name}.h5"
            cell_type = {"cube":"hexahedron","two_tetra":"tetra","pyramid":"pyramid","wedge":"wedge"}[name]
            meshio.write(src, meshio.Mesh(points, [(cell_type, data)]), **kwargs)
            subprocess.check_call([sys.executable,str(BRIDGE),str(src),str(out)])
            import h5py
            with h5py.File(out,"r") as h:
                cells=int(h["cell_offsets"].shape[0]-1)
                faces=int(h["face_offsets"].shape[0]-1)
                internal=int(np.count_nonzero(h["neighbour"][:] >= 0))
                boundary=int(np.count_nonzero(h["neighbour"][:] == -1))
                assert cells==expected["cells"],(name,fmt,cells)
                assert internal==expected["internal"],(name,fmt,internal)
                assert boundary==expected["boundary"],(name,fmt,boundary)
                assert faces==internal+boundary
print("mesh geometry import matrix passed")
