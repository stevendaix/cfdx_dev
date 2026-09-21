#!/usr/bin/env python3
"""Extract a named meshio cell field into CFDX field HDF5."""
from __future__ import annotations
import argparse
import pathlib
import sys
import h5py
import meshio
import numpy as np

VOLUME_TYPES={"tetra","hexahedron","wedge","pyramid","voxel","tetra10","hexahedron20","hexahedron27",
              "wedge15","wedge18","pyramid13","pyramid14","polyhedron"}

def main() -> int:
    p=argparse.ArgumentParser()
    p.add_argument("input",type=pathlib.Path)
    p.add_argument("output",type=pathlib.Path)
    p.add_argument("field")
    a=p.parse_args()
    try:
        mesh=meshio.read(a.input)
        values=[]
        found=False
        cell_data=getattr(mesh,"cell_data",{})
        if a.field in cell_data:
            for block, arr in zip(mesh.cells, cell_data[a.field]):
                if block.type in VOLUME_TYPES:
                    x=np.asarray(arr,dtype=np.float64).reshape(-1)
                    values.extend(x.tolist())
                    found=True
        if not found and a.field in getattr(mesh,"point_data",{}):
            raise ValueError("point_data fields are not supported for cell-field import")
        if not found:
            raise KeyError(f"cell field {a.field!r} not found")
        if not values:
            raise ValueError(f"cell field {a.field!r} contains no supported volume values")
        with h5py.File(a.output,"w") as h:
            grp=h.create_group("fields")
            grp.attrs["name"]=a.field
            grp.attrs["unit"]=""
            grp.attrs["dimension"]="1"
            grp.create_dataset("values",data=np.asarray(values,dtype=np.float64))
        return 0
    except Exception as exc:
        print(f"CFDX scalar import failed: {exc}",file=sys.stderr)
        return 2

if __name__=="__main__":
    raise SystemExit(main())
