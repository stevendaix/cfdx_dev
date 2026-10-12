#!/usr/bin/env python3
"""Import a NASA/TMR NACA0012 PLOT3D 2-D C-grid into CFDX-HDF5."""
from __future__ import annotations
import argparse
from pathlib import Path
import re
import h5py
import numpy as np
from generate_2d_validation_meshes import write_mesh

TMR_URL="https://turbmodels.larc.nasa.gov/naca0012_grids.html"
TMR_ARCHIVE_URL="https://www.nasa.gov/wp-content/uploads/2026/02/naca0012-grids.zip"

def read_plot3d_2d(path: Path):
    tokens=path.read_text(encoding="ascii").replace(","," ").split()
    if len(tokens)<2: raise ValueError(f"PLOT3D file is empty: {path}")
    ni,nj=int(tokens[0]),int(tokens[1])
    expected=2+2*ni*nj
    if len(tokens)!=expected: raise ValueError(f"expected {expected} numeric tokens for {ni}x{nj}, got {len(tokens)}")
    values=np.asarray([float(v) for v in tokens[2:]],dtype=np.float64)
    return values[:ni*nj].reshape((nj,ni)),values[ni*nj:].reshape((nj,ni))

def convert(x,y):
    nj,ni=x.shape; n=ni*nj
    points=[(float(x[j,i]),float(y[j,i]),0.0) for j in range(nj) for i in range(ni)]
    points += [(float(x[j,i]),float(y[j,i]),1.0) for j in range(nj) for i in range(ni)]
    def p(j,i,layer=0): return layer*n+j*ni+i
    cells=[]
    for j in range(nj-1):
        for i in range(ni-1):
            a,b,c,d=p(j,i),p(j,i+1),p(j+1,i+1),p(j+1,i)
            cells.append((a,b,c,d,a+n,b+n,c+n,d+n))
    airfoil=[(p(0,i+1),p(0,i),p(0,i,1),p(0,i+1,1)) for i in range(ni-1)]
    farfield=[(p(nj-1,i),p(nj-1,i+1),p(nj-1,i+1,1),p(nj-1,i,1)) for i in range(ni-1)]
    wake_lower=[(p(j,0),p(j+1,0),p(j+1,0,1),p(j,0,1)) for j in range(nj-1)]
    wake_upper=[(p(j+1,ni-1),p(j,ni-1),p(j,ni-1,1),p(j+1,ni-1,1)) for j in range(nj-1)]
    front=[tuple(c[:4]) for c in cells]
    back=[tuple(reversed(tuple(v+n for v in c[:4]))) for c in cells]
    return points,cells,{"airfoil":airfoil,"farfield":farfield,"wake_lower":wake_lower,"wake_upper":wake_upper,"_front_cells":front,"_back_cells":back}

def infer_dimensions(path: Path):
    m=re.search(r"_(\d+)-(\d+)\.p2dfmt(?:\.gz)?$",path.name)
    if not m: raise ValueError(f"cannot infer TMR grid dimensions from {path.name}")
    return int(m.group(1)),int(m.group(2))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--input",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    args=ap.parse_args()
    ni_name,nj_name=infer_dimensions(args.input)
    x,y=read_plot3d_2d(args.input)
    nj,ni=x.shape
    if (ni,nj)!=(ni_name,nj_name): raise ValueError(f"filename declares {ni_name}x{nj_name}, file contains {ni}x{nj}")
    points,cells,boundaries=convert(x,y)
    write_mesh(args.output,points,cells,boundaries)
    with h5py.File(args.output,"r+") as h5:
        h5.attrs["validation_source"]="NASA TMR NACA0012"
        h5.attrs["validation_source_url"]=TMR_URL
        h5.attrs["validation_archive_url"]=TMR_ARCHIVE_URL
        h5.attrs["validation_grid"]=f"{ni}x{nj}"
        h5.attrs["validation_topology"]="C-grid"
        h5.attrs["validation_physics"]="N9-S8 laminar incompressible campaign"
        h5.attrs["validation_status"]="imported-provenance-only"
    print(f"Imported {args.input} -> {args.output} ({ni}x{nj})")
    return 0

if __name__=="__main__": raise SystemExit(main())
