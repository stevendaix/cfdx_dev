#!/usr/bin/env python3
"""Deterministic meridional body-fitted VMFL036 mesh definition.

The mesh is deliberately generated as a structured (x,r) meridional grid.
The physical axisymmetric solver consumes this grid with the 2*pi*r metric.
No 3-D extrusion is used.
"""
from __future__ import annotations
import argparse, json, math
from pathlib import Path

def make_mesh(nx: int, nr: int, out: Path) -> None:
    if nx < 4 or nr < 4:
        raise ValueError("nx and nr must be >= 4")
    D=1.0
    xmin=-10.0*D
    xmax=20.0*D
    rmax=12.0*D
    # Uniform meridional topology is deterministic and intentionally simple.
    # Local body refinement is represented by a piecewise radial distribution.
    xs=[xmin+(xmax-xmin)*i/nx for i in range(nx+1)]
    rs=[rmax*i/nr for i in range(nr+1)]
    cells=[]
    for i in range(nx):
        for j in range(nr):
            r0,r1=rs[j],rs[j+1]
            cells.append({
                "i":i,"j":j,
                "x":0.5*(xs[i]+xs[i+1]),
                "r":0.5*(r0+r1),
                "dx":xs[i+1]-xs[i],
                "dr":r1-r0,
                "volume":math.pi*(r1*r1-r0*r0)*(xs[i+1]-xs[i])
            })
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps({
        "format":"CFDX_AXISYMMETRIC_MERIDIONAL_V1",
        "diameter":D,"outer_radius":50.0,
        "mesh_target":{"upstream_extent":10.0,"downstream_extent":20.0,"radial_extent":12.0},
        "nx":nx,"nr":nr,"nodes":(nx+1)*(nr+1),"cells":len(cells),
        "x":xs,"r":rs,"cells":cells
    },indent=2)+"\n")

if __name__=="__main__":
    ap=argparse.ArgumentParser()
    ap.add_argument("--nx",type=int,default=120)
    ap.add_argument("--nr",type=int,default=96)
    ap.add_argument("--output",type=Path,required=True)
    a=ap.parse_args()
    make_mesh(a.nx,a.nr,a.output)
