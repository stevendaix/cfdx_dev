#!/usr/bin/env python3
"""Generate a deterministic body-fitted meridional VMFL036 axisymmetric mesh.

The mesh is an annular polar grid around the spherical body. The inner
boundary is the sphere meridian and the outer boundary is the documented
circular farfield of radius 50D. No Gmsh dependency is used.
"""
from __future__ import annotations
import argparse, math
from pathlib import Path

def write_mesh(nr: int, nt: int, out: Path) -> None:
    if nr < 8 or nt < 16:
        raise ValueError("nr>=8 and nt>=16 required")
    D = 1.0
    R = 0.5 * D
    Rout = 50.0 * D
    radii = [R + (Rout - R) * (i / nr) ** 1.65 for i in range(nr + 1)]
    theta = [math.pi * j / nt for j in range(nt + 1)]
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", encoding="utf-8") as f:
        f.write("CFDX_AXISYMMETRIC_POLAR_V2\n")
        f.write(f"{nr} {nt} {D:.17g} {Rout:.17g}\n")
        for s in radii:
            for t in theta:
                f.write(f"{s*math.cos(t):.17g} {s*math.sin(t):.17g}\n")
    print(f"VMFL036_AXISYMESH: PASS nr={nr} nt={nt} cells={nr*nt} nodes={(nr+1)*(nt+1)}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--nr", type=int, default=80)
    ap.add_argument("--nt", type=int, default=120)
    ap.add_argument("--output", type=Path, required=True)
    a = ap.parse_args()
    write_mesh(a.nr, a.nt, a.output)
