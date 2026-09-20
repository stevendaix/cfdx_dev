#!/usr/bin/env python3
"""
Génération des maillages carrés pour la validation MMS de Poisson 2D.
Usage : python3 scripts/generate_mms_meshes.py --output-dir /tmp/mms_meshes
"""
import gmsh, os, argparse, math, sys

def generate_square_mesh(n_cells_per_side, output_path):
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 1)
    gmsh.model.add("mms_square")
    p0 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0)
    p1 = gmsh.model.geo.addPoint(1.0, 0.0, 0.0)
    p2 = gmsh.model.geo.addPoint(1.0, 1.0, 0.0)
    p3 = gmsh.model.geo.addPoint(0.0, 1.0, 0.0)
    l_bottom = gmsh.model.geo.addLine(p0, p1)
    l_right  = gmsh.model.geo.addLine(p1, p2)
    l_top    = gmsh.model.geo.addLine(p2, p3)
    l_left   = gmsh.model.geo.addLine(p3, p0)
    loop = gmsh.model.geo.addCurveLoop([l_bottom, l_right, l_top, l_left])
    surface = gmsh.model.geo.addPlaneSurface([loop])
    gmsh.model.geo.mesh.setTransfiniteCurve(l_bottom, n_cells_per_side + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_right,  n_cells_per_side + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_top,    n_cells_per_side + 1)
    gmsh.model.geo.mesh.setTransfiniteCurve(l_left,   n_cells_per_side + 1)
    gmsh.model.geo.mesh.setTransfiniteSurface(surface)
    gmsh.model.setPhysicalName(1, l_left, "left")
    gmsh.model.setPhysicalName(1, l_right, "right")
    gmsh.model.setPhysicalName(1, l_bottom, "bottom")
    gmsh.model.setPhysicalName(1, l_top, "top")
    gmsh.model.addPhysicalGroup(2, [surface], 10)
    gmsh.model.setPhysicalName(2, 10, "fluid")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_path)
    gmsh.finalize()
    n_cells = n_cells_per_side * n_cells_per_side
    print(f"[OK] {output_path} : {n_cells} cellules")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", default="/tmp/mms_meshes")
    parser.add_argument("--sizes", nargs="+", type=int, default=[16, 32, 64])
    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    for n in args.sizes:
        generate_square_mesh(n, os.path.join(args.output_dir, f"mms_square_{n}x{n}.msh"))

if __name__ == "__main__":
    main()
