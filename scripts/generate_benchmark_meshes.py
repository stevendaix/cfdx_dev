#!/usr/bin/env python3
"""Génération des maillages cubiques 3D pour benchmark."""
import gmsh, os, argparse

def generate_cube_mesh(n_per_side, output_path):
    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", 0)
    gmsh.model.add("bench_cube")
    pts = []
    for z in [0.0, 1.0]:
        for y in [0.0, 1.0]:
            for x in [0.0, 1.0]:
                pts.append(gmsh.model.geo.addPoint(x, y, z))
    edges = []
    for i in range(4):
        a = pts[i]; b = pts[(i+1)%4]; edges.append(gmsh.model.geo.addLine(a, b))
    for i in range(4):
        a = pts[i+4]; b = pts[((i+1)%4)+4]; edges.append(gmsh.model.geo.addLine(a, b))
    for i in range(4):
        edges.append(gmsh.model.geo.addLine(pts[i], pts[i+4]))
    loops = []
    for face_idx, (e1, e2, e3, e4) in enumerate([(0,1,2,3), (4,5,6,7), (0,9,4,8), (2,11,6,10), (3,8,7,10), (1,9,5,11)]):
        loop = gmsh.model.geo.addCurveLoop([edges[e1], edges[e2], -edges[e3], -edges[e4]])
        surface = gmsh.model.geo.addPlaneSurface([loop])
        gmsh.model.setPhysicalName(2, surface, f"face_{face_idx}")
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(3)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.write(output_path)
    gmsh.finalize()
    n_cells = n_per_side**3
    print(f"[OK] {output_path} : {n_cells} cellules 3D")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--output-dir", default="/tmp/bench_meshes")
    p.add_argument("--sizes", nargs="+", type=int, default=[22, 47])
    a = p.parse_args()
    os.makedirs(a.output_dir, exist_ok=True)
    for n in a.sizes:
        generate_cube_mesh(n, os.path.join(a.output_dir, f"bench_cube_{n}x{n}x{n}.msh"))

if __name__ == "__main__":
    main()
