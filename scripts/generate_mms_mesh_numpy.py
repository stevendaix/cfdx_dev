#!/usr/bin/env python3
"""
Génération maillage 2D structuré en NumPy (durable - sans Gmsh).
Écrit un .msh v2.2 avec éléments 1D (frontières) EXPLICITES.
Usage:
  python3 scripts/generate_mms_mesh_numpy.py --n 16 --output /tmp/mms_meshes/mms_square_16x16.msh
"""
import numpy as np
import argparse
import os


def generate_structured_mesh_2d(n: int, output_path: str) -> None:
    # ============================================================
    # 1. POINTS [0,1]² : (n+1) × (n+1) sommets
    # ============================================================
    x = np.linspace(0.0, 1.0, n + 1)
    y = np.linspace(0.0, 1.0, n + 1)
    X, Y = np.meshgrid(x, y)
    points = np.column_stack([X.ravel(), Y.ravel(), np.zeros((n + 1)**2)])
    n_points = len(points)

    def pt(i, j):
        return j * (n + 1) + i + 1  # 1-based Gmsh

    # ============================================================
    # 2. CELLULES : triangles (type 2) par carré
    # ============================================================
    triangles = []  # (tags: 10 = fluid)
    for j in range(n):
        for i in range(n):
            triangles.append([pt(i, j), pt(i + 1, j), pt(i, j + 1)])
            triangles.append([pt(i + 1, j), pt(i + 1, j + 1), pt(i, j + 1)])

    n_triangles = len(triangles)

    # ============================================================
    # 3. FRONTIÈRES : lignes 1D (type 1) avec PhysicalNames
    # ============================================================
    left_lines = []     # x = 0
    right_lines = []    # x = 1
    bottom_lines = []   # y = 0
    top_lines = []      # y = 1
    for k in range(n):
        left_lines.append([pt(0, k), pt(0, k + 1)])
        right_lines.append([pt(n, k), pt(n, k + 1)])
        bottom_lines.append([pt(k, 0), pt(k + 1, 0)])
        top_lines.append([pt(k, n), pt(k + 1, n)])

    n_lines = len(left_lines) + len(right_lines) + len(bottom_lines) + len(top_lines)

    # ============================================================
    # 4. ÉCRITURE .MSH v2.2
    # ============================================================
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w") as f:
        f.write("$MeshFormat\n")
        f.write("2.2 0 8\n")
        f.write("$EndMeshFormat\n")

        # PhysicalNames (tags: 1=left, 2=right, 3=bottom, 4=top, 10=fluid)
        f.write("$PhysicalNames\n")
        f.write("5\n")
        f.write('1 1 "left"\n')
        f.write('1 2 "right"\n')
        f.write('1 3 "bottom"\n')
        f.write('1 4 "top"\n')
        f.write('2 10 "fluid"\n')
        f.write("$EndPhysicalNames\n")

        # Nodes (1-based Gmsh)
        f.write("$Nodes\n")
        f.write(f"{n_points}\n")
        for idx in range(n_points):
            f.write(f"{idx + 1} {points[idx][0]:.10f} {points[idx][1]:.10f} {points[idx][2]:.10f}\n")
        f.write("$EndNodes\n")

        # Elements
        n_elements = n_lines + n_triangles
        f.write("$Elements\n")
        f.write(f"{n_elements}\n")

        elem_id = 1
        # Frontières (type 1 = ligne, tags physiques: 1,2,3,4)
        for line in left_lines:
            f.write(f"{elem_id} 1 1 1 {line[0]} {line[1]}\n")
            elem_id += 1
        for line in right_lines:
            f.write(f"{elem_id} 1 1 2 {line[0]} {line[1]}\n")
            elem_id += 1
        for line in bottom_lines:
            f.write(f"{elem_id} 1 1 3 {line[0]} {line[1]}\n")
            elem_id += 1
        for line in top_lines:
            f.write(f"{elem_id} 1 1 4 {line[0]} {line[1]}\n")
            elem_id += 1

        # Cellules (type 2 = triangle, tag physique: 10 = fluid)
        for tri in triangles:
            f.write(f"{elem_id} 2 1 10 {tri[0]} {tri[1]} {tri[2]}\n")
            elem_id += 1

        f.write("$EndElements\n")

    print(f"[OK] {output_path} : {n_triangles} triangles, {n_lines} lignes frontières, {n_points} points")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=16)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()
    generate_structured_mesh_2d(args.n, args.output)


if __name__ == "__main__":
    main()
