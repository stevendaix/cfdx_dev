#!/usr/bin/env python3
"""Generate CFDX 2-D qualification meshes without Gmsh.

The solver core is three-dimensional, so each 2-D case is represented by a
single hexahedral layer in z.  The mesh topology itself is generated directly
in Python and written in the CFDX-HDF5 mesh interchange format.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import h5py
import numpy as np


def fnv1a(values: np.ndarray, h: int = 1469598103934665603) -> int:
    for b in np.asarray(values).tobytes(order="C"):
        h ^= b
        h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h


def write_mesh(path: Path, points, cells, boundary_faces):
    """Write a structured hexahedral mesh to the CFDX mesh-v1 HDF5 schema."""
    points = np.asarray(points, dtype=np.float64)
    faces = []
    owner = []
    neighbour = []
    cell_faces = []
    face_map = {}

    # Local hexahedron faces. Orientation is corrected below from owner-cell
    # centroids, so the generator remains independent of block orientation.
    templates = (
        (0, 3, 2, 1), (4, 5, 6, 7),
        (0, 1, 5, 4), (1, 2, 6, 5),
        (2, 3, 7, 6), (3, 0, 4, 7),
    )

    for ci, cell in enumerate(cells):
        refs = []
        for local in templates:
            face = [int(cell[i]) for i in local]
            key = tuple(sorted(face))
            fid = face_map.get(key)
            if fid is None:
                fid = len(faces)
                face_map[key] = fid
                faces.append(face)
                owner.append(ci)
                neighbour.append(-1)
            else:
                if neighbour[fid] != -1:
                    raise ValueError(f"non-manifold face {key}")
                neighbour[fid] = ci
            refs.append(fid)
        cell_faces.append(refs)

    # Orient boundary/interior faces outward from the owner.
    cell_centres = np.asarray(
        [points[np.asarray(c)].mean(axis=0) for c in cells], dtype=float
    )
    for fid, face in enumerate(faces):
        p = points[np.asarray(face)]
        centre = p.mean(axis=0)
        rel = p - centre
        sf = 0.5 * np.sum(np.cross(rel, np.roll(rel, -1, axis=0)), axis=0)
        if np.dot(sf, centre - cell_centres[owner[fid]]) < 0.0:
            faces[fid].reverse()

    fv = np.asarray([v for f in faces for v in f], dtype=np.uint64)
    fo = np.zeros(len(faces) + 1, dtype=np.uint64)
    for i, f in enumerate(faces):
        fo[i + 1] = fo[i] + len(f)
    own = np.asarray(owner, dtype=np.uint64)
    nei = np.asarray(neighbour, dtype=np.int64)
    cf = np.asarray([f for c in cell_faces for f in c], dtype=np.uint64)
    co = np.zeros(len(cell_faces) + 1, dtype=np.uint64)
    for i, c in enumerate(cell_faces):
        co[i + 1] = co[i] + len(c)

    face_lookup = {tuple(sorted(f)): i for i, f in enumerate(faces)}
    patch_ids = []
    patch_offsets = [0]
    patch_meta = []
    for name, rows in boundary_faces.items():
        ids = sorted({face_lookup[tuple(sorted(r))] for r in rows})
        patch_ids.extend(ids)
        patch_offsets.append(len(patch_ids))
        ptype = {"inlet": 1, "outlet": 2, "wall": 0, "airfoil": 0,
                 "farfield": 7, "front": 6, "back": 6}.get(name, 7)
        patch_meta.append(f"{name}:{ids[0] if ids else 0}:{len(ids)}:{ptype}")

    topo_hash = fnv1a(fv)
    topo_hash = fnv1a(fo, topo_hash)
    topo_hash = fnv1a(own, topo_hash)
    topo_hash = fnv1a(nei, topo_hash)
    topo_hash = fnv1a(cf, topo_hash)
    topo_hash = fnv1a(co, topo_hash)
    mesh_hash = fnv1a(points[:, 0].astype(np.float64), topo_hash)
    mesh_hash = fnv1a(points[:, 1].astype(np.float64), mesh_hash)
    mesh_hash = fnv1a(points[:, 2].astype(np.float64), mesh_hash)

    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as h:
        h.attrs["format"] = "CFDX-HDF5-mesh-v1"
        h.attrs["format_version"] = "1"
        h.attrs["schema_version"] = "1"
        h.attrs["cfdx_version"] = "0.7"
        h.attrs["topology_hash"] = f"{topo_hash:016x}"
        h.attrs["mesh_hash"] = f"{mesh_hash:016x}"
        h.attrs["n_points"] = str(len(points))
        h.attrs["n_faces"] = str(len(faces))
        h.attrs["n_cells"] = str(len(cells))
        h.attrs["boundary_patches"] = ";".join(patch_meta)
        h.create_dataset("points", data=points)
        h.create_dataset("face_vertices", data=fv)
        h.create_dataset("face_offsets", data=fo)
        h.create_dataset("owner", data=own)
        h.create_dataset("neighbour", data=nei)
        h.create_dataset("cell_faces", data=cf)
        h.create_dataset("cell_offsets", data=co)
        h.create_dataset("patch_face_ids", data=np.asarray(patch_ids, dtype=np.uint64))
        h.create_dataset("patch_face_offsets", data=np.asarray(patch_offsets, dtype=np.uint64))


def extrude_2d(nodes2d, cells2d, boundary2d):
    n = len(nodes2d)
    points = [[x, y, 0.0] for x, y in nodes2d] + [[x, y, 1.0] for x, y in nodes2d]
    cells = []
    for a, b, c, d in cells2d:
        cells.append((a, b, c, d, a + n, b + n, c + n, d + n))

    bf = {}
    for name, edges in boundary2d.items():
        if name.startswith("_"):
            continue
        bf[name] = []
        for a, b in edges:
            bf[name].append((a, b, b + n, a + n))
    front_cells = boundary2d.get("_front_cells", cells2d)
    bf["front"] = [tuple(cell) for cell in front_cells]
    bf["back"] = [tuple(reversed(tuple(v + n for v in cell))) for cell in front_cells]
    return points, cells, {k: v for k, v in bf.items() if not k.startswith("_")}


def channel(nx: int, ny: int, length: float = 4.0):
    nodes = []
    for j in range(ny + 1):
        y = j / ny
        for i in range(nx + 1):
            nodes.append((length * i / nx, y))
    idx = lambda i, j: j * (nx + 1) + i
    cells = []
    for j in range(ny):
        for i in range(nx):
            cells.append((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)))

    inlet = [(idx(0, j), idx(0, j + 1)) for j in range(ny)]
    outlet = [(idx(nx, j + 1), idx(nx, j)) for j in range(ny)]
    bottom = [(idx(i + 1, 0), idx(i, 0)) for i in range(nx)]
    top = [(idx(i, ny), idx(i + 1, ny)) for i in range(nx)]
    front = [(idx(0, 0), idx(1, 0))]
    return extrude_2d(nodes, cells, {
        "inlet": inlet, "outlet": outlet, "bottom": bottom, "top": top,
        "_front": front,
    })


def bfs(n: int):
    # VMFL064 / Armaly geometry, nondimensionalized by the step height s.
    H_in = 5.2 / 4.9
    H_out = H_in + 1.0
    L_up = 200.0 / 4.9
    L_down = 100.0 / 4.9
    nx_up = 5 * n
    nx_down = max(2 * n, int(round(nx_up * L_down / L_up)))
    nodes, key, cells = [], {}, []
    def nid(x, y):
        k = (round(x, 14), round(y, 14))
        if k not in key: key[k] = len(nodes); nodes.append((x, y))
        return key[k]
    def block(x0, x1, y0, y1, nx, ny):
        ids = [[nid(x0+(x1-x0)*i/nx, y0+(y1-y0)*j/ny) for i in range(nx+1)] for j in range(ny+1)]
        for j in range(ny):
            for i in range(nx): cells.append((ids[j][i], ids[j][i+1], ids[j+1][i+1], ids[j+1][i]))
        return ids
    up = block(-L_up, 0.0, 1.0, 1.0+H_in, nx_up, n)
    lo = block(0.0, L_down, 0.0, 1.0, nx_down, n)
    hi = block(0.0, L_down, 1.0, H_out, nx_down, n)
    inlet = [(up[j+1][0], up[j][0]) for j in range(n)]
    outlet = [(hi[j+1][-1], hi[j][-1]) for j in range(n)]
    wall = []
    wall += [(up[0][i], up[0][i+1]) for i in range(nx_up)]
    wall += [(lo[0][i+1], lo[0][i]) for i in range(nx_down)]
    wall += [(hi[-1][i], hi[-1][i+1]) for i in range(nx_down)]
    wall += [(lo[j][0], lo[j+1][0]) for j in range(n)]
    return extrude_2d(nodes, cells, {'inlet': inlet, 'outlet': outlet, 'wall': wall, '_front_cells': cells})

def naca0012_surface(n: int):
    t = 0.12
    beta = np.linspace(0.0, math.pi, n + 1)
    x = 0.5 * (1.0 - np.cos(beta))
    yt = 5.0 * t * (
        0.2969 * np.sqrt(x) - 0.1260*x - 0.3516*x*x
        + 0.2843*x**3 - 0.1015*x**4
    )
    # Closed polygon uses the finite NACA 0012 trailing-edge thickness.
    upper = list(zip(x[::-1], yt[::-1]))  # TE upper -> LE
    lower = list(zip(x[1:], -yt[1:]))    # after LE -> TE lower
    return upper + lower


def naca_o_grid(n_surface: int, n_radial: int, radius: float = 20.0):
    surf = naca0012_surface(n_surface)
    m = len(surf)
    # Parameterize the outer boundary by the same closed-loop index. The
    # surface itself is not a mesh oracle; this is only a deterministic,
    # Gmsh-free body-fitted starting mesh.
    pts = []
    for j in range(n_radial + 1):
        eta = j / n_radial
        rfac = eta ** 1.8
        for i, (x, y) in enumerate(surf):
            theta = 2.0 * math.pi * i / m
            xo, yo = radius * math.cos(theta), radius * math.sin(theta)
            pts.append(((1-rfac)*x + rfac*xo, (1-rfac)*y + rfac*yo))
    def P(j,i): return j*m + (i % m)
    cells2d = []
    for j in range(n_radial):
        for i in range(m):
            cells2d.append((P(j,i),P(j,i+1),P(j+1,i+1),P(j+1,i)))

    airfoil = [(P(0,i+1),P(0,i)) for i in range(m)]
    farfield = [(P(n_radial,i),P(n_radial,i+1)) for i in range(m)]
    # Periodic seams are interior topology because i is wrapped; the only
    # physical boundaries are the airfoil and farfield.
    nodes2d = pts
    return extrude_2d(nodes2d, cells2d, {
        "airfoil": airfoil, "farfield": farfield,
        "_front_cells": cells2d,
    })


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=("channel", "bfs", "naca0012", "all"))
    ap.add_argument("--output-dir", type=Path, default=Path("build/validation_meshes"))
    ap.add_argument("--quick", action="store_true")
    args = ap.parse_args()

    out = args.output_dir
    if args.case in ("channel", "all"):
        for n in ((16,) if args.quick else (16,32,64)):
            p, c, b = channel(4*n, n)
            write_mesh(out / f"channel_n{n}.h5", p, c, b)
    if args.case in ("bfs", "all"):
        levels = (16,) if args.quick else (16,32,64)
        for n in levels:
            p, c, b = bfs(n)
            write_mesh(out / f"bfs_re200_n{n}.h5", p, c, b)
    if args.case in ("naca0012", "all"):
        levels = ((64,24),) if args.quick else ((64,24),(128,48),(256,96))
        for n, nr in levels:
            p, c, b = naca_o_grid(n, nr, 20.0)
            write_mesh(out / f"naca0012_laminar_n{n}.h5", p, c, b)
    print(f"validation meshes written to {out}")


if __name__ == "__main__":
    main()
