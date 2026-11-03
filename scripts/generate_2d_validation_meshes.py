"""Generate CFDX 2-D qualification meshes without Gmsh.

The solver core is three-dimensional, so each 2-D case is represented by a
single hexahedral layer in z.  The mesh topology itself is generated directly
in Python and written in the CFDX-HDF5 mesh interchange format.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np


def fnv1a(values: np.ndarray, h: int = 1469598103934665603) -> int:
    for b in np.asarray(values).tobytes(order="C"):
        h ^= b
        h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h


def write_mesh(
    path: Path,
    points,
    cells,
    boundary_faces,
    patch_type_overrides=None,
):
    """Write a structured hexahedral mesh to the CFDX mesh-v1 HDF5 schema.

    ``patch_type_overrides`` supports fixtures whose physical boundary
    semantics differ from the default extrusion (e.g. thermal end caps).
    """
    points = np.asarray(points, dtype=np.float64)
    faces = []
    owner = []
    neighbour = []
    cell_faces = []
    face_map = {}

    # Local hexahedron faces. Orientation is corrected below from owner-cell
    # centroids, so the generator remains independent of block orientation.
    templates = (
        (0, 3, 2, 1),
        (4, 5, 6, 7),
        (0, 1, 5, 4),
        (1, 2, 6, 5),
        (2, 3, 7, 6),
        (3, 0, 4, 7),
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
            face.reverse()

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
        ptype = {
            "inlet": 1,
            "outlet": 2,
            "wall": 0,
            "airfoil": 0,
            "farfield": 7,
            "farfield_inlet": 1,
            "farfield_outlet": 2,
            "front": 6,
            "back": 6,
        }.get(name, 7)
        if patch_type_overrides and name in patch_type_overrides:
            ptype = int(patch_type_overrides[name])
        if ptype not in range(8):
            raise ValueError(f"invalid patch type {ptype} for patch {name}")
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
        h.create_dataset(
            "patch_face_offsets", data=np.asarray(patch_offsets, dtype=np.uint64)
        )


def extrude_2d(nodes2d, cells2d, boundary2d, layers: int = 1, span: float = 1.0):
    """Extrude an oriented 2-D quad mesh to hexahedra over *layers*.

    The default (one layer) preserves the legacy 2-D test contract. Multiple
    layers are real 3-D cells; front/back are only the two END caps, not the
    interior z-planes. A constant cross-section is an extrusion verification
    fixture, not a finite-wing validation case.
    """
    if layers < 1 or not math.isfinite(span) or span <= 0.0:
        raise ValueError("layers must be positive and span finite and positive")
    n = len(nodes2d)
    points = [[x, y, span * k / layers] for k in range(layers + 1) for x, y in nodes2d]
    cells = []
    for k in range(layers):
        lo, hi = k * n, (k + 1) * n
        for a, b, c, d in cells2d:
            cells.append(
                (a + lo, b + lo, c + lo, d + lo, a + hi, b + hi, c + hi, d + hi)
            )
    bf = {}
    for name, edges in boundary2d.items():
        if name.startswith("_"):
            continue
        bf[name] = []
        for k in range(layers):
            lo, hi = k * n, (k + 1) * n
            for a, b in edges:
                bf[name].append((a + lo, b + lo, b + hi, a + hi))
    front_cells = boundary2d.get("_front_cells", cells2d)
    bf["front"] = [tuple(cell) for cell in front_cells]
    bf["back"] = [
        tuple(reversed(tuple(v + layers * n for v in cell))) for cell in front_cells
    ]
    return points, cells, bf


def channel(nx: int, ny: int, length: float = 4.0):
    nodes = []
    for j in range(ny + 1):
        y = j / ny
        for i in range(nx + 1):
            nodes.append((length * i / nx, y))

    def idx(i, j):
        return j * (nx + 1) + i

    cells = []
    for j in range(ny):
        for i in range(nx):
            cells.append((idx(i, j), idx(i + 1, j), idx(i + 1, j + 1), idx(i, j + 1)))

    inlet = [(idx(0, j), idx(0, j + 1)) for j in range(ny)]
    outlet = [(idx(nx, j + 1), idx(nx, j)) for j in range(ny)]
    bottom = [(idx(i + 1, 0), idx(i, 0)) for i in range(nx)]
    top = [(idx(i, ny), idx(i + 1, ny)) for i in range(nx)]
    front = [(idx(0, 0), idx(1, 0))]
    return extrude_2d(
        nodes,
        cells,
        {
            "inlet": inlet,
            "outlet": outlet,
            "bottom": bottom,
            "top": top,
            "_front": front,
        },
    )


def bfs(n: int):
    # VMFL064 / Armaly geometry, nondimensionalized by the step height s.
    H_in = 5.2 / 4.9
    H_out = H_in + 1.0
    L_up = 200.0 / 4.9
    L_down = 100.0 / 4.9
    nx_up = 5 * n
    nx_down = max(2 * n, round(nx_up * L_down / L_up))
    nodes, key, cells = [], {}, []

    def nid(x, y):
        k = (round(x, 14), round(y, 14))
        if k not in key:
            key[k] = len(nodes)
            nodes.append((x, y))
        return key[k]

    def block(x0, x1, y0, y1, nx, ny):
        ids = [
            [
                nid(x0 + (x1 - x0) * i / nx, y0 + (y1 - y0) * j / ny)
                for i in range(nx + 1)
            ]
            for j in range(ny + 1)
        ]
        for j in range(ny):
            for i in range(nx):
                cells.append(
                    (ids[j][i], ids[j][i + 1], ids[j + 1][i + 1], ids[j + 1][i])
                )
        return ids

    up = block(-L_up, 0.0, 1.0, 1.0 + H_in, nx_up, n)
    lo = block(0.0, L_down, 0.0, 1.0, nx_down, n)
    hi = block(0.0, L_down, 1.0, H_out, nx_down, n)
    inlet = [(up[j + 1][0], up[j][0]) for j in range(n)]
    # The outlet spans both downstream blocks.
    outlet = [(lo[j][-1], lo[j + 1][-1]) for j in range(n)] + [
        (hi[j + 1][-1], hi[j][-1]) for j in range(n)
    ]
    wall = []
    wall += [(up[0][i], up[0][i + 1]) for i in range(nx_up)]
    wall += [(lo[0][i + 1], lo[0][i]) for i in range(nx_down)]
    wall += [(hi[-1][i], hi[-1][i + 1]) for i in range(nx_down)]
    # The upper wall of the upstream duct is physical wall.
    wall += [(up[-1][i], up[-1][i + 1]) for i in range(nx_up)]
    wall += [(lo[j][0], lo[j + 1][0]) for j in range(n)]
    return extrude_2d(
        nodes,
        cells,
        {"inlet": inlet, "outlet": outlet, "wall": wall, "_front_cells": cells},
    )


def naca0012_surface(n: int):
    t = 0.12
    beta = np.linspace(0.0, math.pi, n + 1)
    x = 0.5 * (1.0 - np.cos(beta))
    yt = (
        5.0
        * t
        * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x * x
            + 0.2843 * x**3
            - 0.1015 * x**4
        )
    )
    # Closed polygon uses the finite NACA 0012 trailing-edge thickness.
    upper = list(zip(x[::-1], yt[::-1]))  # TE upper -> LE
    lower = list(zip(x[1:], -yt[1:]))  # after LE -> TE lower
    return upper + lower


def naca_o_grid(
    n_surface: int,
    n_radial: int,
    radius: float = 20.0,
    n_spanwise: int = 1,
    span: float = 1.0,
    split_farfield: bool = False,
    radial_power: float = 1.8,
):
    """Build the deterministic transfinite-blend NACA0012 O-grid.

    Splitting farfield by the sign of the outer-face x-normal is an
    experimental A/B setup; the standard qualification mesh is unchanged.
    """
    if n_surface < 2:
        raise ValueError("n_surface must be at least 2")
    if n_radial < 1:
        raise ValueError("n_radial must be positive")
    if not math.isfinite(radius) or radius <= 0.0:
        raise ValueError("radius must be finite and positive")
    if not math.isfinite(radial_power) or not 1.0 <= radial_power <= 3.0:
        raise ValueError("radial_power must be finite in [1, 3]")

    surf = naca0012_surface(n_surface)
    m = len(surf)
    # Parameterize the outer boundary by the same closed-loop index. The
    # surface itself is not a mesh oracle; this is only a deterministic,
    # Gmsh-free body-fitted starting mesh.
    # Build a wall-normal first layer and blend it smoothly into the
    # circular farfield.  The surface ordering is clockwise, so (ty, -tx)
    # is the outward unit normal on both sides of the airfoil.
    normals = []
    for i in range(m):
        x0, y0 = surf[(i - 1) % m]
        x1, y1 = surf[(i + 1) % m]
        tx, ty = x1 - x0, y1 - y0
        length = math.hypot(tx, ty)
        if not math.isfinite(length) or length <= 0.0:
            raise ValueError("degenerate NACA0012 surface tangent")
        normals.append((ty / length, -tx / length))

    pts = []
    for j in range(n_radial + 1):
        eta = j / n_radial
        rfac = eta**radial_power
        normal_weight = (1.0 - eta) ** 10
        for i, (x, y) in enumerate(surf):
            theta = 2.0 * math.pi * (i + 0.5) / m
            xo, yo = radius * math.cos(theta), radius * math.sin(theta)
            dx, dy = xo - x, yo - y
            radial_length = math.hypot(dx, dy)
            nx, ny = normals[i]
            vx = (1.0 - normal_weight) * dx + normal_weight * radial_length * nx
            vy = (1.0 - normal_weight) * dy + normal_weight * radial_length * ny
            pts.append((x + rfac * vx, y + rfac * vy))

    def P(j, i):
        return j * m + (i % m)

    cells2d = []
    for j in range(n_radial):
        for i in range(m):
            cells2d.append((P(j, i), P(j, i + 1), P(j + 1, i + 1), P(j + 1, i)))

    airfoil = [(P(0, i + 1), P(0, i)) for i in range(m)]
    farfield = [(P(n_radial, i), P(n_radial, i + 1)) for i in range(m)]
    if split_farfield:
        inflow, outflow = [], []
        for edge in farfield:
            # On the circular outer boundary the x-coordinate of the face
            # midpoint has the same sign as the outward normal x-component.
            x_mid = 0.5 * (pts[edge[0]][0] + pts[edge[1]][0])
            (outflow if x_mid >= 0.0 else inflow).append(edge)
        if not inflow or not outflow or len(inflow) + len(outflow) != len(farfield):
            raise RuntimeError("NACA farfield split produced invalid patches")
        outer_patches = {"farfield_inlet": inflow, "farfield_outlet": outflow}
    else:
        outer_patches = {"farfield": farfield}
    # Periodic seams are interior topology because i is wrapped; the only
    # physical boundaries are the airfoil and farfield.
    nodes2d = pts
    return extrude_2d(
        nodes2d,
        cells2d,
        {
            "airfoil": airfoil,
            **outer_patches,
            "_front_cells": cells2d,
        },
        layers=n_spanwise,
        span=span,
    )


def naca_o_grid_quality(
    n_surface: int, n_radial: int, radius: float = 20.0, radial_power: float = 1.8
):
    """Return geometric quality indicators for the 2-D NACA0012 O-grid.

    These metrics are diagnostics, not a CFD qualification or a force oracle.
    Cell area is evaluated from each quadrilateral's polygon; radial spacing
    is measured between corresponding nodes on adjacent layers.
    """
    points, _cells, _ = naca_o_grid(
        n_surface, n_radial, radius, radial_power=radial_power
    )
    points = np.asarray(points, dtype=np.float64)
    m = 2 * n_surface + 1
    planar = points[: (n_radial + 1) * m, :2].reshape(n_radial + 1, m, 2)

    radial_steps = np.linalg.norm(np.diff(planar, axis=0), axis=2)
    first_layer = radial_steps[0]
    positive_steps = radial_steps[radial_steps > np.finfo(float).eps]
    if positive_steps.size != radial_steps.size:
        raise ValueError("NACA0012 O-grid contains a zero-length radial edge")

    growth = radial_steps[1:] / radial_steps[:-1]
    signed_areas = []
    aspect_ratios = []
    for j in range(n_radial):
        for i in range(m):
            quad = np.asarray(
                (
                    planar[j, i],
                    planar[j, (i + 1) % m],
                    planar[j + 1, (i + 1) % m],
                    planar[j + 1, i],
                )
            )
            x, y = quad[:, 0], quad[:, 1]
            signed_areas.append(0.5 * np.sum(x * np.roll(y, -1) - y * np.roll(x, -1)))
            edge_lengths = np.linalg.norm(np.roll(quad, -1, axis=0) - quad, axis=1)
            aspect_ratios.append(
                float(
                    np.max(edge_lengths)
                    / max(np.min(edge_lengths), np.finfo(float).tiny)
                )
            )

    # Measure orthogonality on both families of internal faces. A positive
    # cell area alone does not detect highly skewed faces that can destabilize
    # the pressure-correction operator.
    cell_centres = np.empty((n_radial, m, 2), dtype=np.float64)
    for j in range(n_radial):
        for i in range(m):
            cell_centres[j, i] = 0.25 * (
                planar[j, i]
                + planar[j, (i + 1) % m]
                + planar[j + 1, (i + 1) % m]
                + planar[j + 1, i]
            )

    nonorthogonality = []
    for j in range(1, n_radial):
        for i in range(m):
            edge = planar[j, (i + 1) % m] - planar[j, i]
            normal = np.array((edge[1], -edge[0]))
            delta = cell_centres[j, i] - cell_centres[j - 1, i]
            denom = float(np.linalg.norm(normal)) * float(np.linalg.norm(delta))
            cosine = abs(float(np.dot(normal, delta))) / denom
            nonorthogonality.append(
                math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
            )
    for j in range(n_radial):
        for i in range(m):
            edge = planar[j + 1, i] - planar[j, i]
            normal = np.array((edge[1], -edge[0]))
            delta = cell_centres[j, i] - cell_centres[j, (i - 1) % m]
            denom = float(np.linalg.norm(normal)) * float(np.linalg.norm(delta))
            cosine = abs(float(np.dot(normal, delta))) / denom
            nonorthogonality.append(
                math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
            )

    nonorthogonality = np.asarray(nonorthogonality)
    signed_areas = np.asarray(signed_areas)
    nonzero_areas = signed_areas[np.abs(signed_areas) > np.finfo(float).eps]
    dominant_sign = np.sign(np.median(nonzero_areas)) if nonzero_areas.size else 0.0
    inverted_or_degenerate = int(
        np.count_nonzero(
            (np.abs(signed_areas) <= np.finfo(float).eps)
            | (np.sign(signed_areas) != dominant_sign)
        )
    )
    return {
        "n_surface": int(n_surface),
        "n_radial": int(n_radial),
        "radius": float(radius),
        "radial_power": float(radial_power),
        "cell_count": int(n_surface * 2 * n_radial + n_radial),
        "inverted_or_degenerate_cells": inverted_or_degenerate,
        "minimum_absolute_cell_area": float(np.min(np.abs(signed_areas))),
        "minimum_first_layer_spacing": float(np.min(first_layer)),
        "mean_first_layer_spacing": float(np.mean(first_layer)),
        "maximum_first_layer_spacing": float(np.max(first_layer)),
        "minimum_radial_growth_ratio": float(np.min(growth)),
        "maximum_radial_growth_ratio": float(np.max(growth)),
        "maximum_cell_edge_aspect_ratio": float(np.max(aspect_ratios)),
        "maximum_nonorthogonality_degrees": float(np.max(nonorthogonality)),
        "p95_nonorthogonality_degrees": float(np.percentile(nonorthogonality, 95.0)),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("case", choices=("channel", "bfs", "naca0012", "naca3d", "all"))
    ap.add_argument("--output-dir", type=Path, default=Path("build/validation_meshes"))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument(
        "--naca-radial-power",
        type=float,
        default=1.8,
        help="experimental O-grid radial spacing power; baseline 1.8",
    )
    ap.add_argument(
        "--split-naca-farfield",
        action="store_true",
        help="experimental inlet/outlet outer patches; distinct mesh filenames",
    )
    args = ap.parse_args()

    out = args.output_dir
    if args.case in ("channel", "all"):
        for n in (16,) if args.quick else (16, 32, 64):
            p, c, b = channel(4 * n, n)
            write_mesh(out / f"channel_n{n}.h5", p, c, b)
    if args.case in ("bfs", "all"):
        levels = (16,) if args.quick else (16, 32, 64)
        for n in levels:
            p, c, b = bfs(n)
            write_mesh(out / f"bfs_re200_n{n}.h5", p, c, b)
    if args.case in ("naca0012", "all"):
        levels = ((64, 24),) if args.quick else ((64, 24), (128, 48), (256, 96))
        quality_reports = []
        for n, nr in levels:
            p, c, b = naca_o_grid(
                n,
                nr,
                20.0,
                split_farfield=args.split_naca_farfield,
                radial_power=args.naca_radial_power,
            )
            suffix = "_split" if args.split_naca_farfield else ""
            if args.naca_radial_power != 1.8:
                suffix += "_rp" + str(args.naca_radial_power).replace(".", "p")
            write_mesh(out / f"naca0012_laminar_n{n}{suffix}.h5", p, c, b)
            quality = naca_o_grid_quality(
                n, nr, 20.0, radial_power=args.naca_radial_power
            )
            quality_reports.append(quality)
            print(
                json.dumps(
                    {"mesh": f"naca0012_laminar_n{n}{suffix}.h5", **quality},
                    sort_keys=True,
                )
            )
        report_path = out / (
            "naca0012_mesh_quality_split.json"
            if args.split_naca_farfield
            else "naca0012_mesh_quality.json"
        )
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report_path.write_text(
            json.dumps(quality_reports, indent=2) + "\n", encoding="utf-8"
        )
    if args.case == "naca3d":
        # Two layers: 3-D topology/planar-extrusion smoke only. Not a
        # finite-span wing or a thermal/turbulence physical qualification.
        p, c, b = naca_o_grid(64, 24, 20.0, n_spanwise=2, span=1.0)
        write_mesh(out / "naca0012_extruded_3d_n64_z2.h5", p, c, b)
        p, c, b = naca_o_grid(16, 8, 20.0, n_spanwise=2, span=1.0)
        write_mesh(out / "naca0012_extruded_3d_n16_z2.h5", p, c, b)
        # Thermal through-span MMS needs physical end caps; the fluid
        # extrusion fixture intentionally retains EMPTY front/back patches.
        write_mesh(
            out / "naca0012_extruded_3d_thermal_n16_z2.h5",
            p,
            c,
            b,
            patch_type_overrides={"front": 0, "back": 0},
        )
        # A four-layer spanwise fixture is required to identify the complete
        # 3-D quadratic basis at wall faces (the two-layer smoke mesh is rank
        # deficient for this purpose).
        p, c, b = naca_o_grid(16, 8, 20.0, n_spanwise=4, span=1.0)
        write_mesh(out / "naca0012_extruded_3d_wall_n16_z4.h5", p, c, b)
    print(f"validation meshes written to {out}")


if __name__ == "__main__":
    main()
