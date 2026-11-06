"""Fail-closed topology/quality inspection of CFDX TH-03 blockMesh or polyMesh.

Independent of solver convergence. Face orientation, connectivity, 2D cell
areas, normal orthogonality, and geometric patch positions are checked.
These geometric gates alone do not establish SST wall resolution or accuracy.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import h5py
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from generate_th03_openfoam_blockmesh import (
    PATCH_TYPES,
    H,
    count_invalid_extrusion_patch_faces,
)


def _area_2d(points: np.ndarray) -> float:
    """Compute positive polygon area from convex hull of a hex prism's vertices."""
    p = sorted(set(map(tuple, points[:, :2])))
    if len(p) < 3:
        return 0.0

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for pt in p:
        while len(lower) > 1 and cross(lower[-2], lower[-1], pt) <= 0:
            lower.pop()
        lower.append(pt)
    upper = []
    for pt in p[::-1]:
        while len(upper) > 1 and cross(upper[-2], upper[-1], pt) <= 0:
            upper.pop()
        upper.append(pt)
    h = lower[:-1] + upper[:-1]
    return abs(
        sum(
            h[i][0] * h[(i + 1) % len(h)][1] - h[(i + 1) % len(h)][0] * h[i][1]
            for i in range(len(h))
        )
        * 0.5
    )


def audit(mesh: dict) -> dict:
    points = np.asarray(mesh["points"], dtype=float)
    fv = np.asarray(mesh["face_vertices"], dtype=np.int64)
    fo = np.asarray(mesh["face_offsets"], dtype=np.int64)
    owners = np.asarray(mesh["owner"], dtype=np.int64)
    neighbours = np.asarray(mesh["neighbour"], dtype=np.int64)
    cf = np.asarray(mesh["cell_faces"], dtype=np.int64)
    co = np.asarray(mesh["cell_offsets"], dtype=np.int64)
    nc = int(mesh["n_cells"])
    nf = int(mesh["n_faces"])
    issues = []
    if not (
        len(fo) == nf + 1
        and len(co) == nc + 1
        and len(owners) == nf
        and len(neighbours) == nf
    ):
        raise ValueError("invalid CSR dimensions")
    if not (fo[0] == 0 and fo[-1] == len(fv) and co[0] == 0 and co[-1] == len(cf)):
        raise ValueError("invalid CSR offsets")
    if np.any(np.diff(fo) < 3) or np.any(np.diff(co) < 4):
        raise ValueError("too few vertices on a face or faces in a cell")
    if not (np.all(np.isfinite(points)) and points.ndim == 2 and points.shape[1] == 3):
        raise ValueError("invalid non-finite mesh points")
    if (
        np.any(fv < 0)
        or np.any(fv >= len(points))
        or np.any(cf < 0)
        or np.any(cf >= nf)
    ):
        raise ValueError("invalid CSR vertex/face index")
    if (
        np.any(owners < 0)
        or np.any(owners >= nc)
        or np.any(neighbours >= nc)
        or np.any(neighbours < -1)
    ):
        raise ValueError("owner or neighbour outside cell range")
    if np.any(owners == neighbours):
        raise ValueError("internal face is attached to itself")

    # Collect cell unique point sets from native CSR topology, not generated block IDs.
    cell_centers = np.zeros((nc, 3), dtype=float)
    areas = np.zeros(nc, dtype=float)
    for c in range(nc):
        faces = cf[co[c] : co[c + 1]]
        ids = np.unique(np.concatenate([fv[fo[f] : fo[f + 1]] for f in faces]))
        pts = points[ids]
        cell_centers[c] = pts.mean(axis=0)
        areas[c] = _area_2d(pts)
        if not np.isfinite(areas[c]) or areas[c] <= 1e-18:
            issues.append(f"degenerate cell {c}, area {areas[c]}")
    centroids = np.zeros((nf, 3))
    normal = np.zeros((nf, 3))
    for f in range(nf):
        poly = points[fv[fo[f] : fo[f + 1]]]
        centroids[f] = poly.mean(axis=0)
        # Newell's polygon normal works also for skewed quads.
        normal[f] = 0.5 * np.sum(np.cross(poly, np.roll(poly, -1, axis=0)), axis=0)
    magnitudes = np.linalg.norm(normal, axis=1)
    if not np.all(np.isfinite(magnitudes)) or np.any(magnitudes <= 0):
        issues.append("zero/non-finite face area")
    interior = np.flatnonzero(neighbours >= 0)
    displacement = cell_centers[neighbours[interior]] - cell_centers[owners[interior]]
    dist = np.linalg.norm(displacement, axis=1)
    dot = np.sum(normal[interior] * displacement, axis=1)
    cosine = dot / np.maximum(magnitudes[interior] * dist, 1e-300)
    flipped = int(np.count_nonzero(cosine <= 0))
    if flipped:
        issues.append(f"{flipped} inward-oriented internal faces")
    angles = np.degrees(np.arccos(np.clip(cosine, -1, 1)))
    skew = np.linalg.norm(
        centroids[interior]
        - (
            cell_centers[owners[interior]]
            + displacement
            * np.sum(
                (centroids[interior] - cell_centers[owners[interior]]) * displacement,
                axis=1,
            )[:, None]
            / np.maximum(dist[:, None] ** 2, 1e-300)
        ),
        axis=1,
    )
    skew_rel = skew / np.maximum(dist, 1e-300)
    # Detect incorrect boundary normals independent of the case patch label.
    boundary = np.flatnonzero(neighbours < 0)
    outdot = np.sum(
        normal[boundary] * (centroids[boundary] - cell_centers[owners[boundary]]),
        axis=1,
    )
    bad_boundary = int(np.count_nonzero(outdot <= 0))
    if bad_boundary:
        issues.append(f"{bad_boundary} boundary normals point into owner")
    patch_faces = {p["name"]: p["face_ids"] for p in mesh["boundary_patch_assignment"]}
    if set(patch_faces) != set(PATCH_TYPES):
        issues.append("missing/unexpected patch names")
    invalid_extrusion_patch_faces = count_invalid_extrusion_patch_faces(mesh)
    if invalid_extrusion_patch_faces:
        issues.append(
            f"{invalid_extrusion_patch_faces} faces have incompatible extrusion patch geometry"
        )
    listed = [int(f) for fids in patch_faces.values() for f in fids]
    if (
        len(listed) != len(boundary)
        or len(set(listed)) != len(boundary)
        or set(listed) != set(map(int, boundary))
    ):
        issues.append("boundary patch face coverage is incomplete/duplicated")
    coordinates_min = points.min(axis=0) / H
    coordinates_max = points.max(axis=0) / H
    if not np.allclose(coordinates_min, [-130, 0, 0], atol=1e-8) or not np.allclose(
        coordinates_max, [50, 9, 1], atol=1e-8
    ):
        issues.append("unexpected domain extents in step heights")
    return {
        "status": "PASS_GEOMETRY_ONLY" if not issues else "FAIL_GEOMETRY",
        "qualification": "NOT_QUALIFIED",
        "issues": issues,
        "n_cells": nc,
        "n_faces": nf,
        "n_internal_faces": len(interior),
        "min_cell_area_m2": float(np.min(areas)),
        "max_cell_area_m2": float(np.max(areas)),
        "max_nonorthogonality_deg": float(np.max(angles)),
        "p99_nonorthogonality_deg": float(np.percentile(angles, 99)),
        "max_relative_face_skew": float(np.max(skew_rel)),
        "p99_relative_face_skew": float(np.percentile(skew_rel, 99)),
        "flipped_internal_faces": flipped,
        "flipped_boundary_faces": bad_boundary,
        "invalid_extrusion_patch_faces": invalid_extrusion_patch_faces,
        "extent_over_step_height": {
            "min": coordinates_min.tolist(),
            "max": coordinates_max.tolist(),
        },
        "boundary_face_counts": {k: len(v) for k, v in patch_faces.items()},
    }


def read_hdf5(path: Path) -> dict:
    with h5py.File(path) as h:
        keys = (
            "points",
            "face_vertices",
            "face_offsets",
            "owner",
            "neighbour",
            "cell_faces",
            "cell_offsets",
        )
        mesh = {k: h[k][:] for k in keys}
        metadata = h.attrs["boundary_patches"]
        raw = str(metadata)
        patch = {}
        for descriptor in raw.split(";"):
            if not descriptor.strip():
                continue
            name, kind, _start, _size = descriptor.split(":")
            patch[name] = {"name": name, "type": kind, "face_ids": []}
        raw_counts = np.asarray(h["patch_face_offsets"], dtype=np.int64)
        face_ids = np.asarray(h["patch_face_ids"], dtype=np.int64)
        for entry, begin, end in zip(patch.values(), raw_counts[:-1], raw_counts[1:]):
            entry["face_ids"] = face_ids[int(begin) : int(end)].tolist()
        mesh.update(
            n_cells=int(h.attrs["n_cells"]),
            n_faces=int(h.attrs["n_faces"]),
            boundary_patch_assignment=list(patch.values()),
        )
    return mesh


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("case", type=Path)
    p.add_argument("--output", type=Path)
    a = p.parse_args()
    result = audit(read_hdf5(a.case))
    output = json.dumps(result, indent=2)
    print(output)
    if a.output:
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(output + "\n")
    return 0 if result["status"] == "PASS_GEOMETRY_ONLY" else 2


if __name__ == "__main__":
    raise SystemExit(main())
