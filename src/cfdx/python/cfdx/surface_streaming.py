"""Bounded, sampled surface preview for CFDX HDF5 mesh metadata.

No volume/cell topology is materialized. This is a geometric LOD preview,
NOT a full-volume streaming VTK pipeline nor a conservative sampled field.
The sampled-face fraction and original face IDs remain explicit.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np


@dataclass(frozen=True)
class SurfacePreview:
    points: np.ndarray
    faces: np.ndarray
    patch_indices: np.ndarray
    original_face_ids: np.ndarray
    n_boundary_faces: int
    n_total_mesh_points: int
    source: str

    @property
    def coverage(self) -> float:
        return float(len(self.original_face_ids)) / self.n_boundary_faces


def _safe_dataset(h5, name, rank, kind, max_chunk_bytes=4 * 1024 * 1024):
    if not isinstance(h5.get(name, getlink=True), h5py.HardLink):
        raise ValueError(f"missing or untrusted HDF5 dataset {name}")
    ds = h5.get(name)
    if (
        not isinstance(ds, h5py.Dataset)
        or ds.is_virtual
        or ds.external
        or ds.ndim != rank
        or ds.dtype.kind not in kind
        or ds.dtype.itemsize > 8
    ):
        raise ValueError(f"unsupported or unsafe HDF5 dataset {name}")
    if ds.chunks and math.prod(ds.chunks) * ds.dtype.itemsize > max_chunk_bytes:
        raise ValueError(f"HDF5 chunk exceeds safe preview memory budget: {name}")
    return ds


def read_surface_preview(
    path: str | Path, *, max_faces: int = 5000, max_vertices: int = 30000
) -> SurfacePreview:
    """Read an evenly decimated subset of real boundary faces, not generated cells.

    Face/point indices are checked against the original datasets and all reads
    are bounded by max_faces/max_vertices. Never use the returned face order
    as a CFDX cell index; only original_face_ids are stable.
    """
    if not 1 <= max_faces <= 20000 or not 3 <= max_vertices <= 100000:
        raise ValueError("unsafe surface preview budget")
    path = Path(path)
    if path.is_symlink():
        raise ValueError("surface preview refuses a symbolic-link input")
    with h5py.File(path, "r") as h5:
        pts = _safe_dataset(h5, "points", 2, "fiu")
        offsets = _safe_dataset(h5, "face_offsets", 1, "iu")
        vertices = _safe_dataset(h5, "face_vertices", 1, "iu")
        bfaces = _safe_dataset(h5, "patch_face_ids", 1, "iu")
        po = _safe_dataset(h5, "patch_face_offsets", 1, "iu")
        if pts.shape[1] != 3 or pts.shape[0] == 0 or po.size < 2 or po.size > 4097:
            raise ValueError("invalid point/patch metadata")
        patch_offsets = np.asarray(po[:], dtype=np.int64)
        if (
            patch_offsets[0] != 0
            or patch_offsets[-1] != bfaces.size
            or np.any(np.diff(patch_offsets) < 0)
            or offsets.size < 2
        ):
            raise ValueError("broken patch/face offsets")
        count = int(bfaces.size)
        if count < 1:
            raise ValueError("no boundary faces available")
        indices = np.unique(
            np.linspace(0, count - 1, num=min(count, max_faces), dtype=np.int64)
        )
        patch_ix = np.searchsorted(patch_offsets, indices, side="right") - 1
        if np.any(patch_ix < 0) or np.any(patch_ix >= len(patch_offsets) - 1):
            raise ValueError("boundary face lies outside patch table")
        ids = np.asarray(bfaces[indices], dtype=np.int64)
        if np.any(ids < 0) or np.any(ids >= offsets.size - 1):
            raise ValueError("boundary face references an invalid face")
        face_ids_unique, inv = np.unique(ids, return_inverse=True)
        starts = np.asarray(offsets[face_ids_unique], dtype=np.int64)
        stops = np.asarray(offsets[face_ids_unique + 1], dtype=np.int64)
        sizes = stops - starts
        if (
            np.any(sizes < 3)
            or np.any(sizes > 16)
            or np.any(starts < 0)
            or np.any(stops > vertices.size)
        ):
            raise ValueError("boundary polygon has unsafe or invalid connectivity")
        if int(sizes[inv].sum()) > max_vertices:
            raise ValueError("surface face budget exceeds vertex-index budget")
        polys = [
            np.asarray(vertices[a:b], dtype=np.int64) for a, b in zip(starts, stops)
        ]
        polygons = [polys[index] for index in inv]
        all_ids = np.concatenate(polygons)
        unique_ids, map_back = np.unique(all_ids, return_inverse=True)
        if (
            len(unique_ids) > max_vertices
            or np.any(unique_ids < 0)
            or np.any(unique_ids >= pts.shape[0])
        ):
            raise ValueError("surface preview point budget or indices invalid")
        point_coordinates = np.asarray(pts[unique_ids], dtype=np.float64)
        if (
            point_coordinates.shape != (len(unique_ids), 3)
            or not np.isfinite(point_coordinates).all()
        ):
            raise ValueError("surface preview contains nonfinite coordinates")
        face_buffer = np.empty(len(map_back) + len(polygons), dtype=np.int64)
        pos = 0
        cursor = 0
        for polygon in polygons:
            length = len(polygon)
            face_buffer[pos] = length
            face_buffer[pos + 1 : pos + 1 + length] = map_back[cursor : cursor + length]
            pos += length + 1
            cursor += length
        return SurfacePreview(
            point_coordinates,
            face_buffer,
            patch_ix,
            ids,
            count,
            int(pts.shape[0]),
            str(path.resolve()),
        )


def main(argv: list[str] | None = None) -> int:
    """Display counts of verified sampled boundary faces, not volume fields."""
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="Bounded CFDX HDF5 boundary surface preview"
    )
    parser.add_argument("file", type=Path)
    parser.add_argument("--max-faces", type=int, default=5000)
    args = parser.parse_args(argv)
    try:
        result = read_surface_preview(args.file, max_faces=args.max_faces)
    except (OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"Surface preview refused: {exc}\n")
    print(
        json.dumps(
            {
                "source": result.source,
                "sampled_boundary_faces": len(result.original_face_ids),
                "total_boundary_faces": result.n_boundary_faces,
                "sample_fraction": result.coverage,
                "sampled_points": len(result.points),
                "full_volume_available": False,
                "convergence_certificate": False,
            },
            indent=2,
            allow_nan=False,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
