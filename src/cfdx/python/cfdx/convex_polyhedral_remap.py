"""First-order volume-conservative remapping between verified *convex* CFDX cells.

For small-to-medium diagnostic meshes only. This is not a production general
polyhedron intersector: non-convex cells, warped faces, ambiguous topology,
non-matching domains and pair-budget overruns fail explicitly. Input cell
averages are intensive quantities. SciPy's Qhull validates the geometry;
intersection vertices are calculated from the two sets of supporting planes.
The transfer is first-order accurate, not the second-order supermesh method.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

try:
    from scipy.spatial import ConvexHull, QhullError
except (
    ImportError
):  # Optional scientific geometry dependency; not needed for standard CLI.
    ConvexHull = None
    QhullError = RuntimeError

from .conservative_remap import ConservativeRemapReport, _unit_and_time
from .progressive_fields import _field
from .simulation_compare import IncompatibleSimulations

_TOPOLOGY = ("points", "face_offsets", "face_vertices", "cell_offsets", "cell_faces")
_MAX_CELLS = 256
_MAX_PAIRS = 20_000
_MAX_BYTES = 48 * 1024 * 1024
_MAX_PLANES = 24  # safety for an exhaustive triple-plane intersection


@dataclass(frozen=True, slots=True)
class _Cell:
    vertices: np.ndarray
    planes: np.ndarray
    lower: np.ndarray
    upper: np.ndarray
    volume: float


def _safe_data(h5: h5py.File, key: str, rank: int, kind: str) -> np.ndarray:
    link = h5.get(key, getlink=True)
    ds = h5.get(key)
    if (
        not isinstance(link, h5py.HardLink)
        or not isinstance(ds, h5py.Dataset)
        or ds.is_virtual
        or ds.external
        or ds.ndim != rank
        or ds.dtype.kind not in kind
        or ds.dtype.itemsize > 8
        or ds.size * ds.dtype.itemsize > _MAX_BYTES
    ):
        raise IncompatibleSimulations(
            f"invalid or untrusted convex mesh dataset: {key}"
        )
    if ds.chunks and math.prod(ds.chunks) * ds.dtype.itemsize > 8 * 1024 * 1024:
        raise IncompatibleSimulations("convex mesh compressed HDF5 chunk exceeds 8 MiB")
    return np.asarray(ds[:], dtype=np.float64 if key == "points" else np.int64)


def _read_cells(h5: h5py.File) -> list[_Cell]:
    if ConvexHull is None:
        raise IncompatibleSimulations(
            "convex polyhedral remapping requires scipy (install cfdx[remap])"
        )
    expected_bytes = 0
    for key in _TOPOLOGY:
        ds = h5.get(key)
        if not isinstance(ds, h5py.Dataset):
            raise IncompatibleSimulations(f"missing convex mesh dataset: {key}")
        expected_bytes += int(ds.size * ds.dtype.itemsize)
    if expected_bytes > _MAX_BYTES:
        raise IncompatibleSimulations("convex geometry exceeds combined memory budget")
    raw = {
        k: _safe_data(h5, k, 2 if k == "points" else 1, "f" if k == "points" else "iu")
        for k in _TOPOLOGY
    }
    p, fo, fv, co, cf = (raw[k] for k in _TOPOLOGY)
    n = len(co) - 1
    if (
        p.ndim != 2
        or p.shape[1] != 3
        or len(p) == 0
        or not np.isfinite(p).all()
        or n < 1
        or n > _MAX_CELLS
        or len(fo) < 2
        or fo[0] != 0
        or fo[-1] != len(fv)
        or co[0] != 0
        or co[-1] != len(cf)
        or (np.diff(fo) < 3).any()
        or (np.diff(co) < 4).any()
        or (fv < 0).any()
        or (fv >= len(p)).any()
        or (cf < 0).any()
        or (cf >= len(fo) - 1).any()
    ):
        raise IncompatibleSimulations(
            "invalid convex polyhedral topology or coordinate array"
        )
    cells = []
    for i in range(n):
        ids = cf[co[i] : co[i + 1]]
        if len(ids) > _MAX_PLANES or len(set(map(int, ids))) != len(ids):
            raise IncompatibleSimulations("too many or repeated convex cell faces")
        faces = [np.asarray(fv[fo[f] : fo[f + 1]], dtype=np.int64) for f in ids]
        if any(
            len(face) > 16 or len(set(map(int, face))) != len(face) for face in faces
        ):
            raise IncompatibleSimulations("convex cell has malformed polygonal faces")
        vertex_ids = np.unique(np.concatenate(faces))
        xyz = p[vertex_ids]
        if len(xyz) < 4:
            raise IncompatibleSimulations("convex cell has fewer than four vertices")
        size = float(np.max(np.ptp(xyz, axis=0)))
        if size <= 0 or not math.isfinite(size):
            raise IncompatibleSimulations("zero-size convex cell")
        coordinate_scale = float(np.max(np.abs(xyz)))
        if 32 * np.finfo(float).eps * coordinate_scale > size * 1e-8:
            raise IncompatibleSimulations(
                "convex geometry coordinate precision insufficient"
            )
        try:
            hull = ConvexHull(xyz)
        except QhullError as exc:
            raise IncompatibleSimulations("degenerate or coplanar convex cell") from exc
        if hull.volume <= (size**3) * 1e-13 or len(hull.vertices) != len(xyz):
            raise IncompatibleSimulations("zero-volume or nonconvex cell vertices")
        center = np.mean(xyz, axis=0)
        eps = max(1e-13, size * 1e-9)
        directed_edges: dict[tuple[int, int], list[int]] = {}
        signed_volume = 0.0
        for face in faces:
            coords = p[face]
            vec = np.zeros(3, dtype=np.float64)
            for k in range(1, len(coords) - 1):
                vec += np.cross(coords[k] - coords[0], coords[k + 1] - coords[0])
            norm = float(np.linalg.norm(vec))
            if norm <= eps * size:
                raise IncompatibleSimulations("degenerate convex cell face")
            outward = vec / norm
            if float(np.dot(outward, np.mean(coords, axis=0) - center)) < 0:
                face = face[::-1]
                coords = coords[::-1]
                outward = -outward
            if max(abs(float(np.dot(outward, v - coords[0]))) for v in coords) > eps:
                raise IncompatibleSimulations("nonplanar convex cell face")
            if np.max((xyz - coords[0]) @ outward) > eps:
                raise IncompatibleSimulations("concave polyhedron or inward face")
            # Face polygon must be convex; a concave polygon would cover a
            # different portion of the hull and invalidate signed volumes.
            for k in range(len(coords)):
                turn = np.dot(
                    np.cross(
                        coords[(k + 1) % len(coords)] - coords[k],
                        coords[(k + 2) % len(coords)] - coords[(k + 1) % len(coords)],
                    ),
                    outward,
                )
                if turn < -(eps * size):
                    raise IncompatibleSimulations(
                        "nonconvex or self-intersecting polygonal face"
                    )
            for k in range(len(face)):
                a, b = int(face[k]), int(face[(k + 1) % len(face)])
                key = (min(a, b), max(a, b))
                directed_edges.setdefault(key, []).append(1 if a < b else -1)
            for k in range(1, len(coords) - 1):
                signed_volume += (
                    float(
                        np.dot(
                            coords[0] - center,
                            np.cross(coords[k] - center, coords[k + 1] - center),
                        )
                    )
                    / 6.0
                )
        if (
            not directed_edges
            or any(sorted(dirs) != [-1, 1] for dirs in directed_edges.values())
            or len(vertex_ids) - len(directed_edges) + len(faces) != 2
        ):
            raise IncompatibleSimulations("open or nonmanifold convex cell surface")
        if not math.isclose(
            signed_volume, float(hull.volume), rel_tol=1e-8, abs_tol=size**3 * 1e-11
        ):
            raise IncompatibleSimulations(
                "polyhedral face volume disagrees with convex hull"
            )
        cells.append(
            _Cell(
                xyz,
                hull.equations.astype(np.float64),
                xyz.min(axis=0),
                xyz.max(axis=0),
                float(hull.volume),
            )
        )
    return cells


def _intersection_volume(a: _Cell, b: _Cell) -> float:
    if np.any(np.minimum(a.upper, b.upper) <= np.maximum(a.lower, b.lower)):
        return 0.0
    planes = np.vstack((a.planes, b.planes))
    m, n = len(a.planes), len(b.planes)
    if m + n > 2 * _MAX_PLANES:
        raise IncompatibleSimulations("convex intersection exceeds plane budget")
    size = max(float(np.max(a.upper - a.lower)), float(np.max(b.upper - b.lower)))
    magnitude = max(
        float(np.max(np.abs(a.vertices))), float(np.max(np.abs(b.vertices)))
    )
    tol = max(1e-13, size * 1e-9, 32 * np.finfo(float).eps * magnitude)
    if tol >= size * 1e-6:
        raise IncompatibleSimulations(
            "convex intersection coordinate precision insufficient"
        )
    candidates: list[np.ndarray] = []
    for vertices, other in ((a.vertices, b.planes), (b.vertices, a.planes)):
        for point in vertices:
            if np.all(other[:, :3] @ point + other[:, 3] <= tol):
                candidates.append(point)
    # Triple-plane vertices. Triple intersections within one hull are already
    # present among source/target vertices, so only mixed triples are needed.
    for triplet in itertools.combinations(range(m + n), 3):
        if all(k < m for k in triplet) or all(k >= m for k in triplet):
            continue
        normals = planes[list(triplet), :3]
        determinant = float(np.linalg.det(normals))
        if abs(determinant) < 1e-10:
            continue
        point = np.linalg.solve(normals, -planes[list(triplet), 3])
        if np.isfinite(point).all() and np.all(
            planes[:, :3] @ point + planes[:, 3] <= tol
        ):
            candidates.append(point)
    if len(candidates) < 4:
        return 0.0
    # Quantise *relative* to the local intersection bounding box. Quantising
    # global coordinates would overflow int64 for far-from-origin meshes.
    origin = np.minimum(a.lower, b.lower)
    grid = np.round((np.asarray(candidates) - origin) / tol)
    if not np.isfinite(grid).all() or np.max(np.abs(grid)) > 9e15:
        raise IncompatibleSimulations(
            "convex intersection quantization budget exceeded"
        )
    keys = np.unique(grid.astype(np.int64), axis=0)
    if len(keys) < 4:
        return 0.0
    coords = origin + keys.astype(np.float64) * tol
    try:
        hull = ConvexHull(coords)
    except QhullError:
        return 0.0  # A shared edge/face has zero intersection volume.
    result = float(hull.volume)
    if not math.isfinite(result) or result < 0:
        raise IncompatibleSimulations("non-finite convex intersection volume")
    if result <= size**3 * 1e-12:
        return 0.0
    if result > min(a.volume, b.volume) * (1 + 1e-7):
        raise IncompatibleSimulations("intersection volume exceeds input cell volume")
    return result


def conservative_convex_remap(
    source: str | Path, target: str | Path, field: str
) -> tuple[np.ndarray, ConservativeRemapReport]:
    """Transfer intensive cell scalar averages on complete convex partitions.

    Validation requires full two-way donor/target coverage, finite data,
    matching time and units, positive cell volumes and conservative integrals.
    It does not certify the CFD solutions or handle concave general polyhedra.
    """
    source, target = Path(source), Path(target)
    if source.is_symlink() or target.is_symlink():
        raise IncompatibleSimulations("convex remap refuses symbolic link files")

    def signature(path: Path):
        st = path.stat()
        if not path.is_file():
            raise IncompatibleSimulations("convex remap requires regular files")
        return st.st_dev, st.st_ino, st.st_size, st.st_mtime_ns

    before = (signature(source), signature(target))
    with h5py.File(source, "r") as hs, h5py.File(target, "r") as ht:
        donors, receivers = _read_cells(hs), _read_cells(ht)
        if len(donors) * len(receivers) > _MAX_PAIRS:
            raise IncompatibleSimulations("convex remap pair budget exceeded")
        a, b = _field(hs, field), _field(ht, field)
        if a.size != len(donors) or b.size != len(receivers):
            raise IncompatibleSimulations(
                "convex remap scalar field/cell count mismatch"
            )
        units, time = _unit_and_time(hs, ht, a, b)
        va, vb = np.asarray(a[:], dtype=float), np.asarray(b[:], dtype=float)
        if not np.isfinite(va).all() or not np.isfinite(vb).all():
            raise IncompatibleSimulations("nonfinite convex remap scalar field")
        av = np.asarray([c.volume for c in donors])
        bv = np.asarray([c.volume for c in receivers])
        if not math.isclose(float(av.sum()), float(bv.sum()), rel_tol=2e-8):
            raise IncompatibleSimulations(
                "convex partitions have unequal domain volume"
            )
        mapped = np.empty(len(receivers))
        used = np.zeros(len(donors))
        for j, target_cell in enumerate(receivers):
            overlaps = np.asarray(
                [_intersection_volume(d, target_cell) for d in donors]
            )
            covered = math.fsum(map(float, overlaps))
            if not math.isclose(
                covered, target_cell.volume, rel_tol=2e-7, abs_tol=1e-12
            ):
                raise IncompatibleSimulations(
                    "convex target cell not fully covered by donor mesh"
                )
            used += overlaps
            mapped[j] = (
                math.fsum(float(x * y) for x, y in zip(overlaps, va))
                / target_cell.volume
            )
        if not np.allclose(used, av, rtol=2e-7, atol=1e-12):
            raise IncompatibleSimulations(
                "convex donor cells overlap or target partition incomplete"
            )
        source_int = math.fsum(float(x * y) for x, y in zip(av, va))
        mapped_int = math.fsum(float(x * y) for x, y in zip(bv, mapped))
        target_int = math.fsum(float(x * y) for x, y in zip(bv, vb))
        denom = max(math.fsum(float(x * abs(y)) for x, y in zip(av, va)), 1e-300)
        err = abs(source_int - mapped_int) / denom
        if not math.isfinite(err) or err > 1e-7:
            raise IncompatibleSimulations("convex remap integral not conserved")
        difference = mapped - vb
        if before != (signature(source), signature(target)):
            raise IncompatibleSimulations("convex remap HDF5 changed during evaluation")
        return mapped, ConservativeRemapReport(
            field,
            units,
            len(donors),
            len(receivers),
            "convex_polyhedron_overlap_first_order_qhull",
            time,
            source_int,
            target_int,
            mapped_int,
            err,
            float(np.linalg.norm(difference)),
            float(np.max(np.abs(difference))),
        )


def main(argv: list[str] | None = None) -> int:
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="CFDX convex-polyhedron conservative scalar remapping"
    )
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("field")
    args = parser.parse_args(argv)
    try:
        _, report = conservative_convex_remap(args.source, args.target, args.field)
    except (ValueError, OSError, IncompatibleSimulations) as exc:
        parser.exit(2, f"Convex remapping refused: {exc}\n")
    print(json.dumps(report.to_dict(), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
