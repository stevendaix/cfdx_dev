"""Exact cell-volume overlap remapping on verified axis-aligned hexahedra.

This is an explicitly limited conservative method, not a general polyhedral
intersector. Only closed, 6-face, 8-corner Cartesian boxes are accepted.
Unknown, skewed, non-closed and incompatible grids fail rather than falling
back to nearest-neighbour interpolation. Intensive cell-average scalars only.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from pathlib import Path

import h5py
import numpy as np

from .progressive_fields import _field
from .simulation_compare import IncompatibleSimulations

_MAX_CELLS = 4096
_MAX_PAIRS = 2_000_000
_MAX_INPUT_BYTES = 64 * 1024 * 1024
_TOPOLOGY = ("points", "face_offsets", "face_vertices", "cell_offsets", "cell_faces")


@dataclass(frozen=True)
class ConservativeRemapReport:
    field: str
    units: str
    source_cells: int
    target_cells: int
    method: str
    physical_time: float | None
    source_integral: float
    target_integral: float
    mapped_integral: float
    relative_conservation_error: float
    l2_difference: float
    linf_difference: float
    certificate: bool = False

    def to_dict(self) -> dict:
        return asdict(self)


def _safe_geometry(h5: h5py.File) -> tuple[np.ndarray, np.ndarray]:
    datasets: dict[str, h5py.Dataset] = {}
    byte_total = 0
    for key in _TOPOLOGY:
        if not isinstance(h5.get(key, getlink=True), h5py.HardLink):
            raise IncompatibleSimulations(f"untrusted geometry: {key}")
        dset = h5.get(key)
        if not isinstance(dset, h5py.Dataset) or dset.is_virtual or dset.external:
            raise IncompatibleSimulations(f"untrusted geometry dataset: {key}")
        if (
            dset.dtype.kind not in ("f" if key == "points" else "iu")
            or dset.dtype.itemsize > 8
        ):
            raise IncompatibleSimulations(f"unsupported geometry type: {key}")
        if (
            dset.chunks is not None
            and math.prod(dset.chunks) * dset.dtype.itemsize > 8 * 1024 * 1024
        ):
            raise IncompatibleSimulations("geometry HDF5 chunk exceeds memory budget")
        byte_total += dset.size * dset.dtype.itemsize
        datasets[key] = dset
    if byte_total > _MAX_INPUT_BYTES:
        raise IncompatibleSimulations(
            "geometry exceeds conservative remap memory budget"
        )
    points = np.asarray(datasets["points"][:], dtype=np.float64)
    if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
        raise IncompatibleSimulations("invalid 3D mesh points")
    fo = np.asarray(datasets["face_offsets"][:], dtype=np.int64)
    fv = np.asarray(datasets["face_vertices"][:], dtype=np.int64)
    co = np.asarray(datasets["cell_offsets"][:], dtype=np.int64)
    cf = np.asarray(datasets["cell_faces"][:], dtype=np.int64)
    n_cells = len(co) - 1
    if n_cells < 1 or n_cells > _MAX_CELLS:
        raise IncompatibleSimulations(
            "conservative remapping limited to 4096 cells per mesh"
        )
    if (
        len(fo) < 2
        or fo[0] != 0
        or fo[-1] != len(fv)
        or co[0] != 0
        or co[-1] != len(cf)
        or (np.diff(fo) < 0).any()
        or (np.diff(co) < 0).any()
        or (fv < 0).any()
        or (fv >= len(points)).any()
        or (cf < 0).any()
        or (cf >= len(fo) - 1).any()
    ):
        raise IncompatibleSimulations("broken cell/face connectivity")
    low = np.empty((n_cells, 3))
    high = np.empty((n_cells, 3))
    for cell in range(n_cells):
        faces = cf[co[cell] : co[cell + 1]]
        if len(faces) != 6 or len(set(map(int, faces))) != 6:
            raise IncompatibleSimulations(
                "remapping requires six distinct quadrilateral faces"
            )
        verts = []
        for fid in faces:
            ids = fv[fo[fid] : fo[fid + 1]]
            if len(ids) != 4 or len(set(map(int, ids))) != 4:
                raise IncompatibleSimulations(
                    "cell face must have four distinct vertices"
                )
            verts.extend(map(int, ids))
        unique = sorted(set(verts))
        if len(unique) != 8:
            raise IncompatibleSimulations("cell must have exactly eight corners")
        coordinates = points[unique]
        lo, hi = coordinates.min(axis=0), coordinates.max(axis=0)
        if not np.all(hi > lo):
            raise IncompatibleSimulations("zero-volume or inverted Cartesian cell")
        size = hi - lo
        eps = max(1e-12, float(np.max(size)) * 1e-10)
        corners = set()
        for point in coordinates:
            bits = []
            for axis in range(3):
                if abs(point[axis] - lo[axis]) <= eps:
                    bits.append(0)
                elif abs(point[axis] - hi[axis]) <= eps:
                    bits.append(1)
                else:
                    raise IncompatibleSimulations(
                        "non-Cartesian hexahedron; no conservative map"
                    )
            corners.add(tuple(bits))
        if len(corners) != 8:
            raise IncompatibleSimulations("duplicate hexahedron corner")
        checked_faces = set()
        for fid in faces:
            face_points = points[fv[fo[fid] : fo[fid + 1]]]
            axes = [
                (axis, side)
                for axis in range(3)
                for side, value in ((0, lo[axis]), (1, hi[axis]))
                if np.all(np.abs(face_points[:, axis] - value) <= eps)
            ]
            if len(axes) != 1:
                raise IncompatibleSimulations("invalid Cartesian face plane")
            checked_faces.add(axes[0])
        if len(checked_faces) != 6:
            raise IncompatibleSimulations("missing or duplicate Cartesian face")
        low[cell], high[cell] = lo, hi
    return low, high


def _unit_and_time(a: h5py.File, b: h5py.File, da, db) -> tuple[str, float | None]:
    u1, u2 = da.attrs.get("units"), db.attrs.get("units")
    if isinstance(u1, bytes):
        u1 = u1.decode("utf-8")
    if isinstance(u2, bytes):
        u2 = u2.decode("utf-8")
    if not isinstance(u1, str) or not u1 or u1 != u2:
        raise IncompatibleSimulations("field units missing or incompatible")
    times = [h.attrs.get("physical_time") for h in (a, b)]
    if any(
        v is not None and (isinstance(v, (str, bytes, bool)) or not np.isfinite(v))
        for v in times
    ):
        raise IncompatibleSimulations("physical time must be finite")
    if (times[0] is None) != (times[1] is None) or (
        times[0] is not None and float(times[0]) != float(times[1])
    ):
        raise IncompatibleSimulations("physical times incompatible")
    return u1, float(times[0]) if times[0] is not None else None


def conservative_box_remap(
    source: str | Path, target: str | Path, field: str
) -> tuple[np.ndarray, ConservativeRemapReport]:
    """Volume-weighted exact overlap; target is never overwritten.

    Returns target cell averages mapped from source, and comparison against
    the target's existing field (same name/units/time). Both meshes must cover
    the same domain and be non-overlapping Cartesian box partitions.
    """
    source, target = Path(source), Path(target)
    if source.is_symlink() or target.is_symlink():
        raise IncompatibleSimulations("conservative remap refuses symbolic-link files")

    def signature(path):
        stat = path.stat()
        if not path.is_file():
            raise IncompatibleSimulations(
                "conservative remap inputs must be regular files"
            )
        return stat.st_dev, stat.st_ino, stat.st_size, stat.st_mtime_ns

    before = signature(source), signature(target)
    with h5py.File(source, "r") as a, h5py.File(target, "r") as b:
        a_lo, a_hi = _safe_geometry(a)
        b_lo, b_hi = _safe_geometry(b)
        n, m = len(a_lo), len(b_lo)
        if n * m > _MAX_PAIRS:
            raise IncompatibleSimulations(
                "conservative remapping pair-count budget exceeded"
            )
        fa, fb = _field(a, field), _field(b, field)
        if fa.size != n or fb.size != m:
            raise IncompatibleSimulations("scalar size does not match mesh cells")
        units, time = _unit_and_time(a, b, fa, fb)
        va = np.asarray(fa[:], dtype=np.float64)
        vb = np.asarray(fb[:], dtype=np.float64)
        if not np.isfinite(va).all() or not np.isfinite(vb).all():
            raise IncompatibleSimulations(
                "nonfinite field cannot be conservatively remapped"
            )
        if not (
            np.allclose(a_lo.min(axis=0), b_lo.min(axis=0), rtol=0, atol=1e-10)
            and np.allclose(a_hi.max(axis=0), b_hi.max(axis=0), rtol=0, atol=1e-10)
        ):
            raise IncompatibleSimulations("meshes have different bounding domains")
        vol_a = np.prod(a_hi - a_lo, axis=1)
        vol_b = np.prod(b_hi - b_lo, axis=1)
        if not np.isclose(vol_a.sum(), vol_b.sum(), rtol=1e-10, atol=0):
            raise IncompatibleSimulations("source and target domain volumes differ")
        mapped = np.empty(m, dtype=np.float64)
        used_a = np.zeros(n, dtype=np.float64)
        for j in range(m):
            overlap = np.prod(
                np.maximum(0.0, np.minimum(a_hi, b_hi[j]) - np.maximum(a_lo, b_lo[j])),
                axis=1,
            )
            covered = math.fsum(float(x) for x in overlap)
            if not math.isclose(covered, float(vol_b[j]), rel_tol=1e-10, abs_tol=0):
                raise IncompatibleSimulations(
                    "target cell is not fully covered by source boxes"
                )
            used_a += overlap
            mapped[j] = math.fsum(float(x * y) for x, y in zip(overlap, va)) / vol_b[j]
        if not np.allclose(used_a, vol_a, rtol=1e-10, atol=0):
            raise IncompatibleSimulations(
                "source cells are not fully covered or overlap"
            )
        source_integral = math.fsum(float(x * y) for x, y in zip(vol_a, va))
        mapped_integral = math.fsum(float(x * y) for x, y in zip(vol_b, mapped))
        target_integral = math.fsum(float(x * y) for x, y in zip(vol_b, vb))
        scale = math.fsum(float(x * abs(y)) for x, y in zip(vol_a, va))
        conservation = abs(source_integral - mapped_integral) / max(scale, 1e-300)
        if not math.isfinite(conservation) or conservation > 1e-10:
            raise IncompatibleSimulations("volume integral not conserved by the remap")
        difference = mapped - vb
        report = ConservativeRemapReport(
            field,
            units,
            n,
            m,
            "exact_axis_aligned_hexahedron_volume_overlap",
            time,
            source_integral,
            target_integral,
            mapped_integral,
            conservation,
            float(np.linalg.norm(difference)),
            float(np.max(np.abs(difference))),
        )
        if before != (signature(source), signature(target)):
            raise IncompatibleSimulations("input HDF5 file changed during remapping")
        return mapped, report


def main(argv: list[str] | None = None) -> int:
    """Headless, fail-closed conservative diagnostic for CI and batch users."""
    import argparse
    import json

    parser = argparse.ArgumentParser(
        description="Exact conservative Cartesian-box field remapping"
    )
    parser.add_argument("source", type=Path)
    parser.add_argument("target", type=Path)
    parser.add_argument("field", help="absolute /fields/... scalar path")
    args = parser.parse_args(argv)
    try:
        _, report = conservative_box_remap(args.source, args.target, args.field)
    except (IncompatibleSimulations, OSError, ValueError, KeyError) as exc:
        parser.exit(2, f"Remapping refused: {exc}\n")
    print(json.dumps(report.to_dict(), indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
