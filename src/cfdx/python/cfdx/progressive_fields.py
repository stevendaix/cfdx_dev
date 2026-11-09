"""Chunk-limited HDF5 field inspection for industrial CFDX projects.

This is *data streaming* (metadata, numeric statistics and previews). It is not
a VTK volume streaming renderer, which requires piece-aware topology readers.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import h5py
import numpy as np

MAX_BLOCK_BYTES = 8 * 1024 * 1024
MAX_PREVIEW = 2048


@dataclass(frozen=True, slots=True)
class FieldStatistics:
    path: str
    shape: tuple[int, ...]
    units: str
    count: int
    finite_count: int
    minimum: float | None
    maximum: float | None
    mean: float | None
    preview: tuple[float, ...]


def _field(h5, name: str):
    if not name.startswith("/fields/") or ".." in Path(name).parts:
        raise ValueError("field must be an absolute /fields/... HDF5 dataset path")
    link = h5.get(name, getlink=True)
    data = h5.get(name)
    if not isinstance(link, h5py.HardLink) or not isinstance(data, h5py.Dataset):
        raise ValueError("field is not an internal hard-linked dataset")
    if data.is_virtual or data.external:
        raise ValueError("external and virtual fields are not streamed")
    if data.dtype.kind not in "fiu" or data.dtype.itemsize > 8 or data.ndim != 1:
        raise ValueError("only one-dimensional numeric cell fields are supported")
    if data.size == 0:
        raise ValueError("empty field")
    if (
        data.chunks is not None
        and math.prod(data.chunks) * data.dtype.itemsize > MAX_BLOCK_BYTES
    ):
        raise ValueError("HDF5 chunk decompression exceeds the streaming budget")
    return data


def field_catalog(path: str | Path, *, max_fields: int = 4096) -> tuple[str, ...]:
    """Enumerate bounded direct scalar fields; no volume arrays are loaded."""
    found = []
    with h5py.File(path, "r") as h5:
        root = h5.get("/fields", getlink=True)
        if not isinstance(root, h5py.HardLink):
            return ()
        for prefix in ("/fields/scalar", "/fields/cell"):
            link = h5.get(prefix, getlink=True)
            if link is None:
                continue
            if not isinstance(link, h5py.HardLink):
                raise ValueError("field catalog group must be internal and hard-linked")
            group = h5[prefix]
            if not isinstance(group, h5py.Group):
                continue
            if len(group) > max_fields:
                raise ValueError("field catalog exceeds metadata budget")
            for name in group.keys():
                full = f"{prefix}/{name}"
                try:
                    _field(h5, full)
                except (ValueError, KeyError):
                    continue
                found.append(full)
                if len(found) > max_fields:
                    raise ValueError("field catalog exceeds metadata budget")
    return tuple(sorted(found))


def stream_field_statistics(
    path: str | Path,
    name: str,
    *,
    max_block_bytes: int = MAX_BLOCK_BYTES,
    preview_count: int = 64,
) -> FieldStatistics:
    if not (1024 <= max_block_bytes <= MAX_BLOCK_BYTES) or not (
        1 <= preview_count <= MAX_PREVIEW
    ):
        raise ValueError("invalid streaming memory or preview budget")
    with h5py.File(path, "r") as h5:
        field = _field(h5, name)
        before = (field.shape, field.dtype.str)
        block = max(1, max_block_bytes // field.dtype.itemsize)
        n, finite_count, total, minimum, maximum = (
            int(field.size),
            0,
            0.0,
            math.inf,
            -math.inf,
        )
        samples = []
        interval = max(1, math.ceil(n / preview_count))
        for start in range(0, n, block):
            end = min(n, start + block)
            values = np.asarray(field[start:end], dtype=np.float64)
            finite = values[np.isfinite(values)]
            if finite.size:
                finite_count += int(finite.size)
                total += math.fsum(float(v) for v in finite)
                minimum = min(minimum, float(finite.min()))
                maximum = max(maximum, float(finite.max()))
            first = ((start + interval - 1) // interval) * interval
            for idx in range(first, end, interval):
                if len(samples) < preview_count:
                    samples.append(float(values[idx - start]))
        if (field.shape, field.dtype.str) != before:
            raise RuntimeError("HDF5 field changed during inspection")
        units = field.attrs.get("units", "unknown")
        if isinstance(units, bytes):
            units = units.decode("utf-8")
        if not isinstance(units, str) or len(units) > 128:
            raise ValueError("invalid field units")
        return FieldStatistics(
            name,
            (n,),
            units,
            n,
            finite_count,
            minimum if finite_count else None,
            maximum if finite_count else None,
            total / finite_count if finite_count else None,
            tuple(samples),
        )
