"""Memory-mapped, bounded CFDX native per-cell residual maps (CFDXMAP1).

Diagnostic only. No normalization to a solver convergence certificate.
Every row follows the original CFDX cell ordering on the originating mesh.
"""

from __future__ import annotations

import math
import mmap
import os
import stat
import struct
from dataclasses import dataclass
from pathlib import Path

import numpy as np

MAX_MAP_BYTES = 512 * 1024 * 1024 + 8192
ROW_DTYPE = np.dtype(
    [("momentum", "<f8"), ("mass", "<f8"), ("streak", "<u2"), ("valid", "u1")],
    align=False,
)
assert ROW_DTYPE.itemsize == 19


class CellMapError(ValueError):
    """Binary output is malformed, inconsistent, or too expensive to read."""


@dataclass(frozen=True)
class CellMapFrame:
    iteration: int
    offset: int
    operator_comparable: bool


class CellFieldMap:
    """Read-only, memory-mapped observations; uncommitted tail is omitted."""

    def __init__(self, path: str | Path, max_bytes: int = MAX_MAP_BYTES):
        self.path = Path(path)
        self._mapping = None
        flags = (
            os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
        )
        try:
            fd = os.open(self.path, flags)
        except OSError as exc:
            raise CellMapError(f"cannot safely open map: {exc}") from exc
        try:
            info = os.fstat(fd)
            if (
                not stat.S_ISREG(info.st_mode)
                or info.st_size < 36
                or info.st_size > max_bytes
            ):
                raise CellMapError("map must be a regular file within its size budget")
            self._mapping = mmap.mmap(fd, 0, access=mmap.ACCESS_READ)
        finally:
            os.close(fd)
        try:
            self._parse()
        except Exception:
            self.close()
            raise

    def __enter__(self):
        return self

    def __exit__(self, *_exc):
        self.close()

    def close(self):
        if self._mapping is not None:
            self._mapping.close()
            self._mapping = None

    def _parse(self):
        raw = self._mapping
        if raw[:8] != b"CFDXMAP1":
            raise CellMapError("unsupported CFDX cell map format")
        self.n_cells, self.n_faces, size, self.reference_tolerance = struct.unpack_from(
            "<QQId", raw, 8
        )
        if not (
            0 < self.n_cells <= 5_000_000
            and self.n_faces > 0
            and 0 < size <= 4096
            and math.isfinite(self.reference_tolerance)
            and self.reference_tolerance > 0
        ):
            raise CellMapError("invalid cell map header or tolerance")
        pos = 36 + size
        if pos > len(raw):
            raise CellMapError("truncated cell map header")
        self.mesh_source = raw[36:pos].decode("utf-8")
        if not self.mesh_source:
            raise CellMapError("empty mesh provenance")
        self.frames = []
        self.complete = False
        self.accepted_iterations = None
        previous = 0
        record_bytes = self.n_cells * ROW_DTYPE.itemsize
        while pos < len(raw):
            marker = raw[pos : pos + 4]
            if marker == b"CFND":
                if pos + 12 > len(raw):
                    break  # live, unfinished end marker
                (accepted,) = struct.unpack_from("<Q", raw, pos + 4)
                if accepted < previous or pos + 12 != len(raw):
                    raise CellMapError("invalid end marker")
                self.complete = True
                self.accepted_iterations = accepted
                pos += 12
                break
            if marker != b"FRME":
                # Only an incomplete trailing fragment can be ignored.
                if pos + 4 > len(raw):
                    break
                raise CellMapError("invalid cell frame marker")
            end = pos + 13 + record_bytes + 4
            if end > len(raw):
                break  # half-written frame is never visible
            (iteration,) = struct.unpack_from("<Q", raw, pos + 4)
            comparable = raw[pos + 12]
            if (
                iteration <= previous
                or comparable not in (0, 1)
                or raw[end - 4 : end] != b"DONE"
            ):
                raise CellMapError("invalid or uncommitted cell frame")
            self.frames.append(CellMapFrame(iteration, pos + 13, bool(comparable)))
            previous = iteration
            pos = end
        # Closing a frame stream does NOT establish physical or solver convergence.

    def frame(self, iteration: int | None = None) -> CellMapFrame:
        if not self.frames:
            raise CellMapError("no committed full-cell map frame")
        if iteration is None:
            return self.frames[-1]
        for frame in self.frames:
            if frame.iteration == iteration:
                return frame
        raise CellMapError("requested iteration has no full-cell snapshot")

    def rows(self, iteration: int | None = None):
        frame = self.frame(iteration)
        rows = np.frombuffer(
            self._mapping, dtype=ROW_DTYPE, count=self.n_cells, offset=frame.offset
        )
        if np.any(rows["valid"] > 1):
            raise CellMapError("invalid per-cell valid flag")
        valid = rows["valid"] == 1
        if np.any(~np.isfinite(rows["momentum"][valid])) or np.any(
            rows["momentum"][valid] < 0
        ):
            raise CellMapError("invalid full-cell momentum residual")
        if np.any(~np.isfinite(rows["mass"][valid])) or np.any(rows["mass"][valid] < 0):
            raise CellMapError("invalid full-cell mass defect")
        return rows

    def scalar(self, metric: str, iteration: int | None = None):
        """Return one finite-or-NaN numeric field and a fixed display range.

        Momentum: log10(local residual / fixed recorded reference tolerance).
        This is a diagnostic threshold, not a local solution-error estimate.
        """
        rows = self.rows(iteration)
        valid = rows["valid"] == 1
        values = np.full(self.n_cells, np.nan, dtype=np.float64)
        if metric == "momentum":
            values[valid] = np.log10(
                np.maximum(rows["momentum"][valid] / self.reference_tolerance, 1e-12)
            )
            return values, (-2.0, 4.0), "log10(momentum defect / reference)"
        if metric == "mass":
            values[valid] = np.log10(np.maximum(rows["mass"][valid], 1e-18))
            return (
                values,
                (-12.0, 0.0),
                "log10(|mass imbalance|) [units depend on case]",
            )
        if metric == "persistence":
            values[valid] = rows["streak"][valid]
            return values, (0.0, 10.0), "non-decreasing accepted transitions"
        raise CellMapError("unknown cell map metric")

    def summary(self, iteration: int | None = None):
        rows = self.rows(iteration)
        valid = rows["valid"] == 1
        scores = rows["momentum"][valid]
        return {
            "format": "CFDXMAP1",
            "iteration": self.frame(iteration).iteration,
            "n_cells": self.n_cells,
            "assessed": int(valid.sum()),
            "above_reference": int(np.sum(scores > self.reference_tolerance)),
            "persistent": int(np.sum(rows["streak"][valid] >= 3)),
            "unassessed": int((~valid).sum()),
            "reference_tolerance": self.reference_tolerance,
            "operator_comparable": self.frame(iteration).operator_comparable,
            "mesh_source": self.mesh_source,
            "convergence_certificate": False,
        }
