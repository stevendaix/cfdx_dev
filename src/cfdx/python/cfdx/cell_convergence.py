"""Strict reader for opt-in per-cell CFDX solver diagnostics.

Records describe a bounded union of cell hotspots, NOT the full residual field.
The 'nondecrease_streak' is computed on the same cell at *every accepted*
nonlinear iteration, never from a sequence of different worst-cell IDs.
None of these diagnostics is a convergence or independent V&V certificate.
"""
from __future__ import annotations

import argparse
import json
import math
import stat
from dataclasses import dataclass
from pathlib import Path

MAX_TRACE_BYTES = 32 * 1024 * 1024
MAX_CELLS_PER_FRAME = 192  # 3 diagnostic top-K collections, each K <= 64


class CellConvergenceError(ValueError):
    """Malformed, oversized, stale or incompatible trace."""


def cell_trace_metric_name(metric: str) -> str:
    """Human-readable scalar name used by the UI and exports."""
    names = {
        "momentum": "Momentum / RHS",
        "mass": "|mass defect|",
        "persistence": "Non-decrease streak",
    }
    try:
        return names[metric]
    except KeyError as exc:
        raise ValueError("metric must be momentum, mass, or persistence") from exc


def sampled_cell_scalar_map(trace: "CellTrace", frame: "CellFrame", metric: str = "momentum") -> tuple[str, tuple[float, ...]]:
    """Build a sparse per-cell scalar map from a sampled frame.

    Cells not present in the bounded top-K snapshot are returned as NaN, so
    callers can render them in a neutral color instead of inventing a value.
    """
    title = cell_trace_metric_name(metric)
    values = [math.nan] * trace.n_cells
    field = {
        "momentum": "momentum_scaled",
        "mass": "mass_imbalance_abs",
        "persistence": "nondecrease_streak",
    }[metric]
    for cell in frame.cells:
        values[cell.cell_id] = float(getattr(cell, field))
    return title, tuple(values)


def _unique_pairs(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise CellConvergenceError(f"duplicate trace JSON key {key}")
        result[key] = value
    return result


def _nonfinite(value: str) -> object:
    raise CellConvergenceError(f"non-finite JSON numeric value: {value}")


def _integer(value: object, name: str, *, min_value: int = 0) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < min_value:
        raise CellConvergenceError(f"{name} must be an integer >= {min_value}")
    return value


def _number(value: object, name: str, *, nonnegative: bool = True) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise CellConvergenceError(f"{name} must be numeric")
    result = float(value)
    if not math.isfinite(result) or (nonnegative and result < 0):
        raise CellConvergenceError(f"{name} must be finite{' and nonnegative' if nonnegative else ''}")
    return result


def _triple(value: object, name: str, *, nonnegative: bool = False) -> tuple[float, float, float]:
    if not isinstance(value, list) or len(value) != 3:
        raise CellConvergenceError(f"{name} must contain three components")
    return tuple(_number(part, name, nonnegative=nonnegative) for part in value)  # type: ignore[return-value]


@dataclass(frozen=True, slots=True)
class CellSample:
    cell_id: int
    center: tuple[float, float, float]
    momentum_abs: tuple[float, float, float]
    momentum_scaled: float
    mass_imbalance_abs: float
    nondecrease_streak: int

    @property
    def stagnant(self) -> bool:
        return self.nondecrease_streak >= 3

    def to_dict(self) -> dict[str, object]:
        return {
            "cell_id": self.cell_id, "center": list(self.center),
            "momentum_abs": list(self.momentum_abs),
            "momentum_scaled": self.momentum_scaled,
            "mass_imbalance_abs": self.mass_imbalance_abs,
            "nondecrease_streak": self.nondecrease_streak,
            "classification": "nondecreasing_3+" if self.stagnant else "not_established",
        }


@dataclass(frozen=True, slots=True)
class CellFrame:
    iteration: int
    assessed: int
    nonfinite: int
    momentum_above_reference: int
    persistent_momentum_cells: int
    operator_comparable: bool
    rhs_scale: tuple[float, float, float]
    cells: tuple[CellSample, ...]

    def ranked(self, by: str = "momentum", limit: int = 30) -> tuple[CellSample, ...]:
        cell_trace_metric_name(by)
        if not isinstance(limit, int) or not 1 <= limit <= 192:
            raise ValueError("limit must be within 1..192")
        if by == "momentum":
            key = lambda c: (c.momentum_scaled, c.nondecrease_streak)
        elif by == "mass":
            key = lambda c: (c.mass_imbalance_abs, c.momentum_scaled)
        else:
            key = lambda c: (c.nondecrease_streak, c.momentum_scaled)
        return tuple(sorted(self.cells, key=key, reverse=True)[:limit])


@dataclass(frozen=True, slots=True)
class CellTrace:
    mesh_source: str
    n_cells: int
    n_faces: int
    every: int
    top_k: int
    frames: tuple[CellFrame, ...]
    complete: bool
    accepted_iterations: int | None
    sha256: str  # integrity of bytes read, NOT a solver/mesh certification
    file_path: str

    @property
    def latest(self) -> CellFrame | None:
        return self.frames[-1] if self.frames else None

    def observed_history(self, cell_id: int) -> tuple[tuple[int, CellSample], ...]:
        """Only recorded top-K snapshots, NOT every accepted cell state.

        The solver itself maintains per-cell streaks on *all* accepted steps.
        Missing snapshots are unknown, never interpolated or assigned zero.
        """
        if not isinstance(cell_id, int) or isinstance(cell_id, bool) or not 0 <= cell_id < self.n_cells:
            raise ValueError("cell ID outside the originating mesh")
        return tuple(
            (frame.iteration, sample)
            for frame in self.frames
            for sample in frame.cells
            if sample.cell_id == cell_id
        )

    def to_dict(self) -> dict[str, object]:
        frame = self.latest
        return {
            "format": "CFDX-CELL-CONVERGENCE", "schema_version": 1,
            "mesh_source": self.mesh_source, "n_cells": self.n_cells,
            "n_faces": self.n_faces, "every": self.every,
            "snapshots": len(self.frames), "complete": self.complete,
            "accepted_iterations": self.accepted_iterations,
            "last_iteration": frame.iteration if frame else None,
            "last_assessed": frame.assessed if frame else None,
            "last_nonfinite": frame.nonfinite if frame else None,
            "last_persistent_cells": frame.persistent_momentum_cells if frame else None,
            "latest_hotspots": [cell.to_dict() for cell in frame.ranked(limit=25)] if frame else [],
            "trace_sha256": self.sha256,
            "solver_convergence_certificate": False,
            "spatial_coverage": "bounded top-K snapshots; full-cell persistence counters computed within solver",
        }


def _read_bounded_file(path: Path, max_bytes: int) -> bytes:
    import os
    if max_bytes < 1024 or max_bytes > 512 * 1024 * 1024:
        raise ValueError("trace read budget outside supported bounds")
    # No symlink or special-file reads: a FIFO could hang the GUI forever.
    flags = os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_NONBLOCK", 0)
    try:
        fd = os.open(path, flags)
    except OSError as exc:
        raise CellConvergenceError(f"cannot open trace safely: {exc}") from exc
    try:
        before = os.fstat(fd)
        if not stat.S_ISREG(before.st_mode) or before.st_size > max_bytes:
            raise CellConvergenceError("trace is not a regular file or exceeds the read budget")
        with os.fdopen(fd, "rb", closefd=False) as stream:
            data = stream.read(max_bytes + 1)
        after = os.fstat(fd)
        if len(data) > max_bytes or not stat.S_ISREG(after.st_mode):
            raise CellConvergenceError("trace exceeds the byte budget")
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
            after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns
        ):
            raise CellConvergenceError("trace changed during the bounded read; retry")
        if len(data) != before.st_size:
            raise CellConvergenceError("trace changed or was truncated during read")
        return data
    finally:
        os.close(fd)


def read_cell_trace(path: str | Path, *, max_bytes: int = MAX_TRACE_BYTES) -> CellTrace:
    """Read live JSONL or finished output. Ignore *only* a trailing partial line."""
    import hashlib

    source = Path(path)
    raw = _read_bounded_file(source, max_bytes)
    if not raw:
        raise CellConvergenceError("empty diagnostic stream")
    lines = raw.split(b"\n")
    # A producer may still be writing the last record. Such a record has NOT
    # been committed to JSONL; never interpret it as a valid frame.
    lines.pop()  # empty final terminator, or incomplete record while writing
    if not lines:
        raise CellConvergenceError("diagnostic header not committed")

    def parse(row: bytes, number: int) -> dict:
        if len(row) > 128 * 1024:
            raise CellConvergenceError(f"trace line {number} exceeds 128 KiB")
        try:
            value = json.loads(row.decode("utf-8"), object_pairs_hook=_unique_pairs,
                               parse_constant=_nonfinite)
        except (UnicodeError, json.JSONDecodeError) as exc:
            raise CellConvergenceError(f"invalid trace line {number}: {exc}") from exc
        if not isinstance(value, dict):
            raise CellConvergenceError(f"trace line {number} must be an object")
        return value

    header = parse(lines[0], 1)
    if (header.get("kind") != "header" or header.get("format") != "CFDX-CELL-CONVERGENCE"
            or header.get("schema_version") != 1):
        raise CellConvergenceError("not a supported CFDX cell convergence trace")
    mesh_source = header.get("mesh_source")
    if not isinstance(mesh_source, str) or not mesh_source or len(mesh_source) > 4096:
        raise CellConvergenceError("missing or invalid mesh source")
    n_cells = _integer(header.get("n_cells"), "n_cells", min_value=1)
    n_faces = _integer(header.get("n_faces"), "n_faces", min_value=1)
    every = _integer(header.get("every"), "every", min_value=1)
    top_k = _integer(header.get("top_k"), "top_k", min_value=1)
    if n_cells > 5_000_000 or top_k > 64 or every > 100_000:
        raise CellConvergenceError("diagnostic exceeds native recorder limits")
    frames: list[CellFrame] = []
    finished = False
    accepted: int | None = None
    previous_iteration = 0
    for index, line in enumerate(lines[1:], 2):
        data = parse(line, index)
        kind = data.get("kind")
        if finished:
            raise CellConvergenceError("unexpected data after end-of-run marker")
        if kind == "end":
            accepted = _integer(data.get("iterations"), "end.iterations")
            final_frame = _integer(data.get("last_frame_iteration"), "end.last_frame_iteration")
            if (accepted < previous_iteration or final_frame != previous_iteration or
                    (accepted > 0 and not frames) or
                    (accepted == 0 and frames)):
                raise CellConvergenceError("end marker disagrees with sampled accepted iterations")
            finished = True
            continue
        if kind != "frame":
            raise CellConvergenceError(f"unexpected record kind at line {index}")
        iteration = _integer(data.get("iteration"), "frame.iteration", min_value=1)
        if iteration <= previous_iteration:
            raise CellConvergenceError("frame iterations must increase strictly")
        previous_iteration = iteration
        assessed = _integer(data.get("assessed"), "frame.assessed")
        nonfinite = _integer(data.get("nonfinite"), "frame.nonfinite")
        above = _integer(data.get("momentum_above_reference"), "momentum_above_reference")
        persistent = _integer(data.get("persistent_momentum_cells"), "persistent_momentum_cells")
        if assessed + nonfinite != n_cells or above > assessed or persistent > above:
            raise CellConvergenceError("invalid frame cell counts")
        comparable = data.get("operator_comparable")
        if not isinstance(comparable, bool):
            raise CellConvergenceError("operator comparability flag must be boolean")
        scale = _triple(data.get("rhs_scale"), "rhs_scale", nonnegative=True)
        if any(value <= 0 for value in scale):
            raise CellConvergenceError("RHS scales must be positive")
        entries = data.get("cells")
        if not isinstance(entries, list) or len(entries) > min(n_cells, 3 * top_k):
            raise CellConvergenceError("trace frame exceeds bounded top-K collection")
        cells: list[CellSample] = []
        seen: set[int] = set()
        for row in entries:
            if not isinstance(row, dict):
                raise CellConvergenceError("cell record must be an object")
            cell = _integer(row.get("cell_id"), "cell_id")
            if cell >= n_cells or cell in seen:
                raise CellConvergenceError("duplicate or out-of-range cell ID")
            seen.add(cell)
            center = _triple(row.get("center"), "center")
            momentum = _triple(row.get("momentum_abs"), "momentum_abs", nonnegative=True)
            score = _number(row.get("momentum_scaled"), "momentum_scaled")
            mass = _number(row.get("mass_imbalance_abs"), "mass_imbalance_abs")
            streak = _integer(row.get("nondecrease_streak"), "nondecrease_streak")
            actual = max(component / scale[d] for d, component in enumerate(momentum))
            if not math.isclose(score, actual, rel_tol=1e-9, abs_tol=1e-14):
                raise CellConvergenceError("per-cell normalized momentum does not match component residuals")
            if streak > iteration:
                raise CellConvergenceError("persistence cannot exceed accepted transitions")
            cells.append(CellSample(cell, center, momentum, score, mass, streak))
        frames.append(CellFrame(iteration, assessed, nonfinite, above, persistent,
                                comparable, scale, tuple(cells)))
        if len(frames) > 100_000:
            raise CellConvergenceError("trace frame budget exceeded")
    return CellTrace(mesh_source, n_cells, n_faces, every, top_k, tuple(frames),
                     finished, accepted, hashlib.sha256(raw).hexdigest(), str(source))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Inspect actual CFDX per-cell solver residual hotspots")
    parser.add_argument("trace", type=Path, help="cell_convergence.jsonl from --cell-hotspots")
    parser.add_argument("--metric", choices=("momentum", "mass", "persistence"), default="momentum")
    parser.add_argument("--top", type=int, default=20)
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args(argv)
    try:
        trace = read_cell_trace(args.trace)
        if args.json:
            print(json.dumps(trace.to_dict(), indent=2, sort_keys=True, allow_nan=False))
        else:
            frame = trace.latest
            if frame is None:
                print("Trace initialized; no accepted solver frame recorded yet")
            else:
                print(f"Accepted iteration {frame.iteration}, {frame.assessed}/{trace.n_cells} finite cells; "
                      f"{frame.persistent_momentum_cells} nondecreasing residual histories; "
                      f"trace {'COMPLETE' if trace.complete else 'LIVE / incomplete'}")
                print("Diagnostics only — no convergence certificate; only sampled top-K cells are displayed.")
                for cell in frame.ranked(args.metric, args.top):
                    print(f"  cell {cell.cell_id:8d} xyz={cell.center} "
                          f"momentum={cell.momentum_scaled:.6g} "
                          f"mass={cell.mass_imbalance_abs:.6g} "
                          f"nondecrease={cell.nondecrease_streak}")
    except (CellConvergenceError, OSError, ValueError) as exc:
        parser.exit(2, f"Invalid cell trace: {exc}\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
