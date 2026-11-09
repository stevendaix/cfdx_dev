"""Fail-closed comparison of physically compatible CFDX HDF5 scalar fields.

No implicit interpolation, regridding, or time/units conversion is performed.
All numerical reductions are chunked and avoid materialising full datasets.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import h5py
import numpy as np

from .progressive_fields import MAX_BLOCK_BYTES, _field

_MESH_DATASETS = (
    "points",
    "face_vertices",
    "face_offsets",
    "owner",
    "neighbour",
    "cell_faces",
    "cell_offsets",
)


class IncompatibleSimulations(ValueError):
    """Cannot compare without an independently verified common mesh or units."""


def _mesh_sha256(h5: h5py.File) -> str:
    """Stable topology/geometry identity from actual byte content, not names."""
    digest = hashlib.sha256()
    for name in _MESH_DATASETS:
        if name not in h5:
            raise IncompatibleSimulations(f"missing mesh dataset: {name}")
        link = h5.get(name, getlink=True)
        dataset = h5[name]
        if (
            not isinstance(link, h5py.HardLink)
            or not isinstance(dataset, h5py.Dataset)
            or dataset.is_virtual
            or dataset.external
        ):
            raise IncompatibleSimulations("untrusted mesh topology")
        digest.update(name.encode("ascii"))
        digest.update(str(dataset.dtype).encode("ascii"))
        digest.update(str(tuple(dataset.shape)).encode("ascii"))
        if dataset.ndim not in (1, 2):
            raise IncompatibleSimulations("unsupported mesh topology rank")
        if (
            dataset.chunks is not None
            and math.prod(dataset.chunks) * dataset.dtype.itemsize > MAX_BLOCK_BYTES
        ):
            raise IncompatibleSimulations(
                "mesh chunk exceeds HDF5 decompression budget"
            )
        count = max(
            1,
            MAX_BLOCK_BYTES
            // max(
                1,
                dataset.dtype.itemsize
                * (dataset.shape[-1] if dataset.ndim == 2 else 1),
            ),
        )
        for start in range(0, dataset.shape[0], count):
            digest.update(np.asarray(dataset[start : start + count]).tobytes())
    return digest.hexdigest()


@dataclass(frozen=True)
class ComparisonReport:
    field: str
    units: str
    cells: int
    mesh_sha256: str
    time: float | None
    l2_difference: float
    rms_difference: float
    linf_difference: float
    relative_l2_difference: float | None
    compared_points: int
    nonfinite: int
    certificate: bool = False

    def to_dict(self):
        return asdict(self)


def compare_scalar_fields(
    left: str | Path,
    right: str | Path,
    field: str,
    *,
    block_bytes: int = MAX_BLOCK_BYTES,
) -> ComparisonReport:
    if not (1024 <= block_bytes <= MAX_BLOCK_BYTES):
        raise ValueError("comparison block budget must be 1 KiB..8 MiB")
    with h5py.File(left, "r") as a, h5py.File(right, "r") as b:
        mesh_a = _mesh_sha256(a)
        mesh_b = _mesh_sha256(b)
        if mesh_a != mesh_b:
            raise IncompatibleSimulations(
                "mesh topology or geometry differs; conservative remapping required"
            )
        arr_a, arr_b = _field(a, field), _field(b, field)
        # These are cell fields, not arbitrary arrays with a matching name.
        # A mismatch could otherwise make an incorrect pairwise subtraction
        # look like a valid CFD comparison.
        for h5, array in ((a, arr_a), (b, arr_b)):
            cell_offsets = h5["cell_offsets"]
            if (
                cell_offsets.ndim != 1
                or cell_offsets.size < 2
                or array.size != cell_offsets.size - 1
            ):
                raise IncompatibleSimulations(
                    "cell field size does not match mesh cell_offsets"
                )
        if arr_a.shape != arr_b.shape:
            raise IncompatibleSimulations("field shapes differ")
        unit_a = arr_a.attrs.get("units", None)
        unit_b = arr_b.attrs.get("units", None)
        if unit_a is None or unit_b is None or unit_a != unit_b:
            raise IncompatibleSimulations("field units missing or incompatible")
        if isinstance(unit_a, bytes):
            unit_a = unit_a.decode("utf-8")
        time_a, time_b = (
            a.attrs.get("physical_time", None),
            b.attrs.get("physical_time", None),
        )
        if time_a is not None and (
            isinstance(time_a, (str, bytes)) or not math.isfinite(float(time_a))
        ):
            raise IncompatibleSimulations("physical time must be finite")
        if time_b is not None and (
            isinstance(time_b, (str, bytes)) or not math.isfinite(float(time_b))
        ):
            raise IncompatibleSimulations("physical time must be finite")
        if (time_a is None) != (time_b is None) or (
            time_a is not None and float(time_a) != float(time_b)
        ):
            raise IncompatibleSimulations("physical times incompatible")
        sum_sq, reference_sq, max_abs, nonfinite = 0.0, 0.0, 0.0, 0
        n = int(arr_a.size)
        count = max(1, block_bytes // (8 * 3))
        for start in range(0, n, count):
            lhs = np.asarray(arr_a[start : start + count], dtype=np.float64)
            rhs = np.asarray(arr_b[start : start + count], dtype=np.float64)
            valid = np.isfinite(lhs) & np.isfinite(rhs)
            nonfinite += int((~valid).sum())
            if not valid.all():
                raise IncompatibleSimulations(
                    "non-finite comparison values: no trustworthy field difference"
                )
            delta = rhs - lhs
            if not np.isfinite(delta).all():
                raise IncompatibleSimulations("field difference overflows")
            sum_sq += float(np.dot(delta, delta))
            reference_sq += float(np.dot(lhs, lhs))
            max_abs = max(max_abs, float(np.max(np.abs(delta))))
        return ComparisonReport(
            field,
            str(unit_a),
            n,
            mesh_a,
            float(time_a) if time_a is not None else None,
            math.sqrt(sum_sq),
            math.sqrt(sum_sq / n),
            max_abs,
            math.sqrt(sum_sq / reference_sq) if reference_sq > 0 else None,
            n,
            nonfinite,
        )


def _read_report(path: Path, *, max_bytes: int = 16 * 1024 * 1024) -> dict:
    """Bounded, race-checked native JSON snapshot, not a live run certificate."""
    import json
    import stat

    from .convergence_evidence import _reject_constant, _unique_object

    if path.is_symlink():
        raise IncompatibleSimulations("symbolic-link solver report is not supported")
    try:
        before = path.stat()
        if (
            not stat.S_ISREG(before.st_mode)
            or before.st_size < 8
            or before.st_size > max_bytes
        ):
            raise IncompatibleSimulations("report is missing, unsafe or oversized")
        with path.open("rb") as source:
            raw = source.read(max_bytes + 1)
        after = path.stat()
    except OSError as exc:
        raise IncompatibleSimulations(f"cannot read native report: {exc}") from exc
    if (
        len(raw) != before.st_size
        or len(raw) > max_bytes
        or (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns)
        != (after.st_dev, after.st_ino, after.st_size, after.st_mtime_ns)
    ):
        raise IncompatibleSimulations("solver report changed while reading")
    try:
        value = json.loads(
            raw, object_pairs_hook=_unique_object, parse_constant=_reject_constant
        )
    except (json.JSONDecodeError, UnicodeError, ValueError) as exc:
        raise IncompatibleSimulations(f"invalid native report: {exc}") from exc
    if not isinstance(value, dict):
        raise IncompatibleSimulations("native report must be an object")
    return value


def compare_run_reports(left_dir: str | Path, right_dir: str | Path) -> dict:
    """Compare recorded solver telemetry, and optionally dimensional wall forces.

    A comparison of telemetry is NOT a verified comparison of the same
    geometry/physics or a physical validation. Workbench never invents elapsed
    time or force coefficients absent from the native artifacts.
    """
    left_dir, right_dir = Path(left_dir), Path(right_dir)
    reports = []
    for root in (left_dir, right_dir):
        execution = _read_report(root / "execution.json")
        convergence = _read_report(root / "convergence.json")
        if (
            execution.get("format") != "CFDX-EXECUTION"
            or convergence.get("format") != "CFDX-CONVERGENCE"
            or execution.get("schema_version") != 1
            or convergence.get("schema_version") != 1
            or not isinstance(execution.get("converged"), bool)
            or not isinstance(convergence.get("converged"), bool)
            or execution.get("iterations") != convergence.get("iterations")
        ):
            raise IncompatibleSimulations(
                "native solver report schemas or iteration counts disagree"
            )
        rows = convergence.get("history")
        if not isinstance(rows, list) or not rows or not isinstance(rows[-1], dict):
            raise IncompatibleSimulations("native solver history is missing")
        reports.append((execution, rows[-1]))
    fields = set(reports[0][1]) & set(reports[1][1])
    data = {}
    for name in sorted(fields):
        if name == "iteration":
            continue
        lhs, rhs = reports[0][1][name], reports[1][1][name]
        if any(
            isinstance(v, bool)
            or not isinstance(v, (int, float))
            or not math.isfinite(v)
            for v in (lhs, rhs)
        ):
            continue
        delta = float(rhs) - float(lhs)
        if math.isfinite(delta):
            data[name] = {
                "reference": float(lhs),
                "candidate": float(rhs),
                "delta": delta,
            }
    if not data:
        raise IncompatibleSimulations("no comparable recorded final scalar metrics")
    durations = []
    for execution, _ in reports:
        duration = execution.get("elapsed_wall_seconds")
        if duration is None:
            durations.append(None)  # Legacy execution.json has no measured time.
        elif (
            isinstance(duration, bool)
            or not isinstance(duration, (int, float))
            or not math.isfinite(float(duration))
            or duration < 0
        ):
            raise IncompatibleSimulations(
                "elapsed_wall_seconds must be nonnegative and finite"
            )
        else:
            durations.append(float(duration))
    comparison_time = (
        {
            "reference": durations[0],
            "candidate": durations[1],
            "delta": durations[1] - durations[0],
        }
        if all(value is not None for value in durations)
        else "not_recorded_in_native_reports"
    )
    result = {
        "kind": "CFDX-TELEMETRY-COMPARISON",
        "comparison_certificate": False,
        "mesh_compatibility": "not_verified_for_report_only_comparison",
        "reference_iterations": reports[0][0]["iterations"],
        "candidate_iterations": reports[1][0]["iterations"],
        "reference_converged_as_recorded": reports[0][0]["converged"],
        "candidate_converged_as_recorded": reports[1][0]["converged"],
        "final_solver_metrics": data,
        "wall_forces": None,
        "elapsed_wall_seconds": comparison_time,
    }
    files = [root / "wall_forces.json" for root in (left_dir, right_dir)]
    if all(path.is_file() for path in files):
        a, b = (_read_report(path) for path in files)
        if any(
            obj.get("format") != "CFDX-WALL-FORCES"
            or obj.get("schema_version") != 1
            or obj.get("method") != "quadratic_wls_3d"
            for obj in (a, b)
        ):
            raise IncompatibleSimulations(
                "wall force artifact uses an incompatible method"
            )
        for key in (
            "normal_convention",
            "pressure_reference_pa",
            "moment_center_m",
            "force_definition",
            "wall_patches",
        ):
            if a.get(key) != b.get(key):
                raise IncompatibleSimulations(
                    f"wall force reference incompatible: {key}"
                )
        forces = {}
        for field in (
            "pressure_force_n",
            "viscous_force_n",
            "total_force_n",
            "total_moment_nm",
        ):
            left, right = a.get(field), b.get(field)
            if (
                not isinstance(left, list)
                or not isinstance(right, list)
                or len(left) != 3
                or len(right) != 3
                or any(
                    isinstance(value, bool)
                    or not isinstance(value, (int, float))
                    or not math.isfinite(value)
                    for value in (*left, *right)
                )
            ):
                raise IncompatibleSimulations(f"invalid wall force vector {field}")
            forces[field] = {
                "reference": list(map(float, left)),
                "candidate": list(map(float, right)),
                "delta": [float(v - r) for r, v in zip(left, right)],
            }
        result["wall_forces"] = {
            "dimensional_only": True,
            "mesh_and_reference_area_identity": "not_verified",
            "vectors": forces,
        }
    return result
