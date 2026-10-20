"""Optional foamlib-backed inspection and result reader for OpenFOAM cases.

This module intentionally does not translate solver-specific physics into CFDX.
It provides a robust parsing boundary for dictionaries and fields; conversion to
CFDX's canonical case schema remains the responsibility of the I/O adapter.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any


class OpenFOAMDependencyError(ImportError):
    """Raised when the optional foamlib dependency is unavailable."""


def _foamlib():
    try:
        from foamlib import FoamCase
    except ImportError as exc:
        raise OpenFOAMDependencyError(
            "OpenFOAM dictionary/result inspection requires the optional "
            "'foamlib' package. Install it with: python -m pip install foamlib"
        ) from exc
    return FoamCase


@dataclass(frozen=True)
class OpenFOAMFieldSummary:
    """Metadata for a field available at one OpenFOAM time."""

    name: str
    dimensions: tuple[int, ...] | None
    boundary_types: dict[str, str]


@dataclass(frozen=True)
class OpenFOAMTimeSummary:
    """Fields available in a time directory."""

    name: str
    time: float
    fields: tuple[OpenFOAMFieldSummary, ...]


@dataclass(frozen=True)
class OpenFOAMCaseSummary:
    """Discoverable configuration and result metadata for a case."""

    path: Path
    application: str | None
    times: tuple[OpenFOAMTimeSummary, ...]
    dictionaries: tuple[str, ...]


def inspect_openfoam_case(case_path: str | Path) -> OpenFOAMCaseSummary:
    """Inspect case dictionaries, time directories, and field metadata.

    No input files are modified. Unsupported directives or field encodings are
    reported by foamlib as parsing errors rather than silently ignored here.
    """
    path = Path(case_path).expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(f"OpenFOAM case directory does not exist: {path}")

    FoamCase = _foamlib()
    case = FoamCase(path)

    application: str | None = None
    try:
        value = case.control_dict["application"]
        application = str(value) if value is not None else None
    except (FileNotFoundError, KeyError):
        pass

    dictionary_paths = []
    for folder in ("system", "constant"):
        root = path / folder
        if root.is_dir():
            dictionary_paths.extend(
                str(item.relative_to(path))
                for item in sorted(root.rglob("*"))
                if item.is_file() and not item.name.startswith(".")
            )

    time_summaries = []
    for time_directory in case:
        field_summaries = []
        for field in time_directory:
            boundary_types: dict[str, str] = {}
            try:
                for patch, condition in field.boundary_field.items():
                    condition_type = condition.get("type")
                    if condition_type is not None:
                        boundary_types[str(patch)] = str(condition_type)
            except (AttributeError, TypeError):
                # Preserve field discovery even for a boundary representation
                # that the installed foamlib version cannot expose uniformly.
                pass

            dimensions = None
            try:
                raw_dimensions = field.dimensions
                dimensions = tuple(int(value) for value in raw_dimensions)
            except (AttributeError, KeyError, TypeError, ValueError):
                pass

            field_summaries.append(
                OpenFOAMFieldSummary(
                    name=Path(field).name,
                    dimensions=dimensions,
                    boundary_types=boundary_types,
                )
            )

        time_summaries.append(
            OpenFOAMTimeSummary(
                name=time_directory.name,
                time=float(time_directory.time),
                fields=tuple(sorted(field_summaries, key=lambda item: item.name)),
            )
        )

    return OpenFOAMCaseSummary(
        path=path,
        application=application,
        times=tuple(sorted(time_summaries, key=lambda item: item.time)),
        dictionaries=tuple(dictionary_paths),
    )


def read_openfoam_field(
    case_path: str | Path,
    time: str | float,
    field_name: str,
) -> dict[str, Any]:
    """Read one field's internal and boundary values without changing the case.

    The returned values retain foamlib/NumPy representations and OpenFOAM
    dimensions. Callers must explicitly map units, locations, and boundary
    semantics before inserting the data into CFDX's canonical model.
    """
    path = Path(case_path).expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(f"OpenFOAM case directory does not exist: {path}")

    FoamCase = _foamlib()
    case = FoamCase(path)
    time_directory = case[time]
    field = time_directory[field_name]

    try:
        internal = field.internal_field
    except (AttributeError, KeyError) as exc:
        raise ValueError(
            f"Field {field_name!r} at time {time!r} has no readable internalField"
        ) from exc

    boundary: dict[str, Any] = {}
    for patch, condition in field.boundary_field.items():
        boundary[str(patch)] = condition.as_dict() if hasattr(condition, "as_dict") else dict(condition)

    return {
        "name": field_name,
        "time": float(time_directory.time),
        "dimensions": tuple(int(value) for value in field.dimensions),
        "internal": internal,
        "boundary": boundary,
    }


def read_openfoam_time_fields(
    case_path: str | Path,
    time: str | float = "latest",
    *,
    cell_count: int | None = None,
) -> dict[str, dict[str, Any]]:
    """Read scalar and vector volume fields at one time for CFDX import.

    Uniform fields are expanded to cell_count when supplied, so their arrays
    align with cell-centred fields in CFDX. Tensor fields retain their source
    class and raw values; callers must not silently coerce them to vectors.
    """
    import numpy as np

    path = Path(case_path).expanduser().resolve()
    if not path.is_dir():
        raise NotADirectoryError(f"OpenFOAM case directory does not exist: {path}")

    FoamCase = _foamlib()
    case = FoamCase(path)
    if time == "latest":
        if len(case) == 0:
            raise ValueError(f"No OpenFOAM time directories found in {path}")
        time_directory = case[-1]
    else:
        time_directory = case[time]

    fields: dict[str, dict[str, Any]] = {}
    for field in time_directory:
        try:
            source_class = str(field["FoamFile"]["class"])
            raw = np.asarray(field.internal_field, dtype=np.float64)
            dimensions = tuple(int(value) for value in field.dimensions)
        except (KeyError, TypeError, ValueError, AttributeError) as exc:
            raise ValueError(
                f"Cannot decode field {Path(field).name!s} at time "
                f"{time_directory.name}: {exc}"
            ) from exc

        is_uniform = raw.ndim == 0 or (
            raw.ndim == 1 and source_class.endswith("VectorField") and raw.size == 3
        ) or (
            raw.ndim == 1 and source_class.endswith("TensorField") and raw.size in (6, 9)
        )
        if is_uniform and cell_count is not None:
            raw = np.broadcast_to(raw, (cell_count,) + raw.shape).copy()

        boundary = {}
        for patch, condition in field.boundary_field.items():
            try:
                boundary[str(patch)] = condition.as_dict()
            except AttributeError:
                boundary[str(patch)] = dict(condition)

        fields[Path(field).name] = {
            "name": Path(field).name,
            "class": source_class,
            "time": float(time_directory.time),
            "dimensions": dimensions,
            "internal": raw,
            "boundary": boundary,
            "uniform": is_uniform,
        }

    return fields
