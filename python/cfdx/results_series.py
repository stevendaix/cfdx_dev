"""Headless discovery of CFD result time series."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re
import xml.etree.ElementTree as ET

_SUPPORTED={".vtu",".vtk",".vtp",".pvtu"}
_NUMBER=re.compile(r"(?<![A-Za-z])(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![A-Za-z])")

@dataclass(frozen=True)
class ResultFrame:
    path: Path
    sequence: int
    time: float | None = None
    complete: bool = True
    fields: tuple[str,...] = ()
    iteration: int | None = None
    time_source: str | None = None

@dataclass(frozen=True)
class ResultSeries:
    frames: tuple[ResultFrame,...]
    @property
    def is_transient(self) -> bool:
        return len(self.frames)>1
    @property
    def paths(self) -> tuple[Path,...]:
        return tuple(frame.path for frame in self.frames)
    def field_names(self) -> tuple[str,...]:
        return tuple(sorted({field for frame in self.frames for field in frame.fields}))
    def frame(self,index: int) -> ResultFrame:
        try: return self.frames[index]
        except IndexError as exc: raise IndexError(f"result frame index out of range: {index}") from exc

def _sort_key(path: Path) -> tuple[float,int,str]:
    values=_NUMBER.findall(path.stem)
    value=float(values[-1]) if values else float("inf")
    return (value,0,path.name)

def _metadata_scalar(field_data, names: tuple[str, ...]) -> float | None:
    for name in names:
        if name not in field_data:
            continue
        value = field_data[name]
        try:
            if hasattr(value, "reshape"):
                value = value.reshape(-1)
            if len(value) != 1:
                continue
            value = float(value[0])
        except (TypeError, ValueError, IndexError):
            continue
        if value == value and value not in (float("inf"), float("-inf")):
            return value
    return None

def validate_physical_time_provenance(series: ResultSeries) -> None:
    """Require every non-empty result frame to carry authoritative physical time.

    Filename-derived values remain available for browsing legacy output, but a
    solver-facing acceptance path must reject them because filenames are not a
    numerical time provenance contract.
    """
    missing = [
        frame.path.name
        for frame in series.frames
        if frame.complete and frame.time_source != "metadata"
    ]
    if missing:
        raise ValueError(
            "authoritative physical-time metadata missing for: " + ", ".join(missing)
        )
    times = [
        frame.time for frame in series.frames
        if frame.complete and frame.time is not None
    ]
    if any(b < a for a, b in zip(times, times[1:])):
        raise ValueError("physical-time metadata must be monotonic")


def _metadata_scalar_from_mapping(values: dict[str, float], names: tuple[str, ...]) -> float | None:
    for name in names:
        value = values.get(name)
        if value is not None and value == value and value not in (float("inf"), float("-inf")):
            return value
    return None


def _inspect_vtu_xml(path: Path) -> tuple[tuple[str, ...], int | None, float | None]:
    """Inspect VTK XML metadata without requiring the optional PyVista stack."""
    root = ET.parse(path).getroot()
    field_data = root.find("./UnstructuredGrid/FieldData")
    metadata: dict[str, float] = {}
    if field_data is not None:
        for item in field_data.findall("DataArray"):
            name = item.attrib.get("Name")
            text = (item.text or "").strip().split()
            if name and len(text) == 1:
                try:
                    value = float(text[0])
                except ValueError:
                    continue
                if value == value and value not in (float("inf"), float("-inf")):
                    metadata[name] = value

    fields: set[str] = set()
    for parent_path in ("./UnstructuredGrid/Piece/PointData", "./UnstructuredGrid/Piece/CellData"):
        parent = root.find(parent_path)
        if parent is not None:
            fields.update(
                item.attrib["Name"]
                for item in parent.findall("DataArray")
                if "Name" in item.attrib
            )
    iteration_value = _metadata_scalar_from_mapping(
        metadata, ("iteration", "Iteration", "step", "Step")
    )
    time_value = _metadata_scalar_from_mapping(
        metadata, ("physical_time", "PhysicalTime", "time", "Time", "timeValue", "TIME")
    )
    iteration = int(iteration_value) if iteration_value is not None and iteration_value.is_integer() else None
    return tuple(sorted(fields)), iteration, time_value


def discover_result_series(directory: Path, *, inspect_fields: bool = False, require_physical_time: bool = False) -> ResultSeries:
    """Discover VTK-family files; preserve explicit time/iteration metadata when available."""
    directory=Path(directory)
    if not directory.is_dir(): raise NotADirectoryError(directory)
    paths=sorted((p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in _SUPPORTED),key=_sort_key)
    frames=[]
    for p in paths:
        complete=p.stat().st_size > 0
        fields=()
        iteration=None
        metadata_time=None
        if inspect_fields and complete:
            try:
                import pyvista as pv
                dataset=pv.read(p)
                fields=tuple(sorted(set(dataset.point_data.keys()) | set(dataset.cell_data.keys())))
                iteration_value = _metadata_scalar(dataset.field_data, ("iteration","Iteration","step","Step"))
                metadata_time = _metadata_scalar(dataset.field_data, ("physical_time","PhysicalTime","time","Time","timeValue","TIME"))
                iteration = int(iteration_value) if iteration_value is not None and iteration_value.is_integer() else None
            except ImportError:
                if p.suffix.lower() == ".vtu":
                    try:
                        fields, iteration, metadata_time = _inspect_vtu_xml(p)
                    except (ET.ParseError, OSError, ValueError):
                        complete=False
                else:
                    complete=False
            except (OSError,RuntimeError,ValueError):
                complete=False
        filename_time=_sort_key(p)[0]
        if filename_time == float("inf"): filename_time=None
        frame_time=metadata_time if metadata_time is not None else filename_time
        time_source="metadata" if metadata_time is not None else ("filename" if filename_time is not None else None)
        frames.append(ResultFrame(p,len(frames),frame_time,complete,fields,iteration,time_source))
    series = ResultSeries(tuple(frames))
    if require_physical_time:
        validate_physical_time_provenance(series)
    return series
