"""Headless discovery of CFD result time series."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re

_SUPPORTED={".vtu",".vtk",".vtp",".pvtu"}
_NUMBER=re.compile(r"(?<![A-Za-z])(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![A-Za-z])")
_XML_SUFFIXES=frozenset({".vtu",".vtp",".pvtu"})
_XML_HEADER=b"<?xml"
_LEGACY_HEADER=b"# vtk DataFile Version"

@dataclass(frozen=True)
class ResultFrame:
    """One discovered output file.

    ``complete`` reports whether every readability check that discovery
    actually performed succeeded. Without ``inspect_fields`` that is the
    non-empty size plus the format banner, which rejects empty, truncated and
    non-VTK files but cannot judge a file whose body was cut after a valid
    header. With ``inspect_fields`` a VTK reader is additionally asked, so the
    field reflects a real parse. Discovery never guesses beyond the evidence it
    gathered; it does not treat "unverified" as "usable".
    """

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

def _dataset_is_empty(dataset) -> bool:
    """Report whether a reader produced no points and no cells.

    Missing attributes are treated as unknown rather than empty: readers that
    expose neither count cannot be judged, so the frame is not penalised.
    """
    return (
        getattr(dataset, "n_points", None) == 0
        and getattr(dataset, "n_cells", None) == 0
    )


def _has_recognised_header(path: Path) -> bool:
    """Probe the format banner every supported VTK-family file starts with.

    A size check alone reports a truncated or non-VTK file as usable, because a
    partially written output still has a non-zero length. XML-family files open
    with the XML declaration and legacy .vtk with its version banner, so the
    banner is the cheapest available structural evidence of readability.
    """
    suffix=path.suffix.lower()
    if suffix not in _XML_SUFFIXES and suffix!=".vtk": return True
    try:
        with path.open("rb") as handle: head=handle.read(64).lstrip()
    except OSError: return False
    if suffix in _XML_SUFFIXES: return head.startswith(_XML_HEADER)
    return head.startswith(_LEGACY_HEADER)

def validate_physical_time_provenance(series: ResultSeries) -> None:
    """Require every non-empty result frame to carry authoritative physical time.

    Filename-derived values remain available for browsing legacy output, but a
    solver-facing acceptance path must reject them because filenames are not a
    numerical time provenance contract.
    """
    missing = [
        frame.path.name
        for frame in series.frames
        if frame.time_source != "metadata"
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


def discover_result_series(directory: Path, *, inspect_fields: bool = False, require_physical_time: bool = False) -> ResultSeries:
    """Discover VTK-family files; preserve explicit time/iteration metadata when available."""
    directory=Path(directory)
    if not directory.is_dir(): raise NotADirectoryError(directory)
    paths=sorted((p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in _SUPPORTED),key=_sort_key)
    frames=[]
    for p in paths:
        complete=p.stat().st_size > 0 and _has_recognised_header(p)
        fields=()
        iteration=None
        metadata_time=None
        if inspect_fields and complete:
            try:
                import pyvista as pv
                dataset=pv.read(p)
                # pyvista does not raise on a structurally broken file: the VTK
                # reader logs, warns and hands back an empty dataset. An empty
                # dataset is therefore the only usable unreadability signal.
                if _dataset_is_empty(dataset):
                    complete=False
                fields=tuple(sorted(set(dataset.point_data.keys()) | set(dataset.cell_data.keys())))
                iteration_value = _metadata_scalar(dataset.field_data, ("iteration","Iteration","step","Step"))
                metadata_time = _metadata_scalar(dataset.field_data, ("physical_time","PhysicalTime","time","Time","timeValue","TIME"))
                iteration = int(iteration_value) if iteration_value is not None and iteration_value.is_integer() else None
            except (ImportError,OSError,RuntimeError,ValueError):
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
