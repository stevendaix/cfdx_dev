"""Headless discovery of CFD result time series."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import re

_SUPPORTED={" .vtu".strip(),".vtk",".vtp",".pvtu"}
_NUMBER=re.compile(r"(?<![A-Za-z])(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?(?![A-Za-z])")

@dataclass(frozen=True)
class ResultFrame:
    path: Path
    sequence: int
    time: float | None = None
    complete: bool = True
    fields: tuple[str,...] = ()

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

def discover_result_series(directory: Path, *, inspect_fields: bool = False) -> ResultSeries:
    """Discover VTK-family files; keep files being written as incomplete frames."""
    directory=Path(directory)
    if not directory.is_dir(): raise NotADirectoryError(directory)
    paths=sorted((p for p in directory.iterdir() if p.is_file() and p.suffix.lower() in _SUPPORTED),key=_sort_key)
    frames=[]
    for p in paths:
        complete=p.stat().st_size > 0
        fields=()
        if inspect_fields and complete:
            try:
                import pyvista as pv
                dataset=pv.read(p)
                fields=tuple(sorted(set(dataset.point_data.keys()) | set(dataset.cell_data.keys())))
            except (ImportError,OSError,RuntimeError,ValueError):
                complete=False
        value=_sort_key(p)[0]
        frames.append(ResultFrame(p,len(frames),value if value != float("inf") else None,complete,fields))
    return ResultSeries(tuple(frames))
