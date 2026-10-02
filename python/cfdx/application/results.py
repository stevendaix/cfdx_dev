"""Application-level result browsing and display-object contracts."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from ..results_series import ResultSeries


@dataclass(frozen=True)
class ResultFrameState:
    stable_id: str
    path: str
    sequence: int
    time: float | None
    iteration: int | None
    complete: bool
    fields: tuple[str, ...]


@dataclass(frozen=True)
class ResultsState:
    directory: str | None = None
    frames: tuple[ResultFrameState, ...] = ()
    selected_frame_id: str | None = None
    selected_field: str | None = None

    @property
    def field_names(self) -> tuple[str, ...]:
        return tuple(sorted({field for frame in self.frames for field in frame.fields}))

    @property
    def selected_frame(self) -> ResultFrameState | None:
        return next((frame for frame in self.frames if frame.stable_id == self.selected_frame_id), None)


def build_results_state(
    series: ResultSeries | None,
    *,
    directory: str | Path | None = None,
    selected_frame_id: str | None = None,
    selected_field: str | None = None,
) -> ResultsState:
    frames = tuple(
        ResultFrameState(
            stable_id=f"result:{frame.sequence}",
            path=str(frame.path),
            sequence=frame.sequence,
            time=frame.time,
            iteration=frame.iteration,
            complete=frame.complete,
            fields=frame.fields,
        )
        for frame in (series.frames if series is not None else ())
    )
    if selected_frame_id is None and frames:
        selected_frame_id = frames[0].stable_id
    if selected_frame_id not in {frame.stable_id for frame in frames}:
        selected_frame_id = None
    field_names = {field for frame in frames for field in frame.fields}
    if selected_field not in field_names:
        selected_field = next(iter(sorted(field_names)), None)
    return ResultsState(
        directory=str(directory) if directory is not None else None,
        frames=frames,
        selected_frame_id=selected_frame_id,
        selected_field=selected_field,
    )
