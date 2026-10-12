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


_DISPLAY_KINDS = frozenset({"surface", "contour", "slice", "vector", "streamline", "volume", "probe"})


@dataclass(frozen=True)
class ResultDataset:
    """Authoritative result dataset reference, independent of any renderer."""

    stable_id: str
    frame_id: str
    path: str
    fields: tuple[str, ...]
    time: float | None = None
    iteration: int | None = None

    def field_available(self, field: str) -> bool:
        return field in self.fields


@dataclass(frozen=True)
class DisplayObject:
    """Serializable presentation object referencing authoritative result data."""

    stable_id: str
    kind: str
    dataset_id: str
    field: str | None = None
    parameters: tuple[tuple[str, str], ...] = ()

    def __post_init__(self) -> None:
        if not self.stable_id:
            raise ValueError("display object stable_id must not be empty")
        if self.kind not in _DISPLAY_KINDS:
            raise ValueError(f"unsupported display object kind: {self.kind!r}")
        if not self.dataset_id:
            raise ValueError("display object dataset_id must not be empty")
        if self.field is not None and not self.field:
            raise ValueError("display object field must not be empty")

    def parameter_map(self) -> dict[str, str]:
        return dict(self.parameters)

    def to_dict(self) -> dict[str, object]:
        return {
            "stable_id": self.stable_id,
            "kind": self.kind,
            "dataset_id": self.dataset_id,
            "field": self.field,
            "parameters": {key: value for key, value in self.parameters},
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> "DisplayObject":
        raw_parameters = data.get("parameters", {})
        if not isinstance(raw_parameters, dict):
            raise ValueError("display object parameters must be an object")
        return cls(
            stable_id=str(data["stable_id"]),
            kind=str(data["kind"]),
            dataset_id=str(data["dataset_id"]),
            field=None if data.get("field") is None else str(data["field"]),
            parameters=tuple(sorted((str(k), str(v)) for k, v in raw_parameters.items())),
        )


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
