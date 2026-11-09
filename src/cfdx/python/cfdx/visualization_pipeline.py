"""Non-destructive, replayable result-visualization operations.

The pipeline keeps the loaded dataset as its source of truth. Every evaluation
starts from that source and applies an ordered sequence of derived operations.
Inputs are validated before any VTK/PyVista algorithm is invoked.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from numbers import Real
from typing import Any


def _finite_number(value: float, label: str) -> float:
    if isinstance(value, bool) or not isinstance(value, Real):
        raise ValueError(f"{label} must be a finite number")  # noqa: TRY004
    try:
        number = float(value)
    except (ValueError, OverflowError) as exc:
        raise ValueError(f"{label} must be a finite number") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} must be a finite number")
    return number


def _vector3(vector: tuple[float, float, float], label: str) -> tuple[float, ...]:
    try:
        if len(vector) != 3:
            raise ValueError(f"{label} must contain three finite coordinates")
        return tuple(_finite_number(value, label) for value in vector)
    except TypeError as exc:
        raise ValueError(f"{label} must contain three finite coordinates") from exc


@dataclass(frozen=True, slots=True)
class VisualizationOperation:
    """A validated operation that derives a dataset from its input.

    Threshold is explicitly cell-preferred (matching PyVista's documented
    default); ``clip(..., inplace=False)`` preserves authoritative source data.
    """

    kind: str
    scalars: str | None = None
    isosurfaces: int = 10
    normal: tuple[float, float, float] = (1.0, 0.0, 0.0)
    factor: float = 1.0
    minimum: float = 0.0
    maximum: float = 1.0
    origin: tuple[float, float, float] = (0.0, 0.0, 0.0)
    invert: bool = False

    def __post_init__(self) -> None:
        if self.kind not in {"contour", "slice", "glyph", "clip", "threshold"}:
            raise ValueError(f"unsupported visualization operation: {self.kind!r}")
        if self.kind in {"contour", "glyph", "threshold"} and (
            not isinstance(self.scalars, str) or not self.scalars.strip()
        ):
            raise ValueError(f"{self.kind} requires a scalar or vector field name")
        if self.kind == "contour" and (
            isinstance(self.isosurfaces, bool)
            or not isinstance(self.isosurfaces, int)
            or not 1 <= self.isosurfaces <= 1000
        ):
            raise ValueError("isosurfaces must be between 1 and 1000")
        if self.kind == "glyph" and _finite_number(self.factor, "glyph factor") <= 0:
            raise ValueError("glyph factor must be positive")
        if self.kind in {"slice", "clip"}:
            vector = _vector3(self.normal, "plane normal")
            if not any(component != 0 for component in vector):
                raise ValueError("plane normal must be a non-zero 3-vector")
        if self.kind == "clip":
            _vector3(self.origin, "plane origin")
        if self.kind == "threshold":
            lower = _finite_number(self.minimum, "threshold minimum")
            upper = _finite_number(self.maximum, "threshold maximum")
            if lower > upper:
                raise ValueError("threshold minimum must not exceed maximum")
        if not isinstance(self.invert, bool):
            raise ValueError("invert must be a boolean")  # noqa: TRY004

    @classmethod
    def contour(cls, scalars: str, isosurfaces: int = 10) -> VisualizationOperation:
        return cls(kind="contour", scalars=scalars, isosurfaces=isosurfaces)

    @classmethod
    def slice(
        cls, normal: tuple[float, float, float] = (1.0, 0.0, 0.0)
    ) -> VisualizationOperation:
        return cls(kind="slice", normal=normal)

    @classmethod
    def glyph(cls, scalars: str, factor: float = 1.0) -> VisualizationOperation:
        return cls(kind="glyph", scalars=scalars, factor=factor)

    @classmethod
    def clip(
        cls,
        normal: tuple[float, float, float] = (1.0, 0.0, 0.0),
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
        *,
        invert: bool = True,
    ) -> VisualizationOperation:
        return cls(kind="clip", normal=normal, origin=origin, invert=invert)

    @classmethod
    def threshold(
        cls, scalars: str, minimum: float, maximum: float, *, invert: bool = False
    ) -> VisualizationOperation:
        return cls(
            kind="threshold",
            scalars=scalars,
            minimum=minimum,
            maximum=maximum,
            invert=invert,
        )

    def apply(self, dataset: Any) -> Any:
        """Return a derived dataset without mutating the supplied source."""
        fields = set(dataset.point_data.keys()) | set(dataset.cell_data.keys())
        if self.scalars is not None and self.scalars not in fields:
            raise KeyError(self.scalars)
        if self.kind == "contour":
            return dataset.contour(isosurfaces=self.isosurfaces, scalars=self.scalars)
        if self.kind == "slice":
            return dataset.slice(normal=self.normal)
        if self.kind == "clip":
            return dataset.clip(
                normal=self.normal,
                origin=self.origin,
                invert=self.invert,
                inplace=False,
            )
        if self.kind == "threshold":
            return dataset.threshold(
                value=(self.minimum, self.maximum),
                scalars=self.scalars,
                invert=self.invert,
                preference="cell",
            )
        array = dataset.point_data.get(self.scalars)
        if array is None:
            array = dataset.cell_data.get(self.scalars)
        if array is None:
            raise KeyError(self.scalars)
        kwargs: dict[str, Any] = {"scale": self.scalars, "factor": self.factor}
        if getattr(array, "ndim", 1) == 2 and array.shape[1] == 3:
            kwargs["orient"] = self.scalars
        return dataset.glyph(**kwargs)


@dataclass(frozen=True, slots=True)
class VisualizationPipeline:
    """Immutable ordered operations evaluated against an original dataset."""

    operations: tuple[VisualizationOperation, ...] = ()

    def append(self, operation: VisualizationOperation) -> VisualizationPipeline:
        if not isinstance(operation, VisualizationOperation):
            raise TypeError("operation must be a VisualizationOperation")
        return replace(self, operations=(*self.operations, operation))

    def remove(self, index: int) -> VisualizationPipeline:
        if index < 0 or index >= len(self.operations):
            raise IndexError(index)
        return replace(
            self,
            operations=self.operations[:index] + self.operations[index + 1 :],
        )

    def replace(
        self, index: int, operation: VisualizationOperation
    ) -> VisualizationPipeline:
        if index < 0 or index >= len(self.operations):
            raise IndexError(index)
        if not isinstance(operation, VisualizationOperation):
            raise TypeError(
                f"operation must be a VisualizationOperation, got {operation!r}"
            )
        operations = list(self.operations)
        operations[index] = operation
        return replace(self, operations=tuple(operations))

    def move(self, from_index: int, to_index: int) -> VisualizationPipeline:
        if from_index < 0 or from_index >= len(self.operations):
            raise IndexError(from_index)
        if to_index < 0 or to_index >= len(self.operations):
            raise IndexError(to_index)
        operations = list(self.operations)
        operations.insert(to_index, operations.pop(from_index))
        return replace(self, operations=tuple(operations))

    def clear(self) -> VisualizationPipeline:
        return VisualizationPipeline()

    def evaluate(self, source: Any) -> Any:
        """Replay all operations from the unchanged source dataset."""
        result = source
        for operation in self.operations:
            result = operation.apply(result)
        return result
