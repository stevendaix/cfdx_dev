"""Non-destructive, replayable result-visualization operations.

The pipeline keeps the loaded dataset as its source of truth. Every evaluation
starts from that source and applies an ordered sequence of derived operations.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Any


@dataclass(frozen=True, slots=True)
class VisualizationOperation:
    """A validated operation that derives a dataset from its input."""

    kind: str
    scalars: str | None = None
    isosurfaces: int = 10
    normal: tuple[float, float, float] = (1.0, 0.0, 0.0)
    factor: float = 1.0

    def __post_init__(self) -> None:
        if self.kind not in {"contour", "slice", "glyph"}:
            raise ValueError(f"unsupported visualization operation: {self.kind!r}")
        if self.kind in {"contour", "glyph"} and not self.scalars:
            raise ValueError(f"{self.kind} requires a scalar or vector field name")
        if self.kind == "contour" and self.isosurfaces < 1:
            raise ValueError("isosurfaces must be positive")
        if self.kind == "glyph" and self.factor <= 0:
            raise ValueError("glyph factor must be positive")
        if self.kind == "slice" and (
            len(self.normal) != 3
            or not any(float(value) != 0.0 for value in self.normal)
        ):
            raise ValueError("slice normal must be a non-zero 3-vector")

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

    def apply(self, dataset: Any) -> Any:
        """Return a derived dataset without mutating the supplied source."""
        fields = set(dataset.point_data.keys()) | set(dataset.cell_data.keys())
        if self.scalars is not None and self.scalars not in fields:
            raise KeyError(self.scalars)
        if self.kind == "contour":
            return dataset.contour(isosurfaces=self.isosurfaces, scalars=self.scalars)
        if self.kind == "slice":
            return dataset.slice(normal=self.normal)
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
        return replace(self, operations=(*self.operations, operation))

    def remove(self, index: int) -> VisualizationPipeline:
        if index < 0 or index >= len(self.operations):
            raise IndexError(index)
        return replace(
            self,
            operations=self.operations[:index] + self.operations[index + 1 :],
        )

    def clear(self) -> VisualizationPipeline:
        return VisualizationPipeline()

    def evaluate(self, source: Any) -> Any:
        """Replay all operations from the unchanged source dataset."""
        result = source
        for operation in self.operations:
            result = operation.apply(result)
        return result
