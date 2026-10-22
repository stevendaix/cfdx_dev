from __future__ import annotations

import pytest

from cfdx.visualization_pipeline import (
    VisualizationOperation,
    VisualizationPipeline,
)


class FakeDataset:
    def __init__(self, history=(), fields=("p", "U")):
        self.history = tuple(history)
        vector = type("Vector", (), {"ndim": 2, "shape": (4, 3)})()
        self.point_data = {"p": object(), "U": vector}
        self.cell_data = {}
        self.fields = fields

    def contour(self, *, isosurfaces, scalars):
        return FakeDataset((*self.history, ("contour", isosurfaces, scalars)))

    def slice(self, *, normal):
        return FakeDataset((*self.history, ("slice", normal)))

    def glyph(self, **kwargs):
        return FakeDataset((*self.history, ("glyph", kwargs)))


def test_pipeline_replays_ordered_operations_without_mutating_source():
    source = FakeDataset()
    pipeline = VisualizationPipeline().append(
        VisualizationOperation.contour("p", 7)
    ).append(VisualizationOperation.slice((0.0, 1.0, 0.0)))

    result = pipeline.evaluate(source)

    assert source.history == ()
    assert result.history == (
        ("contour", 7, "p"),
        ("slice", (0.0, 1.0, 0.0)),
    )
    assert len(pipeline.operations) == 2


def test_pipeline_remove_and_clear_return_new_pipelines():
    first = VisualizationOperation.contour("p")
    second = VisualizationOperation.slice()
    pipeline = VisualizationPipeline((first, second))

    shortened = pipeline.remove(0)

    assert pipeline.operations == (first, second)
    assert shortened.operations == (second,)
    assert pipeline.clear().operations == ()
    with pytest.raises(IndexError):
        pipeline.remove(2)


@pytest.mark.parametrize(
    "kwargs",
    [
        {"kind": "warp"},
        {"kind": "contour"},
        {"kind": "contour", "scalars": "p", "isosurfaces": 0},
        {"kind": "glyph", "scalars": "U", "factor": 0.0},
        {"kind": "slice", "normal": (0.0, 0.0, 0.0)},
    ],
)
def test_invalid_operations_are_rejected(kwargs):
    with pytest.raises(ValueError):
        VisualizationOperation(**kwargs)


def test_operation_checks_field_availability_and_supports_vector_glyphs():
    source = FakeDataset()

    with pytest.raises(KeyError, match="temperature"):
        VisualizationOperation.contour("temperature").apply(source)

    result = VisualizationOperation.glyph("U", factor=2.0).apply(source)
    assert result.history == (
        ("glyph", {"scale": "U", "factor": 2.0, "orient": "U"}),
    )
