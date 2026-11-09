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

    def clip(self, **kwargs):
        return FakeDataset((*self.history, ("clip", kwargs)))

    def threshold(self, **kwargs):
        return FakeDataset((*self.history, ("threshold", kwargs)))


def test_pipeline_replays_ordered_operations_without_mutating_source():
    source = FakeDataset()
    pipeline = (
        VisualizationPipeline()
        .append(VisualizationOperation.contour("p", 7))
        .append(VisualizationOperation.slice((0.0, 1.0, 0.0)))
    )

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


def test_pipeline_replace_returns_new_pipeline_and_replays():
    source = FakeDataset()
    first = VisualizationOperation.contour("p", 7)
    second = VisualizationOperation.slice((0.0, 1.0, 0.0))
    pipeline = VisualizationPipeline((first, second))

    updated = pipeline.replace(1, VisualizationOperation.contour("p", 3))

    assert pipeline.operations == (first, second)
    assert updated.operations == (first, VisualizationOperation.contour("p", 3))
    assert updated.evaluate(source).history == (
        ("contour", 7, "p"),
        ("contour", 3, "p"),
    )


def test_pipeline_replace_rejects_invalid_index_and_type():
    pipeline = VisualizationPipeline((VisualizationOperation.slice(),))
    with pytest.raises(IndexError):
        pipeline.replace(1, VisualizationOperation.slice())
    with pytest.raises(TypeError, match="contour"):
        pipeline.replace(0, "contour")


def test_pipeline_move_reorders_operations_immutably():
    first = VisualizationOperation.contour("p")
    second = VisualizationOperation.slice()
    third = VisualizationOperation.glyph("U")
    pipeline = VisualizationPipeline((first, second, third))

    moved = pipeline.move(0, 2)

    assert pipeline.operations == (first, second, third)
    assert moved.operations == (second, third, first)


def test_pipeline_move_rejects_out_of_range_indices():
    pipeline = VisualizationPipeline(
        (VisualizationOperation.contour("p"), VisualizationOperation.slice())
    )
    with pytest.raises(IndexError):
        pipeline.move(-1, 0)
    with pytest.raises(IndexError):
        pipeline.move(0, 2)
    with pytest.raises(IndexError):
        pipeline.move(2, 0)
    with pytest.raises(IndexError):
        pipeline.move(0, -1)


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
    assert result.history == (("glyph", {"scale": "U", "factor": 2.0, "orient": "U"}),)


def test_clip_and_threshold_are_replayable_and_source_safe():
    source = FakeDataset()
    pipeline = (
        VisualizationPipeline()
        .append(VisualizationOperation.clip((0, 1, 0), (0.5, 0, 0), invert=True))
        .append(VisualizationOperation.threshold("p", -2.0, 2.0))
    )
    result = pipeline.evaluate(source)
    assert source.history == ()
    assert result.history == (
        (
            "clip",
            {
                "normal": (0, 1, 0),
                "origin": (0.5, 0, 0),
                "invert": True,
                "inplace": False,
            },
        ),
        (
            "threshold",
            {
                "value": (-2.0, 2.0),
                "scalars": "p",
                "invert": False,
                "preference": "cell",
            },
        ),
    )
    assert pipeline.operations[0].kind == "clip"


@pytest.mark.parametrize(
    "operation",
    [
        lambda: VisualizationOperation.threshold("p", 2, 1),
        lambda: VisualizationOperation.threshold("p", float("nan"), 1),
        lambda: VisualizationOperation.threshold("p", 0, float("inf")),
        lambda: VisualizationOperation.threshold("", 0, 1),
        lambda: VisualizationOperation.clip((float("nan"), 0, 1)),
        lambda: VisualizationOperation.clip((0, 0, 0)),
        lambda: VisualizationOperation.clip(origin=(0, float("inf"), 0)),
        lambda: VisualizationOperation.clip(invert="yes"),
        lambda: VisualizationOperation.glyph("U", float("inf")),
        lambda: VisualizationOperation.slice((0, 0, float("nan"))),
    ],
)
def test_new_filters_fail_closed_on_invalid_physics(operation):
    with pytest.raises(ValueError):
        operation()


def test_visualization_pipeline_does_not_accept_unknown_operation_types():
    with pytest.raises(TypeError, match="VisualizationOperation"):
        VisualizationPipeline().append("nonexistent")


def test_invalid_filter_cannot_mutate_original_pipeline_or_source():
    source = FakeDataset()
    good = VisualizationPipeline().append(VisualizationOperation.clip())
    bad = good.append(VisualizationOperation.threshold("missing", 0, 1))
    with pytest.raises(KeyError, match="missing"):
        bad.evaluate(source)
    assert source.history == ()
    assert good.evaluate(source).history[0][0] == "clip"
