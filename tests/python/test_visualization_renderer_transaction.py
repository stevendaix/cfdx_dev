"""Renderer stage/commit tests independent of the optional Qt/VTK runtime."""

from __future__ import annotations

import pytest
from cfdx.gui_3d import PyVistaQtView
from cfdx.visualization_pipeline import VisualizationOperation, VisualizationPipeline


class Data:
    point_data = {"p": object()}
    cell_data = {}

    def __init__(self, history=()) -> None:
        self.history = history

    def clip(self, **kwargs):
        return Data((*self.history, ("clip", kwargs)))


class Plotter:
    def __init__(self):
        self.clears = 0
        self.calls = []
        self.fail_next_add = False

    def disable_picking(self):
        pass

    def clear(self):
        self.clears += 1

    def add_mesh(self, data, *, scalars=None):
        if self.fail_next_add:
            self.fail_next_add = False
            raise RuntimeError("injected VTK failure")
        self.calls.append((data, scalars))

    def render(self):
        pass


class ViewDouble:
    def __init__(self):
        self.data = Data()
        self.plotter = Plotter()
        self._visualization_pipeline = VisualizationPipeline()
        self._rendered_dataset = self.data
        self._rendered_scalars = None
        self._selected_field = None
        self._selection_highlight = None
        self._selection_id = None

    def _dataset(self):
        return self.data

    def _available_fields(self, data):
        return set(data.point_data) | set(data.cell_data)


def test_failed_filter_keeps_pipeline_and_existing_scene_unchanged():
    view = ViewDouble()
    candidate = VisualizationPipeline().append(
        VisualizationOperation.threshold("missing", 0, 1)
    )
    with pytest.raises(KeyError, match="missing"):
        PyVistaQtView._render_visualization(view, candidate)
    assert view.plotter.clears == 0
    assert view._visualization_pipeline.operations == ()
    assert view._rendered_dataset is view.data


def test_renderer_failure_restores_previous_view_and_pipeline_state():
    view = ViewDouble()
    candidate = VisualizationPipeline().append(VisualizationOperation.clip())
    view.plotter.fail_next_add = True
    with pytest.raises(RuntimeError, match="injected VTK failure"):
        PyVistaQtView._render_visualization(view, candidate)
    assert view.plotter.clears == 2
    assert view.plotter.calls[-1] == (view.data, None)
    assert view._rendered_dataset is view.data
    assert view._visualization_pipeline.operations == ()


def test_successful_filter_commits_pipeline_after_renderer_success():
    view = ViewDouble()
    candidate = VisualizationPipeline().append(VisualizationOperation.clip())
    PyVistaQtView._render_visualization(view, candidate)
    assert view._visualization_pipeline == candidate
    assert view._rendered_dataset is not view.data
    assert view._rendered_dataset.history[0][0] == "clip"
    assert view.plotter.clears == 1
