import importlib.util
from pathlib import Path

import pytest


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="optional Qt stack")
def test_browser_widgets_import():
    from cfdx.mesh_browser_panel import MeshBrowserPanel
    from cfdx.results_series_panel import ResultsSeriesPanel
    assert MeshBrowserPanel is not None and ResultsSeriesPanel is not None


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="optional Qt stack")
def test_results_series_panel_exposes_field_and_playback_controls():
    from cfdx import create_application
    from cfdx.results_series import ResultFrame, ResultSeries
    from cfdx.results_series_panel import ResultsSeriesPanel

    app = create_application(["cfdx-results-panel-test"])
    panel = ResultsSeriesPanel()
    series = ResultSeries((
        ResultFrame(Path("step_0.vtu"), 0, 0.0, True, ("U", "p")),
        ResultFrame(Path("step_1.vtu"), 1, 1.0, True, ("U", "p")),
    ))
    fields = []
    frames = []
    panel.field_changed.connect(fields.append)
    panel.frame_changed.connect(frames.append)
    panel.set_series(series)
    assert [panel.field.itemText(i) for i in range(panel.field.count())] == ["U", "p"]
    assert panel.slider.maximum() == 1
    panel.field.setCurrentText("p")
    panel.slider.setValue(1)
    assert fields == ["p"]
    assert frames[-1].path.name == "step_1.vtu"
    panel.timer.stop()
    panel.deleteLater()
    app.quit()


@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="optional Qt stack")
def test_3d_patch_actor_id_uses_stable_identity():
    from cfdx.gui_3d import PyVistaQtView

    assert PyVistaQtView.patch_actor_id("patch:inlet") == "patch:inlet"
    with pytest.raises(ValueError, match="invalid patch stable ID"):
        PyVistaQtView.patch_actor_id("patch:")
