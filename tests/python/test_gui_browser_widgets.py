import importlib.util
import pytest

@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="optional Qt stack")
def test_browser_widgets_import():
    from cfdx.mesh_browser_panel import MeshBrowserPanel
    from cfdx.results_series_panel import ResultsSeriesPanel
    assert MeshBrowserPanel is not None and ResultsSeriesPanel is not None