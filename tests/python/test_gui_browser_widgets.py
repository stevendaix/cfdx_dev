import importlib.util
from pathlib import Path

import pytest


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="optional Qt stack"
)
def test_browser_widgets_import():
    from cfdx.mesh_browser_panel import MeshBrowserPanel
    from cfdx.results_series_panel import ResultsSeriesPanel
    assert MeshBrowserPanel is not None and ResultsSeriesPanel is not None


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="optional Qt stack"
)
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


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="optional Qt stack"
)
def test_3d_patch_actor_id_uses_stable_identity():
    from cfdx.gui_3d import PyVistaQtView

    assert PyVistaQtView.patch_actor_id("patch:inlet") == "patch:inlet"
    with pytest.raises(ValueError, match="invalid patch stable ID"):
        PyVistaQtView.patch_actor_id("patch:")


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="optional Qt stack"
)
def test_mesh_browser_drills_into_regions_patches_and_faces():
    from cfdx import create_application
    from cfdx.mesh_browser_panel import MeshBrowserPanel
    from cfdx.mesh_model import MeshCatalog, MeshPatch, MeshRegion
    from cfdx.setup_model import MeshSelection

    app = create_application(["cfdx-browser-tree-test"])
    panel = MeshBrowserPanel()
    emitted = []
    panel.selection_changed.connect(emitted.append)
    inlet = MeshPatch(0, "inlet", "inlet", 2, "fluid")
    wall = MeshPatch(1, "wall", "wall", 1, "fluid")
    region = MeshRegion("fluid", 8, 19, 3, (inlet, wall), ("c0", "c1"))
    catalog = MeshCatalog(16, 19, 5, (inlet, wall), (region,))
    panel.set_catalog(catalog)
    assert panel.status.text() == "5 cells / 19 faces / 16 points"
    assert panel.patches.topLevelItemCount() == 1
    region_item = panel.patches.topLevelItem(0)
    assert region_item.text(0) == "fluid  [3 cells]"
    panel.patches.setCurrentItem(region_item)
    assert emitted[-1] == MeshSelection("region", 0, "fluid", "region:fluid")
    assert region_item.childCount() == 2
    patch_item = region_item.child(0)
    assert patch_item.text(0) == "inlet  [inlet]  faces=2"
    panel.patches.setCurrentItem(patch_item)
    assert emitted[-1].kind == "patch"
    assert emitted[-1].stable_id == "region:fluid/patch:inlet"
    patch_item.setExpanded(True)
    assert patch_item.childCount() == 2
    assert patch_item.child(0).text(0) == "Face 0"
    assert patch_item.child(1).text(0) == "Face 1"
    panel.patches.setCurrentItem(patch_item.child(1))
    face_selection = emitted[-1]
    assert face_selection == MeshSelection(
        "face", 1, "inlet", "region:fluid/patch:inlet/face:1"
    )
    panel.clear()
    assert panel.status.text() == "No mesh loaded"
    assert panel.patches.topLevelItemCount() == 0
    panel.deleteLater()
    app.quit()


@pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="optional Qt stack"
)
def test_mesh_browser_keeps_flat_patch_list_for_default_region():
    from cfdx import create_application
    from cfdx.mesh_browser_panel import MeshBrowserPanel
    from cfdx.mesh_model import MeshCatalog, MeshPatch

    app = create_application(["cfdx-browser-flat-test"])
    panel = MeshBrowserPanel()
    emitted = []
    panel.selection_changed.connect(emitted.append)
    inlet = MeshPatch(0, "inlet", "inlet", 2)
    catalog = MeshCatalog(16, 19, 5, (inlet,))
    panel.set_catalog(catalog)
    assert panel.patches.topLevelItemCount() == 1
    patch_item = panel.patches.topLevelItem(0)
    assert patch_item.text(0) == "inlet  [inlet]  faces=2"
    panel.patches.setCurrentItem(patch_item)
    assert emitted[-1].kind == "patch"
    assert emitted[-1].stable_id == "patch:inlet"
    patch_item.setExpanded(True)
    assert patch_item.childCount() == 2
    panel.patches.setCurrentItem(patch_item.child(0))
    assert emitted[-1].kind == "face"
    assert emitted[-1].stable_id == "patch:inlet/face:0"
    panel.clear()
    assert panel.status.text() == "No mesh loaded"
    assert panel.patches.topLevelItemCount() == 0
    panel.deleteLater()
    app.quit()
