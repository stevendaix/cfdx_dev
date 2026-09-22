import importlib.util

import pytest


@pytest.mark.skipif(importlib.util.find_spec("pyvista") is None, reason="PyVista optional")
def test_pyvista_renderer_loads_vtu(tmp_path) -> None:
    import pyvista as pv

    from cfdx.renderer_pyvista import PyVistaRenderer

    mesh = pv.Sphere()
    path = tmp_path / "sphere.vtu"
    mesh.save(path)
    renderer = PyVistaRenderer()
    objects = renderer.load(str(path))
    assert objects[0].object_id == "dataset"
    renderer.select("dataset")
    assert renderer.slice().n_points > 0
    renderer.clear()
    assert renderer.dataset is None


def test_pyvista_renderer_has_clear_dependency_error() -> None:
    if importlib.util.find_spec("pyvista") is not None:
        pytest.skip("PyVista is installed")
    from cfdx.renderer_pyvista import PyVistaRenderer

    with pytest.raises(RuntimeError, match="PyVista is required"):
        PyVistaRenderer()
