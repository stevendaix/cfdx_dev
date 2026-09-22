from cfdx.renderer import NullRenderer, RenderObject


def test_null_renderer_is_headless_and_selectable() -> None:
    renderer = NullRenderer()
    assert renderer.load("result.vtu") == (RenderObject("dataset", "result.vtu"),)
    renderer.select("dataset")
    assert renderer.selected == "dataset"


def test_null_renderer_rejects_unknown_selection() -> None:
    renderer = NullRenderer()
    renderer.load("result.vtu")
    try:
        renderer.select("missing")
    except KeyError:
        pass
    else:
        raise AssertionError("unknown object must be rejected")
