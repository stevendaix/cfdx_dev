import subprocess
import sys


def test_renderer_contract_import_does_not_load_optional_gui_stack() -> None:
    code = """
import sys
from cfdx.renderer import NullRenderer
assert "PySide6" not in sys.modules
assert "pyvista" not in sys.modules
assert "pyvistaqt" not in sys.modules
NullRenderer().load("result.vtu")
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr


def test_renderer_contract_symbols_are_headless() -> None:
    from cfdx.renderer import NullRenderer, RenderObject, Renderer

    assert Renderer is not None
    assert RenderObject is not None
    assert NullRenderer is not None
