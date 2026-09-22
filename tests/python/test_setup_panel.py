import importlib.util
import pytest

@pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="optional Qt stack")
def test_professional_setup_panel_imports():
    from cfdx.setup_panel import CaseSetupPanel
    assert CaseSetupPanel is not None
