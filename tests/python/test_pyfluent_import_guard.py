"""The optional PyFluent import must not prevent ASCII adapter import."""

import builtins
import warnings

from cfdx.io.adapters.fluent import _try_pyfluent


def test_broken_optional_reader_does_not_break_ascii_adapter(monkeypatch):
    original_import = builtins.__import__

    def incompatible_import(name, *args, **kwargs):
        if name == "ansys.fluent.core.filereader.case_file":
            raise AttributeError("VariableCatalog.fluent")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", incompatible_import)
    with warnings.catch_warnings(record=True) as recorded:
        warnings.simplefilter("always")
        assert not _try_pyfluent()
    assert any("incompatible" in str(w.message) for w in recorded)
