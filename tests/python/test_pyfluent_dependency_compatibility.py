"""An installed, broken optional reader is not an acceptable CI skip."""
from importlib import metadata

import pytest


def test_installed_pyfluent_has_compatible_units_and_case_reader():
    try:
        version = metadata.version("ansys-fluent-core")
    except metadata.PackageNotFoundError:
        pytest.skip("PyFluent is not installed in this environment")

    try:
        from ansys.units import VariableCatalog

        if not hasattr(VariableCatalog, "fluent"):
            pytest.fail(
                "ansys-units is incompatible with PyFluent: "
                "VariableCatalog.fluent is missing",
                pytrace=False,
            )
        from ansys.fluent.core.filereader.case_file import CaseFile
    except (ImportError, AttributeError) as exc:
        pytest.fail(
            f"ansys-fluent-core {version} installed but unusable: {exc}",
            pytrace=False,
        )
    assert CaseFile is not None
