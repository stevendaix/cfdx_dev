import pytest

PySide6 = pytest.importorskip("PySide6", reason="GUI tests require the optional [gui] dependencies")

from cfdx import CFDXSession


def test_setup_panel_exposes_physics_controls_and_applies_model():
    from PySide6.QtWidgets import QApplication
    from cfdx.setup_panel import CaseSetupPanel

    app = QApplication.instance() or QApplication([])
    session = CFDXSession()
    panel = CaseSetupPanel(session.case)
    panel.physics_model.setCurrentText("incompressible")
    panel.physics_enabled.setChecked(True)
    panel.physics_fields["density"].setValue(1.225)
    panel.physics_fields["kinematic_viscosity"].setValue(1.5e-5)
    panel.physics_button.click()

    assert session.case.physics["incompressible"]["enabled"] is True
    assert session.case.physics["incompressible"]["density"] == pytest.approx(1.225)
    assert session.case.physics["incompressible"]["kinematic_viscosity"] == pytest.approx(1.5e-5)
    panel.close()
    app.processEvents()


def test_setup_panel_turbulence_filters_to_solver_ready_models():
    from PySide6.QtWidgets import QApplication
    from cfdx.setup_panel import CaseSetupPanel

    app = QApplication.instance() or QApplication([])
    session = CFDXSession()
    panel = CaseSetupPanel(session.case)
    panel.physics_model.setCurrentText("turbulence")
    models = [panel.physics_fields["model"].itemText(i) for i in range(panel.physics_fields["model"].count())]

    assert models
    assert "SST" in models
    assert "SMAGORINSKY" not in models
    panel.close()
    app.processEvents()


def test_setup_panel_routes_physics_diagnostic_to_field_widget():
    from PySide6.QtWidgets import QApplication
    from cfdx.setup_model import SetupDiagnostic
    from cfdx.setup_panel import CaseSetupPanel

    app = QApplication.instance() or QApplication([])
    panel = CaseSetupPanel(CFDXSession().case)
    panel.set_diagnostics(
        [SetupDiagnostic("error", "PHYSICS_VALUE", "density is invalid", "physics.incompressible.density")]
    )

    assert "1 error(s)" in panel.validation_summary.text()
    assert "density is invalid" in panel.physics_fields["density"].toolTip()
    panel.close()
    app.processEvents()


def test_setup_panel_selects_tab_for_material_boundary_and_initialization_diagnostics():
    from PySide6.QtWidgets import QApplication
    from cfdx.setup_model import SetupDiagnostic
    from cfdx.setup_panel import CaseSetupPanel

    app = QApplication.instance() or QApplication([])
    session = CFDXSession()
    session.case.materials["air"] = {"density": -1.0}
    session.case.boundaries["inlet"] = {"type": "inlet", "fields": ["pressure"]}
    panel = CaseSetupPanel(session.case)

    panel.set_diagnostics([SetupDiagnostic("error", "MATERIAL_VALUE", "density invalid", "materials.air.density")])
    assert panel.tabs.currentIndex() == 1
    assert "density invalid" in panel.material_density.toolTip()

    panel.set_diagnostics([SetupDiagnostic("error", "BOUNDARY_SCHEMA", "boundary invalid", "boundaries.inlet.type")])
    assert panel.tabs.currentIndex() == 2
    assert panel.boundary_list.currentItem().text() == "inlet"
    assert "boundary invalid" in panel.boundary_type.toolTip()

    panel.set_diagnostics([SetupDiagnostic("error", "INITIALIZATION", "initialization invalid", "physics.initialization")])
    assert panel.tabs.currentIndex() == 3
    assert "initialization invalid" in panel.initialization_mode.toolTip()
    panel.close()
    app.processEvents()
