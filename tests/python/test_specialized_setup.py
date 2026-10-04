import pytest

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
