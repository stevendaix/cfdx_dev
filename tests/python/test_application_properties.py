import pytest

from cfdx import CFDXSession
from cfdx.application import Application, WorkflowStatus
from cfdx.mesh_model import MeshCatalog, MeshPatch, MeshRegion
from cfdx.selection import face_selection_id, patch_selection_id
from cfdx.setup_model import ParameterType, SetupDiagnostic


def test_selection_resolves_properties_and_edit_updates_authoritative_case() -> None:
    session = CFDXSession()
    session.case.numerics["cfl"] = 1.0
    application = Application(session)

    application.select("numerics", "numerics", "Numerics")
    properties = application.properties()
    assert [prop.key for prop in properties] == ["numerics.cfl"]

    state = application.set_property("numerics.cfl", 2.5)
    assert session.case.numerics["cfl"] == 2.5
    assert state.project.dirty is True
    assert state.selection.stable_id == "numerics"


def test_validation_diagnostics_map_paths_to_workflow_status() -> None:
    session = CFDXSession()
    application = Application(session)
    diagnostic = SetupDiagnostic("error", "CFL", "CFL must be positive", "numerics.cfl")

    state = application.set_diagnostics([diagnostic])
    numerics = next(step for step in state.workflow if step.id == "numerics")

    assert numerics.status is WorkflowStatus.ERROR
    assert numerics.message == "CFL must be positive"


def test_solver_property_marks_restart_required() -> None:
    session = CFDXSession()
    application = Application(session)

    application.select("solver", "solver", "Solver")
    state = application.set_property("execution.solver", "/tmp/solver")

    assert session.case.execution.solver == "/tmp/solver"
    assert state.requires_restart is True


def test_mesh_selection_resolves_read_only_typed_properties() -> None:
    application = Application(CFDXSession())
    inlet = MeshPatch(0, "inlet", "inlet", 2, "fluid")
    region = MeshRegion("fluid", 16, 19, 5, (inlet,), ("c0", "c1"))
    application.mesh_catalog = MeshCatalog(16, 19, 5, (inlet,), (region,))

    application.select("region:fluid", "region", "fluid")
    properties = application.properties()
    assert [prop.key for prop in properties] == [
        "mesh.object",
        "mesh.cells",
        "mesh.faces",
        "mesh.points",
        "mesh.patches",
    ]
    cells = properties[1]
    assert cells.value == 5
    assert cells.kind is ParameterType.INTEGER
    assert cells.editable is False

    application.select("region:fluid/patch:inlet", "patch", "inlet")
    properties = application.properties()
    assert [prop.key for prop in properties] == [
        "mesh.object",
        "mesh.patch_type",
        "mesh.faces",
        "mesh.region",
    ]
    assert properties[1].value == "inlet"
    assert properties[3].value == "fluid"

    face_id = face_selection_id(patch_selection_id("inlet", "fluid"), 1)
    application.select(face_id, "face", "inlet")
    properties = application.properties()
    assert [prop.key for prop in properties] == [
        "mesh.object",
        "mesh.face",
        "mesh.patch",
        "mesh.region",
    ]
    assert properties[1].value == 1


def test_mesh_selection_properties_fail_closed_on_unknown_objects() -> None:
    application = Application(CFDXSession())
    inlet = MeshPatch(0, "inlet", "inlet", 2, "fluid")
    region = MeshRegion("fluid", 16, 19, 5, (inlet,), ("c0", "c1"))
    application.mesh_catalog = MeshCatalog(16, 19, 5, (inlet,), (region,))

    application.select("region:air", "region", "air")
    with pytest.raises(KeyError, match="mesh region not found"):
        application.properties()

    face_id = face_selection_id(patch_selection_id("inlet", "fluid"), 2)
    application.select(face_id, "face", "inlet")
    with pytest.raises(ValueError, match="mesh face out of range"):
        application.properties()


def test_mesh_selection_without_catalog_has_no_properties() -> None:
    application = Application(CFDXSession())

    application.select("region:fluid", "region", "fluid")

    assert application.properties() == ()


def test_set_property_rejects_untyped_values() -> None:
    session = CFDXSession()
    application = Application(session)

    with pytest.raises(ValueError, match="expected finite real"):
        application.set_property("numerics.cfl", "0.5")
    with pytest.raises(ValueError, match="physics enabled must be boolean"):
        application.set_property("physics.energy.enabled", "off")
