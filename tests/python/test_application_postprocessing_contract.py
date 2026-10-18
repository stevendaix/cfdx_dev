"""Post-processing contract qualification against real production output.

The registry declares capabilities; the Workbench must never present a
field as selectable when the production result series cannot actually
produce it. This module closes that loop by asserting the registry's
availability claims against the real VTU frames emitted by the production
solver, so a regression that re-advertises an unavailable field is caught
by the test suite rather than by a user noticing a greyed-out control.
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import h5py
import pytest
from xml.etree import ElementTree

from cfdx.application.postprocessing import build_post_processing_registry
from cfdx.application.results import ResultDataset, ResultsState
from cfdx.application.workflow import build_workflow_state
from cfdx.session import CFDXSession, SimulationState


SOLVER_ENV = "CFDX_PRODUCTION_SOLVER"
MESH_ENV = "CFDX_PRODUCTION_MESH"


def _production_solver() -> str:
    solver = os.environ.get(SOLVER_ENV)
    if not solver:
        pytest.skip(f"{SOLVER_ENV} is provided by ctest")
    return solver


def _production_mesh() -> str:
    mesh = os.environ.get(MESH_ENV)
    if not mesh:
        pytest.skip(f"{MESH_ENV} is provided by ctest")
    return mesh


def _run_solver(solver: str, mesh: str, output_dir: Path,
                iterations: int = 2) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [solver, "--mesh", mesh, "--output-dir", str(output_dir),
         "--iterations", str(iterations)],
        capture_output=True, text=True, timeout=60, check=False,
    )


def _vtu_fields(vtu_path: Path) -> tuple[str, ...]:
    """Return the cell-data field names declared in a VTU frame."""
    tree = ElementTree.parse(vtu_path)
    root = tree.getroot()
    # CellData and PointData both declare arrays; the authoritative set for
    # field display is the union, deduplicated, excluding VTK-internal
    # connectivity/offsets/types bookkeeping.
    names: list[str] = []
    for container in root.iter():
        tag = container.tag.split("}")[-1] if "}" in container.tag else container.tag
        if tag in {"CellData", "PointData"}:
            for arr in container.findall("DataArray"):
                name = arr.get("Name")
                if name and name not in names and not name.startswith("connectivity"):
                    names.append(name)
    return tuple(names)


def test_registry_availability_matches_real_vtu_output(tmp_path: Path) -> None:
    """A spec marked available must actually be present in the VTU frames."""
    solver = _production_solver()
    mesh = _production_mesh()

    output_dir = tmp_path / "run"
    completed = _run_solver(solver, mesh, output_dir, iterations=2)
    assert completed.returncode in (0, 1), completed.stderr

    vtus = sorted(output_dir.glob("result_*.vtu"))
    assert vtus, "production solver must emit at least one VTU frame"
    actual_fields = set()
    for vtu in vtus:
        actual_fields.update(_vtu_fields(vtu))

    registry = build_post_processing_registry()
    available_specs = [s for s in registry.all() if s.availability == "available"]

    # Every "available" native field must be present in at least one frame.
    native_fields = {s.stable_id.removeprefix("field:") for s in available_specs
                     if s.category == "native_field"}
    for field in native_fields:
        assert field in actual_fields, (
            f"registry marks field:{field} as available but the production "
            f"solver's VTU frames only contain {sorted(actual_fields)}"
        )


def test_conditional_fields_are_not_selectable_from_vtu(tmp_path: Path) -> None:
    """A conditional field must not be selectable from the result series."""
    solver = _production_solver()
    mesh = _production_mesh()

    output_dir = tmp_path / "run"
    completed = _run_solver(solver, mesh, output_dir, iterations=2)
    assert completed.returncode in (0, 1), completed.stderr

    vtus = sorted(output_dir.glob("result_*.vtu"))
    assert vtus
    actual_fields = set()
    for vtu in vtus:
        actual_fields.update(_vtu_fields(vtu))

    # Build a ResultDataset from the real frame and prove U is not selectable.
    fields = _vtu_fields(vtus[0])
    dataset = ResultDataset(
        stable_id="frame_1",
        frame_id="frame_1",
        path=str(vtus[0]),
        fields=fields,
    )
    assert dataset.field_available("p") is True
    assert dataset.field_available("U") is False, (
        "field:U is conditional in the registry and must not be selectable "
        "from the production result series, but the dataset reports it present"
    )


def test_conditional_fields_render_as_conditional_in_workflow(tmp_path: Path) -> None:
    """The Workbench must not show a green check for conditional capabilities."""
    solver = _production_solver()
    mesh = _production_mesh()

    output_dir = tmp_path / "run"
    completed = _run_solver(solver, mesh, output_dir, iterations=2)
    assert completed.returncode in (0, 1), completed.stderr

    session = CFDXSession()
    session.case.name = "workflow-qualification"
    session.case.enable("incompressible")
    session.case.materials["air"] = {
        "density": 1.2, "dynamic_viscosity": 1.8e-5,
        "cp": 1005.0, "conductivity": 0.026,
    }
    session.case.execution.solver = solver
    session.case.execution.mesh = mesh
    session.case.execution.restart_option = "--restart"
    session.state = SimulationState.STOPPED

    state = build_workflow_state(session)
    by_id = {step.id: step for step in state}

    # The RESULTS group must be CONDITIONAL because field:U and the
    # mesh-quality metrics are conditional: the group cannot be reported
    # complete even though every declared node is present.
    assert by_id["workflow.results"].status.value == "conditional"
    assert by_id["fields"].status.value == "conditional"

    # A conditional node must never be reported complete.
    assert by_id["fields"].status.value != "complete"