"""Quantitative headless application E2E for the GUI workflow.

The test deliberately exercises the real application orchestration stack:
mesh catalog import, case validation, initialization state, asynchronous solver
execution, live monitoring, pause/resume, DAT persistence/reload, restart and
VTU result-series inspection. The solver is a deterministic executable fixture
so CI does not depend on an external installation.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from textwrap import dedent
from xml.etree import ElementTree

import h5py
import pytest

from cfdx.case_io import read_case_with_dat, save_case_with_dat
from cfdx.dat_io import DatField, DatRestart, write_dat_hdf5
from cfdx.execution import ExecutionController
from cfdx.mesh_model import read_mesh_catalog
from cfdx.results_series import discover_result_series
from cfdx.session import CFDXSession, SimulationState
from cfdx.validation import validate_case
from cfdx.runner import SolverRunner


def _write_mesh(path: Path) -> None:
    """Create the smallest real CFDX mesh/catalog accepted by the application."""
    with h5py.File(path, "w") as h5:
        h5.create_dataset("points", data=[
            0.0, 0.0, 0.0,
            1.0, 0.0, 0.0,
            1.0, 1.0, 0.0,
            0.0, 1.0, 0.0,
            0.0, 0.0, 1.0,
            1.0, 0.0, 1.0,
            1.0, 1.0, 1.0,
            0.0, 1.0, 1.0,
        ])
        # One quad face; the application catalog only needs consistent topology
        # and patch identity for this orchestration acceptance test.
        h5.create_dataset("face_vertices", data=[0, 1, 2, 3])
        h5.create_dataset("face_offsets", data=[0, 4])
        h5.create_dataset("owner", data=[0])
        h5.create_dataset("neighbour", data=[-1])
        h5.create_dataset("cell_faces", data=[0])
        h5.create_dataset("cell_offsets", data=[0, 1])
        h5.create_dataset("patch_face_ids", data=[0])
        h5.create_dataset("patch_face_offsets", data=[0, 1])
        h5.create_dataset("cell_ids", data=["0"])
        h5.attrs["n_cells"] = 1
        h5.attrs["boundary_patches"] = "inlet:inlet:0:1"


def _write_solver(path: Path, results: Path) -> None:
    path.write_text(
        dedent(
            f"""
            import sys
            import time
            from pathlib import Path

            results = Path(r"{results}")
            results.mkdir(parents=True, exist_ok=True)
            restart = "--restart" in sys.argv
            start = 3 if restart else 1
            for iteration in range(start, start + 2):
                physical_time = 0.1 * iteration
                (results / f"result_{{iteration:04d}}.vtu").write_text(
                    f'''<VTKFile type="UnstructuredGrid" version="0.1">
                    <UnstructuredGrid>
                      <FieldData>
                        <DataArray type="Float64" Name="physical_time" NumberOfTuples="1">\
{{physical_time}}</DataArray>
                        <DataArray type="Int64" Name="iteration" NumberOfTuples="1">\
{{iteration}}</DataArray>
                      </FieldData>
                      <Piece NumberOfPoints="0" NumberOfCells="0"/>
                    </UnstructuredGrid>
                    </VTKFile>''',
                    encoding="utf-8",
                )
                print(f"Iteration {{iteration}} Time = {{physical_time}} CFL: 0.5", flush=True)
                if not restart and iteration == 1:
                    (results / "ready").write_text("ready", encoding="utf-8")
                    time.sleep(2.0)
            """,
        ),
        encoding="utf-8",
    )


def _wait_for(path: Path, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if path.exists():
            return
        time.sleep(0.02)
    raise AssertionError(f"timeout waiting for {path}")


def _wait_for_state(session: CFDXSession, state: SimulationState, timeout: float = 5.0) -> None:
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if session.state is state:
            return
        time.sleep(0.02)
    raise AssertionError(f"timeout waiting for state {state}; got {session.state}")


def test_full_application_workflow_quantitative(tmp_path: Path) -> None:
    mesh_path = tmp_path / "case.cfdx.h5"
    solver_path = tmp_path / "solver.py"
    results_dir = tmp_path / "results"
    _write_mesh(mesh_path)
    _write_solver(solver_path, results_dir)

    # 1. Mesh import through the same application catalog used by the GUI.
    catalog = read_mesh_catalog(mesh_path)
    assert catalog.n_cells == 1
    assert catalog.patch_names == ("inlet",)

    # 2. Setup + initialization + case check.
    session = CFDXSession()
    session.case.name = "e2e-channel"
    session.case.enable("incompressible")
    session.case.materials["air"] = {
        "density": 1.2,
        "dynamic_viscosity": 1.8e-5,
        "cp": 1005.0,
        "conductivity": 0.026,
    }
    session.case.set_numerics(cfl=0.5)
    session.case.set_boundary("inlet", type="inlet", value="1.0")
    session.case.physics["initialization"] = {"mode": "uniform", "value": 0.0}
    session.case.execution.solver = sys.executable
    session.case.execution.restart_option = "--restart"
    report = validate_case(session.case, catalog)
    assert report.ok, report.diagnostics
    assert session.case.physics["initialization"]["mode"] == "uniform"

    # 3. Run + monitor + pause/resume using the real asynchronous runner.
    runner = SolverRunner([sys.executable, str(solver_path), str(mesh_path)])
    controller = ExecutionController(session, runner)
    controller.start()
    _wait_for(results_dir / "ready")
    _wait_for_state(session, SimulationState.RUNNING)
    _wait_for_state(session, SimulationState.RUNNING)

    # Give the first metric line time to cross the stdout reader.
    deadline = time.monotonic() + 5.0
    while time.monotonic() < deadline and session.iteration < 1:
        time.sleep(0.02)
    assert session.iteration == 1
    assert session.time == pytest.approx(0.1)
    assert controller.monitor_series.samples
    assert controller.monitor_series.samples[-1].iteration == 1

    controller.pause()
    assert session.state is SimulationState.PAUSED
    assert runner.paused

    # 4. Save a numerical DAT artifact while paused.
    source_dat = tmp_path / "solver.dat"
    write_dat_hdf5(
        source_dat,
        DatRestart(
            version=2,
            cells=1,
            iteration=session.iteration,
            time=session.time,
            fields={
                "U": DatField("U", 3, [1.0, 0.0, 0.0]),
                "p": DatField("p", 1, [101325.0]),
            },
            cell_ids=(0,),
        ),
    )
    case_path, dat_path = save_case_with_dat(
        session, mesh_path, source_dat
    )
    assert dat_path.is_file()

    controller.resume()
    assert session.state is SimulationState.RUNNING
    runner._thread.join(timeout=5.0)
    assert session.state is SimulationState.CONVERGED
    assert session.iteration == 2
    assert session.time == pytest.approx(0.2)

    # 5. Reload case + DAT and verify numerical state/metadata.
    reloaded, reloaded_dat = read_case_with_dat(case_path)
    assert reloaded.case.name == "e2e-channel"
    assert reloaded.iteration == 1
    assert reloaded.time == pytest.approx(0.1)
    assert reloaded_dat == dat_path

    # 6. Restart from the persisted DAT through the real execution controller.
    restart_session = reloaded
    restart_runner = SolverRunner([sys.executable, str(solver_path), str(mesh_path)])
    restart_controller = ExecutionController(restart_session, restart_runner)
    restart_controller.restart(reloaded_dat)
    restart_runner._thread.join(timeout=5.0)
    assert restart_session.state is SimulationState.CONVERGED
    assert restart_session.iteration == 4
    assert restart_session.time == pytest.approx(0.4)
    assert restart_controller.monitor_series.samples[-1].iteration == 4

    # 7. Discover the complete result series and inspect authoritative metadata.
    series = discover_result_series(results_dir)
    assert len(series.frames) == 4
    assert [frame.sequence for frame in series.frames] == [0, 1, 2, 3]

    metadata = []
    for frame in series.frames:
        root = ElementTree.parse(frame.path).getroot()
        values = {
            item.attrib["Name"]: item.text
            for item in root.iter("DataArray")
            if item.attrib.get("Name") in {"physical_time", "iteration"}
        }
        metadata.append((int(values["iteration"]), float(values["physical_time"])))

    assert metadata == [
        (1, pytest.approx(0.1)),
        (2, pytest.approx(0.2)),
        (3, pytest.approx(0.3)),
        (4, pytest.approx(0.4)),
    ]
