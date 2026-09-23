from __future__ import annotations

import math
from pathlib import Path
import sys
import time
import xml.etree.ElementTree as ET

import h5py

from cfdx import CFDXSession, ExecutionController, SolverRunner
from cfdx.case_io import read_case, read_case_with_dat, save_case, save_case_with_dat, validate_case_bundle
from cfdx.dat_io import read_dat_restart
from cfdx.mesh_model import read_mesh_catalog
from cfdx.results_series import discover_result_series
from cfdx.validation import validate_case


def _write_unit_cube_mesh(path: Path) -> None:
    with h5py.File(path, "w") as h5:
        h5.attrs["format"] = "CFDX"
        h5.attrs["schema_version"] = 1
        h5.create_dataset("points", data=[
            [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
            [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.]])
        h5.create_dataset("face_vertices", data=[
            0,3,2,1, 4,5,6,7, 0,1,5,4, 3,7,6,2, 0,4,7,3, 1,2,6,5], dtype="u8")
        h5.create_dataset("face_offsets", data=[0,4,8,12,16,20,24], dtype="u8")
        h5.create_dataset("owner", data=[0]*6, dtype="u8")
        h5.create_dataset("neighbour", data=[-1]*6, dtype="i8")
        h5.create_dataset("cell_faces", data=[0,1,2,3,4,5], dtype="u8")
        h5.create_dataset("cell_offsets", data=[0,6], dtype="u8")
        h5.attrs["boundary_patches"] = "wall:0:6:0"
        h5.create_dataset("patch_face_ids", data=[0,1,2,3,4,5], dtype="u8")
        h5.create_dataset("patch_face_offsets", data=[0,6], dtype="u8")
        fields = h5.create_group("fields")
        fields.create_dataset("values", data=[0.])
        fields.attrs["name"] = "p"
        fields.attrs["unit"] = "Pa"
        fields.attrs["dimension"] = "1"


def _make_session(driver: Path) -> CFDXSession:
    session = CFDXSession()
    session.case.name = "issue168-e2e"
    session.case.physics["incompressible"] = {
        "enabled": True, "density": 1.0,
        "kinematic_viscosity": 1.0e-3, "algorithm": "SIMPLE"}
    session.case.physics["initialization"] = {
        "mode": "uniform", "field": "U", "value": 0.0}
    session.case.numerics.update({"cfl": 1.0, "max_iterations": 2})
    session.case.materials["water"] = {
        "density": 1000.0, "dynamic_viscosity": 1.0e-3,
        "cp": 4180.0, "conductivity": 0.6}
    session.case.boundaries["wall"] = {
        "type": "wall", "fields": ["velocity"],
        "velocity_type": "FIXED_VALUE",
        "value_x": 0.0, "value_y": 0.0, "value_z": 0.0}
    session.case.execution.solver = str(driver)
    session.case.execution.mpi_ranks = 1
    session.case.execution.deterministic = True
    return session


def _wait_running(controller: ExecutionController) -> None:
    deadline = time.monotonic() + 5.0
    while not controller.runner.running and time.monotonic() < deadline:
        time.sleep(0.01)
    assert controller.runner.running


def _assert_vtu_frame(path: Path, expected_iteration: int) -> None:
    root = ET.parse(path).getroot()
    field_data = root.find("./UnstructuredGrid/FieldData")
    assert field_data is not None
    metadata = {item.attrib["Name"]: float((item.text or "").strip())
                for item in field_data.findall("DataArray")}
    assert metadata["iteration"] == expected_iteration
    assert math.isfinite(metadata["physical_time"])
    piece = root.find("./UnstructuredGrid/Piece")
    assert piece is not None
    assert piece.attrib["NumberOfPoints"] == "8"
    assert piece.attrib["NumberOfCells"] == "1"
    types = piece.find("./Cells/DataArray[@Name='types']")
    assert types is not None
    assert (types.text or "").strip() == "42"
    cell_data = piece.find("./CellData")
    assert cell_data is not None
    names = {item.attrib["Name"] for item in cell_data.findall("DataArray")}
    assert {"Ux", "Uy", "Uz", "p"} <= names


def _assert_field_equivalent(reference: Path, candidate: Path) -> None:
    ref = read_dat_restart(reference)
    got = read_dat_restart(candidate)
    assert ref.cells == got.cells == 1
    assert ref.fields.keys() == got.fields.keys()
    for name in ref.fields:
        assert ref.fields[name].dimension == got.fields[name].dimension
        assert got.fields[name].values == ref.fields[name].values


def main(driver: str) -> int:
    driver_path = Path(driver).resolve()
    assert driver_path.is_file()
    root = Path.cwd() / "cfdx_issue168_e2e"
    root.mkdir(parents=True, exist_ok=True)
    case_path = root / "issue168.cfdx.h5"
    run1 = root / "run1"
    run2 = root / "run2"
    for directory in (run1, run2):
        directory.mkdir(exist_ok=True)
        for child in directory.iterdir():
            if child.is_file():
                child.unlink()

    _write_unit_cube_mesh(case_path)
    session = _make_session(driver)
    catalog = read_mesh_catalog(case_path)
    assert catalog.n_cells == 1 and catalog.patch_names == ("wall",)
    report = validate_case(session.case, catalog)
    assert report.ok, [d.message for d in report.diagnostics]
    session.validate()
    assert session.state.value == "READY"
    save_case(session, case_path)
    assert validate_case_bundle(case_path)["mesh"]

    runner = SolverRunner(
        [str(driver_path), str(case_path), str(run1),
         "--iterations", "2", "--delay-ms", "100"], cwd=root)
    controller = ExecutionController(session, runner)
    controller.start()
    _wait_running(controller)
    controller.pause()
    assert session.state.value == "PAUSED" and runner.paused
    controller.resume()
    assert session.state.value == "RUNNING" and not runner.paused
    assert runner._thread is not None
    runner._thread.join(timeout=15)
    assert session.state.value == "CONVERGED"
    assert session.iteration == 2
    assert len(controller.monitor_series.samples) == 2

    step1, step2 = run1 / "step_0001.dat", run1 / "step_0002.dat"
    assert step1.is_file() and step2.is_file()
    checkpoint = read_dat_restart(step1)
    assert checkpoint.iteration == 1 and checkpoint.time == 0.0

    series = discover_result_series(run1, inspect_fields=True, require_physical_time=True)
    assert [f.path.name for f in series.frames] == ["step_0001.vtu", "step_0002.vtu"]
    assert series.frames[-1].iteration == 2
    _assert_vtu_frame(series.frames[-1].path, 2)

    case_path, paired_dat = save_case_with_dat(session, case_path, step1)
    loaded_case = read_case(case_path)
    assert loaded_case.iteration == session.iteration
    loaded_session, loaded_dat = read_case_with_dat(case_path, paired_dat)
    assert loaded_dat == paired_dat
    assert read_dat_restart(paired_dat).iteration == 1
    assert loaded_session.case.name == "issue168-e2e"

    restart_runner = SolverRunner(
        [str(driver_path), str(case_path), str(run2),
         "--restart", str(step1), "--iterations", "1"], cwd=root)
    restart_controller = ExecutionController(loaded_session, restart_runner)
    restart_controller.start()
    assert restart_runner._thread is not None
    restart_runner._thread.join(timeout=15)
    assert loaded_session.state.value == "CONVERGED"
    restart_step = run2 / "step_0001.dat"
    assert restart_step.is_file()
    _assert_field_equivalent(step2, restart_step)
    _assert_vtu_frame(run2 / "step_0001.vtu", 1)

    result_only = read_case(case_path)
    assert result_only.case.name == "issue168-e2e"
    assert paired_dat.is_file()
    print("CFDX issue #168 E2E: mesh/setup/check/initialize/run/pause/monitor/"
          "DAT save/reload/restart/time-series/3D inspection PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1]))
