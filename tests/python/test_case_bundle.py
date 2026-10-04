from __future__ import annotations

import json
from pathlib import Path

import h5py
import pytest

from cfdx import CFDXSession
from cfdx.case_io import read_case, save_case, validate_case_bundle


def make_session() -> CFDXSession:
    session = CFDXSession()
    session.case.name = "bundle"
    session.case.physics["model"] = "incompressible"
    session.case.set_numerics(cfl=1.0, scheme="linear")
    session.case_revision = 3
    session.mesh_revision = 2
    session.physics_revision = 4
    session.numerics_revision = 5
    return session


def write_native_mesh(h5: h5py.File) -> None:
    h5.attrs["format"] = "CFDX"
    h5.attrs["schema_version"] = 1
    h5.create_dataset("points", data=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    h5.create_dataset("face_vertices", data=[0, 2, 1, 0, 1, 3, 1, 3, 2, 2, 3, 0], dtype="u8")
    h5.create_dataset("face_offsets", data=[0, 3, 6, 9, 12], dtype="u8")
    h5.create_dataset("owner", data=[0, 0, 0, 0], dtype="u8")
    h5.create_dataset("neighbour", data=[-1, -1, -1, -1], dtype="i8")
    h5.create_dataset("cell_faces", data=[0, 1, 2, 3], dtype="u8")
    h5.create_dataset("cell_offsets", data=[0, 4], dtype="u8")
    h5.create_group("fields").create_dataset("values", data=[101.0])
    h5["fields"].attrs["name"] = "p"
    h5["fields"].attrs["unit"] = "Pa"
    h5["fields"].attrs["dimension"] = "1"


def test_case_bundle_is_single_self_contained_hdf5_artifact(tmp_path: Path) -> None:
    path = tmp_path / "bundle.cfdx.h5"
    with h5py.File(path, "w") as h5:
        write_native_mesh(h5)

    save_case(make_session(), path)

    loaded = read_case(path)
    assert loaded.case.physics["model"] == "incompressible"
    assert loaded.case.numerics["cfl"] == 1.0

    report = validate_case_bundle(path)
    assert report["mesh"] is True
    assert report["fields"] is True
    assert report["physics"] is True
    assert report["numerics"] is True
    assert report["runtime_separated"] is True


def test_case_bundle_rejects_missing_mesh(tmp_path: Path) -> None:
    path = save_case(make_session(), tmp_path / "incomplete.cfdx.h5")
    with pytest.raises(ValueError, match="mesh"):
        validate_case_bundle(path)


def test_case_bundle_rejects_missing_field_data(tmp_path: Path) -> None:
    path = tmp_path / "incomplete.cfdx.h5"
    with h5py.File(path, "w") as h5:
        write_native_mesh(h5)
        del h5["fields"]
    save_case(make_session(), path)
    with pytest.raises(ValueError, match="fields"):
        validate_case_bundle(path)


def bundle_with_parasite(tmp_path: Path, populate) -> Path:
    """Write a valid bundle and then inject runtime state into the case artifact."""
    path = tmp_path / "leaky.cfdx.h5"
    with h5py.File(path, "w") as h5:
        write_native_mesh(h5)
    save_case(make_session(), path)
    with h5py.File(path, "a") as h5:
        populate(h5)
    return path


@pytest.mark.parametrize(
    "populate",
    [
        pytest.param(lambda h5: h5.create_group("runtime/checkpoint"), id="checkpoint"),
        pytest.param(lambda h5: h5.create_group("runtime/restart"), id="restart"),
        pytest.param(lambda h5: h5.create_dataset("runtime/foo", data=[1.0]), id="arbitrary-runtime-group"),
        pytest.param(lambda h5: h5.create_dataset("iteration", data=[7]), id="top-level-iteration"),
        pytest.param(lambda h5: h5.create_dataset("time", data=[0.5]), id="top-level-time"),
        pytest.param(lambda h5: h5.create_dataset("solution/u", data=[1.0]), id="top-level-solution"),
    ],
)
def test_runtime_state_parasites_are_reported_as_not_separated(tmp_path: Path, populate) -> None:
    """Absent checkpoint/restart groups alone never proved runtime separation.

    Each parasite is runtime state that leaked into the case artifact under a
    name other than the two historical ones, so it must flip the flag.
    """
    path = bundle_with_parasite(tmp_path, populate)
    assert validate_case_bundle(path)["runtime_separated"] is False


def test_case_layout_members_are_not_mistaken_for_runtime_state(tmp_path: Path) -> None:
    """The check must not fire on the mesh and field layout it has to allow."""
    path = tmp_path / "bundle.cfdx.h5"
    with h5py.File(path, "w") as h5:
        write_native_mesh(h5)
        h5.create_dataset("patch_face_ids", data=[0, 1], dtype="u8")
        h5.create_group("fields/scalar")
    save_case(make_session(), path)
    assert validate_case_bundle(path)["runtime_separated"] is True
