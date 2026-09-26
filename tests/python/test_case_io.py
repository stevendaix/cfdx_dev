from pathlib import Path

import h5py
import pytest

from cfdx import CFDXSession
from cfdx.case_io import read_case, read_case_with_dat, save_case, save_case_with_dat
from cfdx.dat_io import DatField, DatRestart, write_dat_hdf5


def make_session() -> CFDXSession:
    session = CFDXSession()
    session.case.name = "channel"
    session.case.enable("energy")
    session.case.materials["air"] = {"density": 1.2, "dynamic_viscosity": 1.8e-5, "cp": 1005.0, "conductivity": 0.026}
    session.case.set_numerics(cfl=2.5, scheme="upwind")
    session.case.set_boundary("inlet", type="inlet", value="1.0")
    session.case.execution.solver = "cfdx-solver"
    session.case.execution.mpi_ranks = 4
    session.case_revision = 7
    session.mesh_revision = 3
    session.physics_revision = 5
    session.numerics_revision = 6
    session.iteration = 120
    session.time = 2.5
    return session


def test_case_hdf5_roundtrip_preserves_configuration_without_numerical_state(tmp_path: Path) -> None:
    original = make_session()
    path = save_case(original, tmp_path / "channel.cfdx.h5")

    loaded = read_case(path)

    assert loaded.case.as_dict() == original.case.as_dict()
    assert loaded.case_revision == 7
    assert loaded.mesh_revision == 3
    assert loaded.physics_revision == 5
    assert loaded.numerics_revision == 6
    assert loaded.iteration == 0
    assert loaded.time == pytest.approx(0.0)


def test_case_hdf5_contains_no_separate_checkpoint_file(tmp_path: Path) -> None:
    path = save_case(make_session(), tmp_path / "channel.cfdx.h5")

    assert path.is_file()
    assert not path.with_suffix(".json").exists()
    with h5py.File(path, "r") as h5:
        assert "case/config" in h5
        assert "runtime/checkpoint" not in h5
        assert "runtime/restart" not in h5


def test_save_and_read_case_with_dat(tmp_path: Path) -> None:
    source = tmp_path / "solver.dat"
    write_dat_hdf5(
        source,
        DatRestart(
            version=2,
            cells=1,
            iteration=120,
            time=2.5,
            fields={"U": DatField("U", 3, [1.0, 2.0, 3.0]), "p": DatField("p", 1, [4.0])},
        ),
    )
    case_path, dat_path = save_case_with_dat(
        make_session(), tmp_path / "channel.cfdx.h5", source
    )

    assert dat_path.name == "channel.dat.h5"
    loaded, loaded_dat = read_case_with_dat(case_path)

    assert loaded.case.as_dict() == make_session().case.as_dict()
    assert loaded.iteration == 120
    assert loaded.time == pytest.approx(2.5)
    assert loaded_dat == dat_path
    assert loaded_dat.read_bytes() == source.read_bytes()


def test_read_case_with_dat_is_independent_of_case_state(tmp_path: Path) -> None:
    source = tmp_path / "solver.dat"
    write_dat_hdf5(source, DatRestart(2, 1, 120, 2.5, {"p": DatField("p", 1, [1.0])}))
    case_path, dat_path = save_case_with_dat(
        make_session(), tmp_path / "channel.cfdx.h5", source
    )
    loaded, loaded_dat = read_case_with_dat(case_path)
    assert loaded.iteration == 0
    assert loaded.time == pytest.approx(0.0)
    assert loaded_dat == dat_path


def test_read_case_with_dat_requires_paired_dat(tmp_path: Path) -> None:
    case_path = save_case(make_session(), tmp_path / "channel.cfdx.h5")

    with pytest.raises(FileNotFoundError):
        read_case_with_dat(case_path)


def test_save_case_preserves_existing_hdf5_mesh_data(tmp_path: Path) -> None:
    path = tmp_path / "channel.cfdx.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset("points", data=[0.0, 1.0, 2.0])

    save_case(make_session(), path)

    with h5py.File(path, "r") as h5:
        assert list(h5["points"][()]) == [0.0, 1.0, 2.0]
