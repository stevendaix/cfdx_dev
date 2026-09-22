from pathlib import Path

import h5py
import pytest

from cfdx import CFDXSession
from cfdx.case_io import read_case, read_case_with_dat, save_case, save_case_with_dat


def make_session() -> CFDXSession:
    session = CFDXSession()
    session.case.name = "channel"
    session.case.enable("energy")
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


def test_case_hdf5_roundtrip_preserves_configuration_and_checkpoint(tmp_path: Path) -> None:
    original = make_session()
    path = save_case(original, tmp_path / "channel.cfdx.h5")

    loaded = read_case(path)

    assert loaded.case.as_dict() == original.case.as_dict()
    assert loaded.case_revision == 7
    assert loaded.mesh_revision == 3
    assert loaded.physics_revision == 5
    assert loaded.numerics_revision == 6
    assert loaded.iteration == 120
    assert loaded.time == pytest.approx(2.5)


def test_case_hdf5_contains_no_separate_checkpoint_file(tmp_path: Path) -> None:
    path = save_case(make_session(), tmp_path / "channel.cfdx.h5")

    assert path.is_file()
    assert not path.with_suffix(".json").exists()
    with h5py.File(path, "r") as h5:
        assert "case/config" in h5
        assert "runtime/checkpoint" in h5


def test_save_and_read_case_with_dat(tmp_path: Path) -> None:
    source = tmp_path / "solver.dat"
    source.write_text("CFDX DAT restart\niteration=120\ntime=2.5\n", encoding="utf-8")
    case_path, dat_path = save_case_with_dat(
        make_session(), tmp_path / "channel.cfdx.h5", source
    )

    loaded, loaded_dat = read_case_with_dat(case_path)

    assert loaded.case.as_dict() == make_session().case.as_dict()
    assert loaded.iteration == 120
    assert loaded.time == pytest.approx(2.5)
    assert loaded_dat == dat_path
    assert loaded_dat.read_text(encoding="utf-8") == source.read_text(encoding="utf-8")


def test_read_case_with_dat_rejects_modified_dat(tmp_path: Path) -> None:
    source = tmp_path / "solver.dat"
    source.write_text("restart A\n", encoding="utf-8")
    case_path, dat_path = save_case_with_dat(
        make_session(), tmp_path / "channel.cfdx.h5", source
    )
    dat_path.write_text("restart B\n", encoding="utf-8")

    with pytest.raises(ValueError, match="incompatible"):
        read_case_with_dat(case_path)


def test_read_case_with_dat_requires_paired_restart(tmp_path: Path) -> None:
    case_path = save_case(make_session(), tmp_path / "channel.cfdx.h5")

    with pytest.raises(ValueError, match="paired DAT"):
        read_case_with_dat(case_path)


def test_save_case_preserves_existing_hdf5_mesh_data(tmp_path: Path) -> None:
    path = tmp_path / "channel.cfdx.h5"
    with h5py.File(path, "w") as h5:
        h5.create_dataset("points", data=[0.0, 1.0, 2.0])

    save_case(make_session(), path)

    with h5py.File(path, "r") as h5:
        assert list(h5["points"][()]) == [0.0, 1.0, 2.0]
