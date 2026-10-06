from pathlib import Path

import pytest

from cfdx import CFDXSession, Project
from cfdx.case_io import save_case


def make_session() -> CFDXSession:
    session = CFDXSession()
    session.case.name = "channel"
    session.case.enable("incompressible")
    session.case.set_numerics(cfl=1.0)
    session.case.set_boundary("inlet", type="inlet", value="1.0")
    session.case_revision = 3
    session.mesh_revision = 2
    session.physics_revision = 4
    session.numerics_revision = 5
    return session


def test_project_normalizes_path_and_pairs_dat(tmp_path: Path) -> None:
    project = Project(tmp_path / "channel.cfdx.h5")

    assert project.name == "channel"
    assert project.path == tmp_path / "channel.cfdx.h5"
    assert project.dat_path == tmp_path / "channel.dat.h5"
    assert not project.exists()
    assert not project.has_dat()


def test_project_rejects_non_cfdx_case(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="canonical .cfdx.h5 extension"):
        Project(tmp_path / "channel.h5")


def test_project_save_load_is_setup_only(tmp_path: Path) -> None:
    project = Project(tmp_path / "channel.cfdx.h5")
    original = make_session()
    original.iteration = 17
    original.time = 2.5

    saved = project.save(original)
    loaded = project.load()

    assert saved == project.path
    assert project.exists()
    assert loaded.case.as_dict() == original.case.as_dict()
    assert loaded.case_revision == original.case_revision
    assert loaded.iteration == 0
    assert loaded.time == pytest.approx(0.0)
    assert not project.has_dat()


def test_project_validate_delegates_to_case_bundle(tmp_path: Path) -> None:
    project = Project(tmp_path / "channel.cfdx.h5")
    path = save_case(make_session(), project.path)

    # A minimal configuration artifact is not a self-contained mesh bundle,
    # so validation must report the existing case artifact as invalid rather
    # than silently treating persistence as a complete project.
    assert path == project.path
    with pytest.raises(ValueError, match="mesh data is missing"):
        project.validate()