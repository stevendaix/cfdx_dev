import pytest

from cfdx.probe import (
    nearest_probe,
    Probe,
    ProbeCatalog,
    ProbeField,
)


def test_nearest_probe() -> None:
    assert nearest_probe(((0, 0), (2, 0), (0, 2)), (1, 4, 8), (1.8, 0.1)) == 4.0


def test_probe_validates_alignment() -> None:
    with pytest.raises(ValueError):
        nearest_probe(((0,),), (), (0,))


def test_probe_field_cli_matches_solver_option() -> None:
    assert ProbeField.PRESSURE.cli == "p"
    assert ProbeField.U_MAGNITUDE.cli == "u_mag"
    assert ProbeField.from_cli("p") is ProbeField.PRESSURE
    assert ProbeField.from_cli("u_x") is ProbeField.U_X
    assert ProbeField.from_cli("u_mag") is ProbeField.U_MAGNITUDE
    with pytest.raises(ValueError):
        ProbeField.from_cli("bogus")


def test_probe_spec_matches_solver_cli_format() -> None:
    probe = Probe("pressure", 0.5, 0.0, 1.0, ProbeField.PRESSURE)
    assert probe.spec() == "pressure:0.5,0.0,1.0:p"


def test_probe_rejects_empty_name_and_non_finite_coords() -> None:
    with pytest.raises(ValueError):
        Probe("", 0, 0, 0)
    with pytest.raises(ValueError):
        Probe("ok", float("nan"), 0, 0)


def test_probe_rejects_name_with_comma() -> None:
    with pytest.raises(ValueError):
        Probe("bad,name", 0, 0, 0)


def test_probe_catalog_specs_and_csv_path(tmp_path) -> None:
    catalog = ProbeCatalog([
        Probe("p_wall", 1.0, 2.0, 3.0, ProbeField.PRESSURE),
        Probe("u_inlet", 0.0, 0.0, 0.0, ProbeField.U_X),
    ])
    assert catalog.specs() == ["p_wall:1.0,2.0,3.0:p", "u_inlet:0.0,0.0,0.0:u_x"]
    case_file = tmp_path / "case.cfdx.h5"
    assert catalog.csv_path(case_file) == tmp_path / "case.cfdx.probes.csv"
    assert catalog.csv_path(case_file).parent.exists()


def test_probe_catalog_rejects_duplicate_names() -> None:
    catalog = ProbeCatalog()
    catalog.add(Probe("p", 0, 0, 0, ProbeField.PRESSURE))
    with pytest.raises(ValueError):
        catalog.add(Probe("p", 1, 1, 1, ProbeField.U_X))


def test_probe_catalog_remove() -> None:
    catalog = ProbeCatalog([Probe("a", 0, 0, 0), Probe("b", 1, 1, 1)])
    assert catalog.remove("a") is True
    assert catalog.names() == ["b"]
    assert catalog.remove("missing") is False
