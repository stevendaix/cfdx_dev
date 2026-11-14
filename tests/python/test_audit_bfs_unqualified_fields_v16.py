from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "audit_bfs_unqualified_fields_v16.py"
SPEC = importlib.util.spec_from_file_location("bfs_field_audit_v16", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def _sample_vtu(path: Path, *, k: str = "1.0 2.0", omega: str = "4.0 5.0") -> None:
    pairs = {"k": k, "omega": omega, "nu_t": "0.1 0.2",
             "u_x": "-1.0 2.0", "u_y": "0.0 0.0", "u_z": "0.0 0.0",
             "p": "2.0 3.0", "wall_distance": "0.1 0.2"}
    arrays = "".join(f'<DataArray Name="{key}" format="ascii">{value}</DataArray>'
                     for key, value in pairs.items())
    path.write_text(f'<VTKFile type="UnstructuredGrid"><UnstructuredGrid>'
                    f'<Piece NumberOfCells="2"><CellData>{arrays}</CellData>'
                    f'</Piece></UnstructuredGrid></VTKFile>', encoding="utf-8")


def _run(path: Path) -> dict:
    return AUDIT.audit_vtu(path, u_ref=4.0, molecular_nu=0.01,
                           inlet_k=0.01, inlet_omega=1.0)


def test_diagnostics_do_not_certify_physics(tmp_path: Path):
    path = tmp_path / "a.vtu"
    _sample_vtu(path)
    report = _run(path)
    assert report["cells"] == 2
    assert report["physical_validation"] == "NOT_QUALIFIED"
    assert report["algebraic_convergence"] == "NOT_CERTIFIED_BY_VTU"
    assert report["fields"]["k"]["max_cell"] == 1
    assert report["diagnostic_ratios"]["max_k_to_inlet_k"] == 200
    assert report["diagnostic_ratios"]["negative_u_x_cells"] == 1


@pytest.mark.parametrize("k, omega", [("1.0 nan", "4.0 5.0"),
                                     ("-1.0 2.0", "4.0 5.0"),
                                     ("1.0 2.0", "4.0 0.0"),
                                     ("1.0 2.0 3.0", "4.0 5.0")])
def test_invalid_fields_fail_closed(tmp_path: Path, k: str, omega: str):
    path = tmp_path / "bad.vtu"
    _sample_vtu(path, k=k, omega=omega)
    with pytest.raises(ValueError):
        _run(path)


def test_invalid_reference_fails_closed(tmp_path: Path):
    path = tmp_path / "good.vtu"
    _sample_vtu(path)
    with pytest.raises(ValueError):
        AUDIT.audit_vtu(path, u_ref=0.0, molecular_nu=0.01,
                        inlet_k=0.01, inlet_omega=1.0)
