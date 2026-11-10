"""SST hotspot auditing must never pretend to certify numerical convergence."""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[2] / "scripts" / "audit_bfs_sst_hotspots_v17.py"
SPEC = importlib.util.spec_from_file_location("audit_bfs_sst_hotspots_v17", MODULE)
assert SPEC is not None and SPEC.loader is not None
m = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(m)


def _vtu(path: Path, *, k: str = "1 2", p: str = "10 12") -> None:
    names = {"k": k, "omega": "4 5", "nu_t": "0.01 0.02",
             "sst_strain_rate": "20 30", "p": p, "u_x": "2 2",
             "u_y": "0 0", "u_z": "0 0"}
    data = "\n".join(f'<DataArray Name="{name}" format="ascii">{values}</DataArray>'
                     for name, values in names.items())
    path.write_text(f'''<VTKFile type="UnstructuredGrid"><UnstructuredGrid>
<Piece NumberOfCells="2" NumberOfPoints="3">
<Points><DataArray format="ascii">0 0 0 1 0 0 2 0 0</DataArray></Points>
<CellData>{data}</CellData>
<Cells><DataArray Name="connectivity" format="ascii">0 1 1 2</DataArray>
<DataArray Name="offsets" format="ascii">2 4</DataArray></Cells>
</Piece></UnstructuredGrid></VTKFile>''')


def test_sources_and_static_pressure(tmp_path: Path) -> None:
    path = tmp_path / "a.vtu"
    _vtu(path)
    result = m.audit(path, rho=1.0, nu=1e-5, u_ref=2, k_inlet=0.1, top_n=1)
    assert result["status"] == "DIAGNOSTIC_ONLY_NOT_QUALIFIED"
    kmax = next(s for s in result["hotspots"] if s["cell"] == 1)
    assert kmax["k_m2_s2"] == 2
    assert kmax["p_static_Pa"] == pytest.approx(12 - 4/3)
    assert kmax["Cp_star_minus_Cp_static"] == pytest.approx(2/3)
    assert kmax["net_k_source_W_m3"] == pytest.approx(
        min(0.02 * 30**2, 10 * .09 * 2 * 5) - .09 * 2 * 5
    )


def test_rejects_nan_and_invalid_mesh(tmp_path: Path) -> None:
    path = tmp_path / "bad.vtu"
    _vtu(path, k="NaN 2")
    with pytest.raises(ValueError):
        m.audit(path, 1, 1e-5, 2, 0.1)
    _vtu(path)
    txt = path.read_text().replace('2 4</DataArray>', '2 6</DataArray>')
    path.write_text(txt)
    with pytest.raises(ValueError):
        m.audit(path, 1, 1e-5, 2, 0.1)


def test_requires_positive_refs(tmp_path: Path) -> None:
    path = tmp_path / "a.vtu"
    _vtu(path)
    with pytest.raises(ValueError):
        m.audit(path, 0, 1e-5, 2, 0.1)
