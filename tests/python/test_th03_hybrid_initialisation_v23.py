"""V23 hybrid: no cross-grid transfer; conservative initial flux, positive SST."""


import h5py
import numpy as np
import pytest

from scripts.generate_th03_openfoam_blockmesh import make_structured_mesh, write_case
from scripts.initialize_th03_hybrid import _geometry, divergence, initialise_case


@pytest.fixture(scope="module")
def case(tmp_path_factory):
    root = tmp_path_factory.mktemp("hybrid")
    raw = root / "smoke_raw.cfdx.h5"
    write_case(
        raw,
        make_structured_mesh("smoke"),
        mesh_source="hybrid-v23-test",
        max_iterations=2,
    )
    hybrid = root / "smoke_hybrid.cfdx.h5"
    report = initialise_case(raw, hybrid)
    return raw, hybrid, report


def test_projection_conserves_each_cell_and_inlet(case):
    _raw, hybrid, report = case
    with h5py.File(hybrid) as f:
        geo = _geometry(f)
        flux = f["fields/face/mdot"][:]
        div = divergence(geo["owner"], geo["nb"], flux, len(geo["centres"]))
    assert np.abs(div).sum() < 1.0e-10
    assert report["normalized_after_l1"] < 1.0e-9
    assert report["after_l1_kg_s"] < report["before_l1_kg_s"] * 1.0e-8
    assert abs(report["inlet_kg_s"] + report["outlet_kg_s"]) < 1e-10
    assert report["boundary_immutable"]


def test_fields_are_nonuniform_positive_and_same_grid(case):
    raw, hybrid, report = case
    with h5py.File(hybrid) as f:
        u = f["fields/vector/U"][:]
        k = f["fields/scalar/k"][:]
        omega = f["fields/scalar/omega"][:]
        assert len(u) == int(f.attrs["n_cells"]) == 712
        assert u.shape == (712, 3)
        assert np.std(u[:, 0]) > 1.0
        assert np.std(k) > 1.0e-6
        assert np.std(omega) > 100.0
        assert (k > 0).all() and (omega > 0).all()
        assert (
            np.isfinite(u).all() and np.isfinite(k).all() and np.isfinite(omega).all()
        )
        assert f.attrs["hybrid_qualification"] == "NOT_QUALIFIED"
        assert f.attrs["topology_hash"]
    with h5py.File(raw) as f:
        assert "fields" not in f
    assert report["fields_imported"] == ["U", "k", "omega", "mdot"]


def test_reject_accidental_overwrite(case):
    raw, hybrid, _ = case
    with pytest.raises(ValueError, match="new, absent"):
        initialise_case(raw, hybrid)
    with pytest.raises(ValueError, match="new, absent"):
        initialise_case(raw, raw)


def test_reject_bad_bl_length(case, tmp_path):
    raw, _, _ = case
    with pytest.raises(ValueError, match="invalid BL"):
        initialise_case(raw, tmp_path / "bad.cfdx.h5", delta_h=-1)
    assert not (tmp_path / "bad.cfdx.h5").exists()
