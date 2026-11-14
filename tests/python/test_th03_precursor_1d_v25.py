"""V25: solved 1-D SST *channel* precursor, not measured BFS validation."""

from dataclasses import replace
from unittest.mock import patch

import h5py
import numpy as np
import pytest

from scripts.generate_th03_openfoam_blockmesh import (
    VELOCITY,
    H,
    make_structured_mesh,
    write_case,
)
from scripts.initialize_th03_hybrid import (
    _geometry,
    divergence,
    initialise_case,
)
from scripts.solve_th03_precursor_1d import (
    PrecursorControls,
    fv_diffusion,
    grid,
    residual,
    solve_precursor,
)


def test_1d_fv_diffusion_matches_analytic_poiseuille():
    c = replace(PrecursorControls(), n_cells=160, stretch=0.0)
    edges, y = grid(c)
    nu = 1.56e-5
    # -nu * d2U/dy2 = 1, with U(0)=U(H)=0.
    sol, op = fv_diffusion(
        edges, y, np.full(len(y), nu), np.zeros(len(y)), np.ones(len(y)), 0.0, 0.0
    )
    exact = y * (c.height - y) / (2.0 * nu)
    assert np.max(np.abs(sol - exact)) / np.max(exact) < 1.0e-4
    assert np.max(np.abs(residual(sol, op))) < 2.0e-11


def test_fv_diffusion_rejects_negative_diagonal_sink():
    c = PrecursorControls(n_cells=16)
    edges, y = grid(c)
    with pytest.raises(ValueError, match="non-physical"):
        fv_diffusion(
            edges, y, np.ones(len(y)), np.full(len(y), -1.0), np.ones(len(y)), 0, 0
        )


def test_sst_channel_precursor_is_discretely_converged_with_debit():
    case = solve_precursor(PrecursorControls())
    d = case["diagnostics"]
    assert d["converged"]
    for name in (
        "U_equation_defect",
        "k_equation_defect",
        "omega_equation_defect",
        "mean_velocity_error",
        "update_relative",
    ):
        assert d[name] < 2.0e-7
    assert d["maximum_nut_over_nu"] > 1.0
    assert d["minimum_k"] > 0 and d["minimum_omega"] > 0
    assert abs(d["bulk_velocity_m_s"] - VELOCITY) < 1.0e-10
    assert np.isfinite(case["U"]).all()


def test_underresolved_first_wall_cell_is_rejected_even_if_residuals_pass():
    r = solve_precursor(PrecursorControls(n_cells=64, stretch=2.3))
    assert r["diagnostics"]["discrete_equations_converged"]
    assert not r["diagnostics"]["wall_resolved"]
    assert not r["diagnostics"]["converged"]


def test_precursor_grid_sensitivity_requires_stretched_wall_mesh():
    r96 = solve_precursor(PrecursorControls(n_cells=96))
    r192 = solve_precursor(PrecursorControls(n_cells=192))
    assert r96["diagnostics"]["converged"]
    assert r192["diagnostics"]["converged"]
    g96 = r96["diagnostics"]["momentum_pressure_source_m_s2"]
    g192 = r192["diagnostics"]["momentum_pressure_source_m_s2"]
    assert abs(g96 - g192) / g192 < 0.005


@pytest.fixture(scope="module")
def precursor_case(tmp_path_factory):
    base = tmp_path_factory.mktemp("v25-precursor")
    source = base / "source.cfdx.h5"
    output = base / "solved-1d.cfdx.h5"
    write_case(
        source,
        make_structured_mesh("smoke"),
        mesh_source="v25-discrete-precursor",
        max_iterations=2,
    )
    report = initialise_case(source, output, sst_style="precursor_1d", require_native_geometry=False)
    return source, output, report


def test_1d_precursor_integration_preserves_flux_and_disclaims_qualification(
    precursor_case,
):
    source, output, report = precursor_case
    assert report["precursor_1d_evidence"]["converged"]
    assert report["normalized_after_l1"] < 1e-9
    assert report["boundary_immutable"]
    with h5py.File(output) as f:
        assert f.attrs["hybrid_qualification"] == "NOT_QUALIFIED"
        assert (
            f.attrs["hybrid_inlet_profile"]
            == "CONVERGED_1D_CHANNEL_NOT_MEASURED_BFS_INFLOW"
        )
        geo = _geometry(f)
        n = int(f.attrs["n_cells"])
        phi = f["fields/face/mdot"][:]
        assert np.sum(np.abs(divergence(geo["owner"], geo["nb"], phi, n))) < 1e-10
        u = f["fields/vector/U"][:]
        k = f["fields/scalar/k"][:]
        omega = f["fields/scalar/omega"][:]
        assert np.isfinite(u).all()
        assert np.all(k > 0) and np.all(omega > 0)
        near_station = (geo["centres"][:, 0] / H > -10) & (geo["centres"][:, 0] < 0)
        assert np.std(u[near_station, 0]) > 1.0
    with h5py.File(source) as original:
        assert "fields" not in original


def test_refuse_uncertified_1d_initialization_without_output(precursor_case, tmp_path):
    source, _, _ = precursor_case
    target = tmp_path / "refused.cfdx.h5"
    with (
        patch(
            "scripts.solve_th03_precursor_1d.solve_precursor",
            return_value={"diagnostics": {"converged": False}},
        ),
        pytest.raises(ArithmeticError, match="unconverged 1D"),
    ):
        initialise_case(source, target, sst_style="precursor_1d", require_native_geometry=False)
    assert not target.exists()


def test_velocity_only_mode_never_claims_verified_sst_precursor(
    precursor_case, tmp_path
):
    source, _, _ = precursor_case
    legacy = tmp_path / "legacy.cfdx.h5"
    velocity_only = tmp_path / "velocity_only.cfdx.h5"
    initialise_case(source, legacy, sst_style="wall_aware", require_native_geometry=False)
    report = initialise_case(source, velocity_only, sst_style="precursor_1d_u_only", require_native_geometry=False)
    assert report["profile_validation"] == (
        "1D_U_DISCRETE_CERTIFIED_K_OMEGA_UNVERIFIED_NOT_2D_QUALIFIED"
    )
    with h5py.File(legacy) as ref, h5py.File(velocity_only) as new:
        assert np.array_equal(ref["fields/scalar/k"][:], new["fields/scalar/k"][:])
        assert np.array_equal(
            ref["fields/scalar/omega"][:], new["fields/scalar/omega"][:]
        )
        assert not np.array_equal(ref["fields/vector/U"][:], new["fields/vector/U"][:])
        assert new.attrs["hybrid_initialisation"] == (
            "self_grid_channel_1d_U_only_SST_surrogate_v2"
        )
        assert new.attrs["hybrid_qualification"] == "NOT_QUALIFIED"
