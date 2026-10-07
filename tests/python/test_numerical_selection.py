"""Tests for the Python-to-CFDX canonical numerical selection boundary."""

import json
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src" / "cfdx" / "python"))

import h5py

from cfdx.io.converter import convert
from cfdx.io.numerical_selection import build_numerical_selection, validate_numerical_selection
from cfdx.io.schema import NumericalScheme


def test_build_numerical_selection_uses_registry_keys():
    scheme = NumericalScheme(
        gradient_operator="least_squares",
        momentum_scheme="second_order_upwind",
        momentum_interpolation="linear",
        coupled_solver="SIMPLEC",
        transient_scheme="bdf2",
        preconditioner="native_amg",
        linear_solver="fgmres",
    )
    selection = build_numerical_selection(scheme)

    assert [(e.family, e.configuration_key) for e in selection.entries] == [
        ("gradient", "numerics.gradient.least_squares"),
        ("convection", "numerics.convection.second_order_upwind"),
        ("temporal", "numerics.temporal.bdf2"),
        ("pressure_velocity", "pressure_velocity.simplec"),
        ("linear_solver", "linear.fgmres"),
        ("preconditioner", "preconditioner.native_amg"),
    ]
    assert selection.required_families == [
        "gradient",
        "convection",
        "temporal",
        "pressure_velocity",
        "linear_solver",
        "preconditioner",
    ]


def test_coupled_schur_is_a_distinct_canonical_family():
    scheme = NumericalScheme(coupled_solver="coupled", coupled_schur="pcd")
    selection = build_numerical_selection(scheme)
    assert ("pressure_velocity", "pressure_velocity.coupled") in [
        (e.family, e.configuration_key) for e in selection.entries
    ]
    assert ("schur", "schur.pcd") in [
        (e.family, e.configuration_key) for e in selection.entries
    ]
    assert "schur" in selection.required_families
    assert validate_numerical_selection(scheme) == []


def test_unmapped_source_value_is_not_silently_rewritten():
    scheme = NumericalScheme(momentum_scheme="source_specific_unknown")
    selection = build_numerical_selection(scheme)
    assert all(e.family != "convection" for e in selection.entries)
    errors = validate_numerical_selection(scheme)
    assert errors == ["unmapped numerical setting for convection: 'source_specific_unknown'"]


def test_known_source_values_have_no_mapping_gaps():
    scheme = NumericalScheme(
        gradient_operator="green_gauss_cell",
        momentum_scheme="upwind",
        coupled_solver="SIMPLE",
    )
    assert validate_numerical_selection(scheme) == []


def test_adapter_recorded_unmapped_setting_is_reported():
    """An adapter that cannot map a source scheme must not drop it silently.

    The SU2 adapter leaves the canonical convection field empty for families
    with no faithful CFDX counterpart and records the source value instead.
    An empty canonical field is indistinguishable from an unspecified one, so
    the recorded entry is what keeps the conversion blocking.
    """
    scheme = NumericalScheme(
        momentum_scheme="",
        unmapped_settings=["convection: 'JST' (no faithful CFDX convection registry equivalent)"],
    )
    errors = validate_numerical_selection(scheme)
    assert errors == [
        "unmapped numerical setting for convection: 'JST' "
        "(no faithful CFDX convection registry equivalent)"
    ]


def test_converted_case_with_unmappable_scheme_is_blocking(tmp_path):
    data_dir = Path(__file__).parent / ".." / "data" / "su2"
    source = tmp_path / "case.su2"
    shutil.copy(data_dir / "mesh_NACA0012_inv.su2", source)
    shutil.copy(data_dir / "inv_NACA0012_basic.cfg", tmp_path / "case.cfg")

    result = convert(source)

    assert result.success is False
    blocking = [f for f in result.gap_analysis.findings if f.severity.value == "unsupported_blocking"]
    assert any(f.feature == "canonical_selection" and "JST" in f.detail for f in blocking)


def test_converted_case_contains_canonical_selection(tmp_path):
    data_dir = Path(__file__).parent / ".." / "data" / "su2"
    source = tmp_path / "case.su2"
    shutil.copy(data_dir / "mesh_NACA0012_inv.su2", source)
    shutil.copy(data_dir / "inv_NACA0012_basic.cfg", tmp_path / "case.cfg")
    output = tmp_path / "case.cfdx.h5"
    result = convert(source, output=output)

    selection = result.case.numerics.selection
    assert selection.entries
    assert all(e.configuration_key for e in selection.entries)

    with h5py.File(output, "r") as h5:
        payload = json.loads(h5.attrs["case_setup_json"])
        assert payload["numerics"]["selection"]["entries"]
        assert "numerical_selection_report" in h5.attrs
        report = h5.attrs["numerical_selection_report"]
        assert "scheme[" in report
