"""Tests for the Python-to-CFDX canonical numerical selection boundary."""

import json

import h5py

from cfdx.io.converter import convert
from cfdx.io.numerical_selection import build_numerical_selection
from cfdx.io.schema import NumericalScheme


def test_build_numerical_selection_uses_registry_keys():
    scheme = NumericalScheme(
        gradient_operator="least_squares",
        momentum_scheme="second_order_upwind",
        momentum_interpolation="linear",
        coupled_solver="SIMPLEC",
        transient_scheme="bdf2",
        preconditioner="native_amg",
    )
    selection = build_numerical_selection(scheme)

    assert [(e.family, e.configuration_key) for e in selection.entries] == [
        ("gradient", "numerics.gradient.least_squares"),
        ("convection", "numerics.convection.second_order_upwind"),
        ("interpolation", "numerics.interpolation.linear"),
        ("temporal", "numerics.temporal.bdf2"),
        ("pressure_velocity", "pressure_velocity.simplec"),
        ("preconditioner", "preconditioner.native_amg"),
    ]
    assert selection.required_families == [
        "gradient",
        "convection",
        "interpolation",
        "temporal",
        "pressure_velocity",
        "preconditioner",
    ]


def test_unmapped_source_value_is_not_silently_rewritten():
    scheme = NumericalScheme(momentum_scheme="source_specific_unknown")
    selection = build_numerical_selection(scheme)
    assert all(e.family != "convection" for e in selection.entries)


def test_converted_case_contains_canonical_selection(su2_case, tmp_path):
    output = tmp_path / "case.cfdx.h5"
    result = convert(su2_case, output=output)

    selection = result.case.numerics.selection
    assert selection.entries
    assert all(e.configuration_key for e in selection.entries)

    with h5py.File(output, "r") as h5:
        payload = json.loads(h5.attrs["case_setup_json"])
        assert payload["numerics"]["selection"]["entries"]
        assert "numerical_selection_report" in h5.attrs
        report = h5.attrs["numerical_selection_report"]
        assert "scheme[" in report
