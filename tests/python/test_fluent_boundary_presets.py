"""Physical BC presets should preserve the native C++ wire contract."""

import pytest

from cfdx.fluent_boundary_presets import (
    SIMPLE_BOUNDARY_CHOICES,
    UNSUPPORTED_CHOICE_REASONS,
    make_simple_boundary,
)
from cfdx.io.schema import BCType, BCValueType, BoundarySpec


def test_presets_cover_only_supported_physical_meanings():
    names = {name for name, _, _ in SIMPLE_BOUNDARY_CHOICES}
    assert {
        "velocity-inlet",
        "pressure-outlet",
        "no-slip-wall",
        "moving-wall",
        "symmetry",
    } <= names
    assert names.isdisjoint(UNSUPPORTED_CHOICE_REASONS)


def test_velocity_inlet_preserves_vector_and_native_schema():
    spec = make_simple_boundary(" in ", "velocity-inlet", velocity=(0, 1.25, 0))
    assert spec.type == BCType.INLET
    assert spec.value_type == BCValueType.FIXED
    assert spec.velocity_vector == [0, 1.25, 0]
    assert BoundarySpec.model_validate_json(spec.model_dump_json()) == spec


def test_pressure_outlet_preserves_gauge_pressure():
    spec = make_simple_boundary("out", "pressure-outlet", pressure_pa=-20.0)
    assert spec.type == BCType.PRESSURE_OUTLET
    assert spec.value_type == BCValueType.OUTLET_PRESSURE
    assert spec.pressure == -20.0
    assert spec.back_pressure == -20.0


def test_walls_symmetry_empty_outflow():
    wall = make_simple_boundary("wall", "no-slip-wall")
    assert wall.value_type == BCValueType.WALL_NO_SLIP
    moving = make_simple_boundary("wall", "moving-wall", velocity=(2, 0, 0))
    assert moving.value_type == BCValueType.FIXED
    assert moving.velocity_vector == [2, 0, 0]
    heated = make_simple_boundary("wall", "no-slip-wall", wall_temperature_k=320)
    assert heated.value_type == BCValueType.WALL_THERMAL
    assert heated.wall_temperature == 320
    assert make_simple_boundary("sym", "symmetry").type == BCType.SYMMETRY
    assert make_simple_boundary("front", "empty").type == BCType.EMPTY
    outflow = make_simple_boundary("out", "outflow")
    assert outflow.value_type == BCValueType.ZERO_GRADIENT


@pytest.mark.parametrize("kind", sorted(UNSUPPORTED_CHOICE_REASONS))
def test_unimplemented_fluent_types_are_rejected(kind):
    with pytest.raises(ValueError, match="non disponible"):
        make_simple_boundary("patch", kind)


@pytest.mark.parametrize(
    "data",
    [
        (),
        (1,),
        (1, 2),
        (1, 2, float("inf")),
        (1, 2, float("nan")),
        (True, 0, 0),
    ],
)
def test_inlet_rejects_malformed_velocity(data):
    with pytest.raises(ValueError, match="vitesse"):
        make_simple_boundary("in", "velocity-inlet", velocity=data)


def test_parameters_are_not_silently_assigned_to_wrong_bc():
    with pytest.raises(ValueError, match="pression"):
        make_simple_boundary("in", "velocity-inlet", velocity=(1, 0, 0), pressure_pa=0)
    with pytest.raises(ValueError, match="vitesse"):
        make_simple_boundary("in", "velocity-inlet")
    with pytest.raises(ValueError, match="kelvins"):
        make_simple_boundary("wall", "no-slip-wall", wall_temperature_k=-1)
    with pytest.raises(ValueError, match="temperature de paroi"):
        make_simple_boundary(
            "in", "velocity-inlet", velocity=(1, 0, 0), wall_temperature_k=300
        )
