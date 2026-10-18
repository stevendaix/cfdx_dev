from __future__ import annotations

import pytest

from cfdx.application import (
    PostProcessingRegistry,
    PostProcessingSpec,
    build_post_processing_registry,
)


def test_post_processing_registry_exposes_existing_capabilities() -> None:
    registry = build_post_processing_registry()

    assert registry.get("derived:u_magnitude").computation_owner.endswith(
        "execute_derived"
    )
    assert registry.get("derived:u_magnitude").inputs == ("U",)
    assert registry.get("mesh_quality:skewness").category == "mesh_quality"
    assert registry.get("probe:p").validation_status == "validated"


def test_mesh_quality_computation_is_owned_by_the_cpp_kernel_not_python() -> None:
    """Mesh-quality values must not be duplicated in the Python application layer.

    The registry declares mesh-quality as owned by ``cfdx::core::compute_face_quality``
    (the C++ geometry kernel). The Python layer must not re-derive skewness or
    non-orthogonality: that would silently diverge from the discretisation the
    solver actually used. This test pins the ownership string so a future
    refactor cannot swap in a Python duplicate without failing the contract.

    The values are "conditional" rather than "available" because, although the
    C++ kernel computes them, there is not yet an export path that exposes the
    per-face values through the results contract. The registry is capability
    declaration, not executable evidence, and must not imply a path that does
    not exist.
    """
    registry = build_post_processing_registry()
    for spec in registry.by_category("mesh_quality"):
        assert spec.computation_owner == "cfdx::core::compute_face_quality"
        assert spec.availability == "conditional"


def test_post_processing_registry_marks_unavailable_fields() -> None:
    """The registry must not advertise fields the production solver does not emit.

    The production solver's per-iteration VTU frames carry pressure only
    (``apps/production_solver/main.cpp`` writes ``{{"p", pressure}}``). Velocity
    and the velocity magnitude are therefore "conditional": they are real,
    evidenced capabilities, but only from a richer artifact (a DAT checkpoint)
    than the readable result series. The registry must surface that distinction
    instead of silently claiming every quantity is available from a VTU frame.
    """
    registry = build_post_processing_registry()

    assert registry.get("field:U").availability == "conditional"
    assert registry.get("derived:u_magnitude").availability == "conditional"
    assert registry.get("field:p").availability == "available"
    assert all(spec.availability != "unavailable" for spec in registry.all())


def test_post_processing_registry_has_no_duplicate_ids() -> None:
    registry = build_post_processing_registry()
    specs = registry.all()

    assert len(specs) == len({spec.stable_id for spec in specs})
    assert {spec.category for spec in specs} >= {
        "native_field",
        "derived_field",
        "mesh_quality",
        "probe",
    }


def test_registry_rejects_duplicate_and_unknown_capabilities() -> None:
    registry = PostProcessingRegistry()
    spec = PostProcessingSpec(
        stable_id="derived:test",
        name="Test",
        category="derived_field",
        units=None,
        inputs=(),
        availability="available",
        computation_owner="test.owner",
        validation_status="implemented",
    )
    registry.register(spec)

    with pytest.raises(ValueError, match="already registered"):
        registry.register(spec)

    with pytest.raises(KeyError, match="unknown post-processing capability"):
        registry.get("derived:missing")


def test_post_processing_spec_round_trips_as_metadata() -> None:
    spec = build_post_processing_registry().get("mesh_quality:non_orthogonality")

    assert spec.to_dict() == {
        "stable_id": "mesh_quality:non_orthogonality",
        "name": "Face non-orthogonality",
        "category": "mesh_quality",
        "units": "rad",
        "inputs": ["face_centre", "owner_centre", "neighbour_centre", "Sf"],
        # Conditional: the C++ kernel computes the value, but there is not yet
        # an export path exposing per-face values through the results contract.
        "availability": "conditional",
        "computation_owner": "cfdx::core::compute_face_quality",
        "validation_status": "verified",
    }


def test_registry_filters_categories_and_availability() -> None:
    registry = build_post_processing_registry()

    assert all(spec.category == "probe" for spec in registry.by_category("probe"))
    native = {spec.stable_id for spec in registry.by_category("native_field")}
    assert native == {"field:p", "field:U"}


def test_execute_derived_is_gated_by_registry() -> None:
    """A derived field is only computable when the registry declares it executable.

    The production solver stores ``U`` as one cell field with three components
    per cell, so the dataset exposes ``U -> [[u_x, u_y, u_z], ...]`` and the
    magnitude is computed per cell from the vector array.
    """
    from cfdx.application import execute_derived

    class _Dataset:
        def get(self, name: str):
            if name == "U":
                return [[1.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
            return None

        def names(self) -> tuple[str, ...]:
            return ("U",)

    magnitude = execute_derived("derived:u_magnitude", _Dataset())
    assert magnitude == (1.0, 2.0)

    with pytest.raises(KeyError):
        execute_derived("derived:not_registered", _Dataset())

    class _EmptyDataset:
        def get(self, name: str):
            return None

        def names(self) -> tuple[str, ...]:
            return ()

    with pytest.raises(KeyError):
        execute_derived("derived:u_magnitude", _EmptyDataset())
