from __future__ import annotations

import pytest

from cfdx.application import PostProcessingRegistry, PostProcessingSpec, build_post_processing_registry


def test_post_processing_registry_exposes_existing_capabilities() -> None:
    registry = build_post_processing_registry()

    assert registry.get("derived:u_magnitude").computation_owner.endswith(
        "derived_magnitude"
    )
    assert registry.get("mesh_quality:skewness").category == "mesh_quality"
    assert registry.get("probe:p").validation_status == "validated"


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
        "availability": "available",
        "computation_owner": "cfdx::core::compute_face_quality",
        "validation_status": "verified",
    }


def test_registry_filters_categories_and_availability() -> None:
    registry = build_post_processing_registry()

    assert all(spec.category == "probe" for spec in registry.by_category("probe"))
    assert registry.available() == registry.all()
