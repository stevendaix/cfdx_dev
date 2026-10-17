from cfdx.physics_setup import (
    TURBULENCE_MODELS,
    TURBULENCE_CATALOG,
    TurbulenceSelectionContext,
    TurbulenceStatus,
    default_parameters,
    physics_spec,
    recommend_turbulence_models,
    turbulence_model,
    turbulence_model_requirements,
    TURBULENCE_MODEL_GAPS,
)
from cfdx.setup_model import ParameterType


def test_turbulence_catalog_exposes_all_current_model_states():
    keys = {model.key for model in TURBULENCE_CATALOG}
    assert {
        "LAMINAR", "KEPSILON", "RNG_KEPSILON", "REALIZABLE_KEPSILON",
        "KOMEGA", "SST", "SPALART_ALLMARAS", "SMAGORINSKY", "WALE",
        "DYNAMIC_KEQN", "DES", "DDES", "IDDES",
    } <= keys
    assert turbulence_model("SST").status is TurbulenceStatus.SOLVER_READY
    assert turbulence_model("DES").status is TurbulenceStatus.KERNEL_ONLY
    assert turbulence_model("DYNAMIC_KEQN").status is TurbulenceStatus.PLANNED


def test_turbulence_dropdown_contains_only_solver_ready_models():
    params = {p.name: p for p in default_parameters(physics_spec("turbulence"))}
    assert params["model"].kind is ParameterType.CHOICE
    assert params["model"].choices == TURBULENCE_MODELS
    assert "SST" in params["model"].choices
    assert "SMAGORINSKY" not in params["model"].choices
    assert "DES" not in params["model"].choices


def test_turbulence_recommendations_are_explicit_and_non_mutating():
    context = TurbulenceSelectionContext(
        application="external",
        dimensions=3,
        steady=True,
        wall_bounded=True,
        strong_separation=True,
        transition_expected=True,
    )
    recommendations = recommend_turbulence_models(context)
    assert recommendations[0].model == "SST"
    assert any("transition" in warning.lower() for warning in recommendations[0].warnings)
    assert any("separation" in reason.lower() for reason in recommendations[0].reasons)


def test_turbulence_selection_rejects_invalid_context():
    try:
        recommend_turbulence_models(TurbulenceSelectionContext(dimensions=4))
    except ValueError:
        pass
    else:
        raise AssertionError("invalid dimensionality must be rejected")


def test_legacy_turbulence_settings_normalize_to_structured_schema():
    from cfdx.physics_setup import turbulence_model_from_case, validate_turbulence_selection

    normalized = turbulence_model_from_case({"model": "SST"})
    assert normalized["family"] == "RANS"
    assert normalized["model"] == "SST"
    assert normalized["wall_treatment"] == "resolved"
    assert normalized["transition"] == {"model": "none"}
    assert validate_turbulence_selection({"model": "SST"})["model"] == "SST"


def test_non_solver_ready_turbulence_models_are_rejected_for_production_cases():
    from cfdx.physics_setup import validate_turbulence_selection

    # Guard precedence matters here. Every KERNEL_ONLY closure is LES or
    # HYBRID, so the default steady solver type trips the transient guard
    # first and the user is told the actionable reason. Production readiness
    # is the only remaining objection once a transient driver is selected.
    try:
        validate_turbulence_selection({"model": "DES"})
    except ValueError as exc:
        assert "transient 3-D" in str(exc)
        assert "kernel_only" not in str(exc)
    else:
        raise AssertionError("kernel-only turbulence model must not be production-selectable")

    try:
        validate_turbulence_selection({"model": "DES", "solver_type": "transient"})
    except ValueError as exc:
        assert "kernel_only" in str(exc)
    else:
        raise AssertionError("kernel-only turbulence model must not be production-selectable")


def test_steady_scale_resolving_models_are_rejected_explicitly():
    from cfdx.physics_setup import validate_turbulence_selection

    for model in ("SMAGORINSKY", "WALE", "DES", "DDES", "IDDES"):
        try:
            validate_turbulence_selection({"model": model, "solver_type": "steady"})
        except ValueError as exc:
            assert "transient 3-D" in str(exc)
        else:
            raise AssertionError(f"{model} must reject steady selection")


def test_transient_hybrid_guard_is_independent_of_gui_filtering():
    from cfdx.physics_setup import validate_turbulence_selection

    try:
        validate_turbulence_selection({"model": "DDES", "solver_type": "steady"})
    except ValueError:
        pass
    else:
        raise AssertionError("manual steady DDES case must be rejected")

    # The same guard must not reject the intended transient configuration
    # merely because the model remains KERNEL_ONLY.
    try:
        validate_turbulence_selection({"model": "DDES", "solver_type": "transient"})
    except ValueError as exc:
        assert "kernel_only" in str(exc)
    else:
        raise AssertionError("DDES must still be blocked by production readiness")



def test_every_turbulence_model_has_explicit_requirements_and_gaps():
    for model in TURBULENCE_CATALOG:
        required, missing = turbulence_model_requirements(model.key)
        assert model.key in TURBULENCE_MODEL_GAPS
        assert isinstance(required, tuple)
        assert isinstance(missing, tuple)
        assert all(item.strip() for item in required)
        assert all(item.strip() for item in missing)
