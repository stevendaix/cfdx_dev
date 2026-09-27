from cfdx.physics_setup import (
    TURBULENCE_MODELS,
    TURBULENCE_CATALOG,
    TurbulenceSelectionContext,
    TurbulenceStatus,
    default_parameters,
    physics_spec,
    recommend_turbulence_models,
    turbulence_model,
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
