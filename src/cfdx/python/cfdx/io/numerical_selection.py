"""Canonical N1 numerical selections for the Python adapter boundary.

The adapter layer owns translation from source-specific numerical names to the
stable CFDX registry configuration keys.  The solver must consume the canonical
selection, not these source-specific strings.
"""

from __future__ import annotations

from cfdx.io.schema import NumericalScheme, NumericalSelection, NumericalSelectionConfig


_GRADIENT_KEYS = {
    "green_gauss_cell": "numerics.gradient.gauss",
    "green_gauss": "numerics.gradient.gauss",
    "green_gauss_cell_based": "numerics.gradient.gauss",
    "green_gauss_vertex": "numerics.gradient.gauss_vertex",
    "green_gauss_node": "numerics.gradient.gauss_vertex",
    "least_squares": "numerics.gradient.least_squares",
    "weighted_least_squares": "numerics.gradient.weighted_least_squares",
}

_CONVECTION_KEYS = {
    "first_order": "numerics.convection.upwind",
    "upwind": "numerics.convection.upwind",
    "second_order_upwind": "numerics.convection.second_order_upwind",
    "central": "numerics.convection.central",
    "tvd_minmod": "numerics.convection.tvd.minmod",
    "tvd_vanleer": "numerics.convection.tvd.vanleer",
    "tvd_superbee": "numerics.convection.tvd.superbee",
    "tvd_vanalbada": "numerics.convection.tvd.vanalbada",
    "tvd_mc": "numerics.convection.tvd.mc",
}

_INTERPOLATION_KEYS = {
    "linear": "numerics.interpolation.linear",
    "upwind": "numerics.interpolation.upwind",
}

_TEMPORAL_KEYS = {
    "euler_explicit": "numerics.temporal.euler_explicit",
    "explicit_euler": "numerics.temporal.euler_explicit",
    "euler_implicit": "numerics.temporal.euler_implicit",
    "implicit_euler": "numerics.temporal.euler_implicit",
    "crank_nicolson": "numerics.temporal.crank_nicolson",
    "bdf2": "numerics.temporal.bdf2",
    "rk2": "numerics.temporal.rk2",
    "rk3": "numerics.temporal.rk3",
}

_PRESSURE_VELOCITY_KEYS = {
    "simple": "pressure_velocity.simple",
    "simplec": "pressure_velocity.simplec",
    "piso": "pressure_velocity.piso",
    "pimple": "pressure_velocity.pimple",
    "fractional_step": "pressure_velocity.fractional_step",
    "fractional step": "pressure_velocity.fractional_step",
    "coupled": "pressure_velocity.coupled",
}

_PRECONDITIONER_KEYS = {
    "native_amg": "preconditioner.native_amg",
    "amg": "preconditioner.native_amg",
    "smoothed_aggregation_amg": "preconditioner.smoothed_aggregation_amg",
    "native_fieldsplit": "preconditioner.native_fieldsplit",
    "coupled_block_schur": "preconditioner.coupled_block_schur",
}


def _normalise(value: str) -> str:
    return " ".join(str(value).strip().lower().replace("-", "_").split())


def _add(
    cfg: NumericalSelectionConfig,
    family: str,
    value: str,
    mapping: dict[str, str],
) -> None:
    key = mapping.get(_normalise(value))
    if key is not None:
        cfg.entries.append(NumericalSelection(family=family, configuration_key=key))
        cfg.required_families.append(family)


def build_numerical_selection(numerics: NumericalScheme) -> NumericalSelectionConfig:
    """Translate normalized adapter numerics into canonical registry keys.

    Only values with a defined one-to-one mapping are emitted.  An unmapped
    source value is therefore never silently converted into a different
    numerical method.
    """
    cfg = NumericalSelectionConfig()

    _add(cfg, "gradient", numerics.gradient_operator, _GRADIENT_KEYS)
    _add(cfg, "convection", numerics.momentum_scheme, _CONVECTION_KEYS)
    _add(cfg, "interpolation", numerics.momentum_interpolation, _INTERPOLATION_KEYS)

    if _normalise(numerics.transient_scheme) not in ("", "steady"):
        _add(cfg, "temporal", numerics.transient_scheme, _TEMPORAL_KEYS)

    _add(cfg, "pressure_velocity", numerics.coupled_solver, _PRESSURE_VELOCITY_KEYS)

    if numerics.preconditioner:
        _add(cfg, "preconditioner", numerics.preconditioner, _PRECONDITIONER_KEYS)

    # De-duplicate required families while preserving deterministic order.
    cfg.required_families = list(dict.fromkeys(cfg.required_families))
    return cfg
