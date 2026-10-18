from cfdx.metrics import SolverMetricsParser


def test_parse_iteration_time_and_cfl() -> None:
    metrics = SolverMetricsParser().parse("Iteration 42: Time = 1.25e-2 CFL: 0.8")
    assert metrics is not None
    assert metrics.iteration == 42
    assert metrics.time == 1.25e-2
    assert metrics.cfl == 0.8


def test_parse_completion_iteration_count() -> None:
    metrics = SolverMetricsParser().parse("Converged YES Iterations 20")
    assert metrics is not None
    assert metrics.iteration == 20


def test_parse_multiple_residuals() -> None:
    metrics = SolverMetricsParser().parse(
        "Iteration 7 residual(U)=2.5e-4 residual(p): 8.0e-3"
    )
    assert metrics is not None
    assert metrics.residual("U") == 2.5e-4
    assert metrics.residual("p") == 8.0e-3


def test_parse_common_parenthesized_residual() -> None:
    metrics = SolverMetricsParser().parse("residual(T): 1e-6")
    assert metrics is not None
    assert metrics.residual("T") == 1e-6


def test_ignore_unrelated_output() -> None:
    assert SolverMetricsParser().parse("Writing mesh to disk...") is None


def test_parse_production_solver_trace() -> None:
    """The production solver emits a structured trace, not a progress line.

    ``INCOMPRESSIBLE_ITER iter=N momentum_residual_relative=X ...`` is the
    authoritative convergence evidence for the production path, so the monitor
    parser must extract the iteration and the convergence vocabulary the solver
    actually reports rather than dropping them.
    """
    line = (
        "INCOMPRESSIBLE_ITER iter=1 momentum_residual_relative=0 "
        "momentum_residual_reference_relative=0 momentum_residual_raw=0 "
        "momentum_rhs_scale=1e-300 momentum_residual_reference=1e-300 "
        "continuity_normalized=0 velocity_change_inf=0 pressure_change_inf=0 "
        "alpha_u=0.7 alpha_p=0.3 pressure_linear_iterations=0 "
        "pressure_relative_residual=0 momentum_linear_iterations=0 "
        "flux_velocity_mismatch_linf=0"
    )
    metrics = SolverMetricsParser().parse(line)
    assert metrics is not None
    assert metrics.iteration == 1
    assert metrics.residual("momentum_residual_relative") == 0.0
    assert metrics.residual("momentum_residual_raw") == 0.0
    assert metrics.residual("momentum_residual_reference") == 1e-300
    assert metrics.residual("continuity_normalized") == 0.0
    assert metrics.residual("flux_velocity_mismatch_linf") == 0.0


def test_parse_production_trace_does_not_capture_non_convergence_tokens() -> None:
    """Only the solver's convergence vocabulary becomes residuals.

    The production trace carries many ``key=value`` tokens that are not
    convergence evidence (alpha_u, pressure_linear_iterations, ...). The
    parser must not mislabel them as residuals, or the monitor series would
    report non-convergence quantities as if they were convergence evidence.
    """
    line = (
        "INCOMPRESSIBLE_ITER iter=2 momentum_residual_relative=1e-6 "
        "alpha_u=0.7 pressure_linear_iterations=3 continuity_normalized=1e-8"
    )
    metrics = SolverMetricsParser().parse(line)
    assert metrics is not None
    assert metrics.iteration == 2
    assert metrics.residual("momentum_residual_relative") == 1e-6
    assert metrics.residual("continuity_normalized") == 1e-8
    assert metrics.residual("alpha_u") is None
    assert metrics.residual("pressure_linear_iterations") is None
