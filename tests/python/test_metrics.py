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
