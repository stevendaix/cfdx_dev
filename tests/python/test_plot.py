import importlib.util

import pytest

from cfdx.metrics import SolverMetrics
from cfdx.plot import MetricsHistory, create_metrics_plot


def test_metrics_history_is_bounded() -> None:
    history = MetricsHistory(max_points=2)
    for i in range(3):
        history.append(SolverMetrics(iteration=i, residuals=(("p", 1.0 / (i + 1)),)), "p")
    assert history.iterations == [1, 2]
    assert len(history.residuals) == 2


def test_plot_has_clear_optional_dependency_error() -> None:
    if importlib.util.find_spec("pyqtgraph") is not None:
        pytest.skip("pyqtgraph is installed")
    with pytest.raises(RuntimeError, match="pyqtgraph is required"):
        create_metrics_plot(MetricsHistory())
