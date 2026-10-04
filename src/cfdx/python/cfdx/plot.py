"""Optional pyqtgraph monitor plot widget."""
from __future__ import annotations

from .metrics import SolverMetrics


class MetricsHistory:
    """Bounded in-memory history independent of any GUI toolkit."""

    def __init__(self, max_points: int = 2000) -> None:
        if max_points <= 0:
            raise ValueError("max_points must be positive")
        self.max_points = max_points
        self.iterations: list[int] = []
        self.residuals: list[float] = []

    def append(self, metrics: SolverMetrics, residual_name: str) -> None:
        value = metrics.residual(residual_name)
        if metrics.iteration is None or value is None:
            return
        self.iterations.append(metrics.iteration)
        self.residuals.append(value)
        del self.iterations[:-self.max_points]
        del self.residuals[:-self.max_points]


def create_metrics_plot(history: MetricsHistory):
    """Create a pyqtgraph PlotWidget from history, importing it lazily."""
    try:
        import pyqtgraph as pg
    except ImportError as exc:
        raise RuntimeError("pyqtgraph is required for metrics plots") from exc
    plot = pg.PlotWidget()
    plot.plot(history.iterations, history.residuals)
    plot.setLabel("bottom", "Iteration")
    plot.setLabel("left", "Residual")
    return plot
