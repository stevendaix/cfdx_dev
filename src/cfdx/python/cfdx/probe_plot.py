"""Optional interactive plot adapter for native point-probe histories."""

from __future__ import annotations

from collections.abc import Sequence

from .probe_validation import ProbeSeries
from .theme import DEFAULT_THEME, style_scientific_plot, theme_tokens


def _pyqtgraph():
    try:
        import pyqtgraph as pg
    except ImportError as exc:
        raise RuntimeError(
            "pyqtgraph is required for interactive probe plots; "
            'install the optional GUI dependencies with `pip install "cfdx[gui]"`'
        ) from exc
    return pg


def create_probe_plot():
    """Create a pan/zoom plot widget without importing Qt at module import time."""
    pg = _pyqtgraph()
    plot = pg.PlotWidget()
    plot.setObjectName("workbench.probes.history_plot")
    style_scientific_plot(plot, DEFAULT_THEME)
    plot.showGrid(x=True, y=True, alpha=0.25)
    plot.setMouseEnabled(x=True, y=True)
    plot.setMenuEnabled(True)
    plot.setClipToView(True)
    plot.setDownsampling(auto=True, mode="peak")
    plot.setLabel("bottom", "Time", units="s")
    plot.setLabel("left", "Value")
    return plot


def update_probe_plot(
    plot,
    series: Sequence[ProbeSeries],
    *,
    x_axis: str = "time",
    theme: str = DEFAULT_THEME,
) -> None:
    """Redraw selected probe series, keeping unit groups separate in the caller."""
    pg = _pyqtgraph()
    if x_axis not in {"time", "iteration"}:
        raise ValueError("x_axis must be 'time' or 'iteration'")
    tokens = theme_tokens(theme)
    colors = (
        tokens["accent.primary"], tokens["accent.secondary"],
        tokens["status.info"], tokens["status.warning"], tokens["status.error"],
    )
    plot.clear()
    plot.addLegend(offset=(10, 10))
    if x_axis == "time":
        plot.setLabel("bottom", "Time", units="s")
    else:
        plot.setLabel("bottom", "Iteration")
    if series:
        first = series[0]
        field_label = first.field or "Probe value"
        unit_label = first.unit or ""
        plot.setLabel("left", field_label, units=unit_label or None)
    else:
        plot.setLabel("left", "Value")
    for index, probe in enumerate(series):
        x_values = [
            sample.time if x_axis == "time" else sample.iteration
            for sample in probe.samples
        ]
        y_values = [sample.value for sample in probe.samples]
        if not x_values:
            continue
        plot.plot(
            x_values,
            y_values,
            pen=pg.mkPen(colors[index % len(colors)], width=2),
            name=probe.name,
        )
    plot.enableAutoRange()
    plot.autoRange()
