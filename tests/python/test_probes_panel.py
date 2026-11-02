import importlib.util

import pytest
from cfdx.probe_validation import ProbeSample, ProbeSeries

from cfdx import CFDXSession
from cfdx.application import Application

pytestmark = pytest.mark.skipif(
    importlib.util.find_spec("PySide6") is None, reason="PySide6 optional"
)


def test_probe_history_csv_load_and_export_uses_application_contract(
    tmp_path, monkeypatch
) -> None:
    pytest.importorskip("pyqtgraph")
    from cfdx.probes_panel import ProbesPanel
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QApplication, QFileDialog

    app = QApplication.instance() or QApplication(["cfdx-probes-panel-test"])
    application = Application(CFDXSession())
    panel = ProbesPanel(application)
    source = tmp_path / "run.probes.csv"
    destination = tmp_path / "export.probes.csv"
    application.export_probe_results(
        source,
        (
            ProbeSeries(
                name="outlet",
                point=(1.0, 0.0, 0.0),
                field="p",
                unit="Pa",
                samples=(
                    ProbeSample(iteration=2, time=0.1, value=101325.0),
                    ProbeSample(iteration=4, time=0.2, value=101326.0),
                ),
            ),
            ProbeSeries(
                name="inlet",
                point=(0.0, 0.0, 0.0),
                field="p",
                unit="Pa",
                samples=(
                    ProbeSample(iteration=2, time=0.1, value=101300.0),
                    ProbeSample(iteration=4, time=0.2, value=101301.0),
                ),
            ),
            ProbeSeries(
                name="velocity",
                point=(0.5, 0.0, 0.0),
                field="u_x",
                unit="m/s",
                samples=(
                    ProbeSample(iteration=2, time=0.1, value=1.0),
                    ProbeSample(iteration=4, time=0.2, value=1.5),
                ),
            ),
        ),
    )
    monkeypatch.setattr(
        QFileDialog, "getOpenFileName", lambda *_args: (str(source), "")
    )
    monkeypatch.setattr(
        QFileDialog, "getSaveFileName", lambda *_args: (str(destination), "")
    )

    panel._load_results()

    assert panel._results_summary.text() == "3 probes · 6 samples · iterations 2–4"
    assert panel._export_results_button.isEnabled()
    assert panel._plot_field_selector.isEnabled()
    assert panel._plot_field_selector.count() == 2
    assert len(panel._probe_plot.plotItem.listDataItems()) == 2
    panel._plot_field_selector.setCurrentIndex(1)
    assert len(panel._probe_plot.plotItem.listDataItems()) == 1
    panel._plot_field_selector.setCurrentIndex(0)
    panel._plot_probe_list.item(1).setCheckState(Qt.CheckState.Unchecked)
    assert len(panel._probe_plot.plotItem.listDataItems()) == 1
    panel._plot_axis_selector.setCurrentIndex(1)
    x_values, _y_values = panel._probe_plot.plotItem.listDataItems()[0].getData()
    assert list(x_values) == [2, 4]
    panel._export_results()
    assert application.read_probe_results(
        destination
    ) == application.read_probe_results(source)
    panel.close()
    app.processEvents()
