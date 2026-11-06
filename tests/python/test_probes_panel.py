import importlib.util
import json

import pytest
from cfdx.probe import Probe, ProbeField
from cfdx.probe_stream import ProbeStreamSample
from cfdx.probe_validation import ProbeSample, ProbeSeries

from cfdx import CFDXSession
from cfdx.application import Application, ProbeSampleReceived

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


def test_live_probe_samples_update_interactive_plot_from_worker_thread() -> None:
    pytest.importorskip("pyqtgraph")
    from threading import Thread

    from cfdx.probes_panel import ProbesPanel
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication(["cfdx-live-probes-test"])
    session = CFDXSession()
    session.case.probes.append(Probe("pressure", 0.5, 0.0, 0.0, ProbeField.PRESSURE))
    application = Application(session)
    panel = ProbesPanel(application)

    def publish(iteration: int, value: float) -> None:
        application.events.publish(
            ProbeSampleReceived(
                application.state,
                ProbeStreamSample(
                    "pressure", "p", "Pa", iteration, iteration / 10, value
                ),
            )
        )

    for iteration, value in ((1, 101325.0), (2, 101326.0)):
        thread = Thread(target=publish, args=(iteration, value))
        thread.start()
        thread.join(timeout=2)
        assert not thread.is_alive()
        app.processEvents()
    panel._flush_live_plot()

    assert panel._results_summary.text() == (
        "Live · 1 probe · 2 buffered samples · iterations 1–2"
    )
    assert len(panel._probe_plot.plotItem.listDataItems()) == 1
    x_values, y_values = panel._probe_plot.plotItem.listDataItems()[0].getData()
    assert list(x_values) == pytest.approx([1, 2])
    assert list(y_values) == pytest.approx([101325.0, 101326.0])

    # The live axis selector must redraw immediately and retain the chosen
    # physical-time axis across later batched stream updates.
    panel._plot_axis_selector.setCurrentIndex(
        panel._plot_axis_selector.findData("time")
    )
    x_values, _y_values = panel._probe_plot.plotItem.listDataItems()[0].getData()
    assert list(x_values) == pytest.approx([0.1, 0.2])

    thread = Thread(target=publish, args=(3, 101327.0))
    thread.start()
    thread.join(timeout=2)
    assert not thread.is_alive()
    app.processEvents()
    panel._flush_live_plot()
    x_values, y_values = panel._probe_plot.plotItem.listDataItems()[0].getData()
    assert list(x_values) == pytest.approx([0.1, 0.2, 0.3])
    assert list(y_values) == pytest.approx([101325.0, 101326.0, 101327.0])

    panel._plot_axis_selector.setCurrentIndex(
        panel._plot_axis_selector.findData("iteration")
    )
    x_values, _y_values = panel._probe_plot.plotItem.listDataItems()[0].getData()
    assert list(x_values) == pytest.approx([1, 2, 3])
    panel.close()
    app.processEvents()


def test_live_probe_stream_exports_buffered_samples_as_csv_or_json(
    tmp_path, monkeypatch
) -> None:
    pytest.importorskip("pyqtgraph")
    from cfdx.probes_panel import ProbesPanel
    from PySide6.QtWidgets import QApplication, QFileDialog

    app = QApplication.instance() or QApplication(["cfdx-live-probe-export-test"])
    session = CFDXSession()
    session.case.probes.append(Probe("pressure", 0.5, 0.0, 0.0, ProbeField.PRESSURE))
    application = Application(session)
    panel = ProbesPanel(application)
    assert not panel._export_live_button.isEnabled()

    for iteration, value in ((1, 101325.0), (2, 101326.0)):
        application.events.publish(
            ProbeSampleReceived(
                application.state,
                ProbeStreamSample(
                    "pressure", "p", "Pa", iteration, iteration / 10, value
                ),
            )
        )
    panel._flush_live_plot()
    assert panel._export_live_button.isEnabled()

    csv_path = tmp_path / "live_snapshot.csv"
    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(csv_path), "CFDX probe CSV (*.csv)"),
    )
    panel._export_live_stream()
    csv_series = application.read_probe_results(csv_path)
    assert len(csv_series) == 1
    assert csv_series[0].name == "pressure"
    assert [sample.iteration for sample in csv_series[0].samples] == [1, 2]

    json_path = tmp_path / "live_snapshot.json"
    monkeypatch.setattr(
        QFileDialog,
        "getSaveFileName",
        lambda *_args: (str(json_path), "CFDX probe JSON (*.json)"),
    )
    panel._export_live_stream()
    document = json.loads(json_path.read_text(encoding="utf-8"))
    assert document["source"] == "live-stream-buffer"
    assert document["probes"][0]["name"] == "pressure"
    assert [sample["iteration"] for sample in document["probes"][0]["samples"]] == [
        1,
        2,
    ]
    panel.close()
    app.processEvents()


def test_probe_comparison_overlays_same_probe_from_multiple_csv_runs(
    tmp_path, monkeypatch
) -> None:
    pytest.importorskip("pyqtgraph")
    from cfdx.probes_panel import ProbesPanel
    from PySide6.QtWidgets import QApplication, QFileDialog

    app = QApplication.instance() or QApplication(["cfdx-probe-comparison-test"])
    application = Application(CFDXSession())
    first_dir = tmp_path / "baseline"
    second_dir = tmp_path / "candidate"
    first_dir.mkdir()
    second_dir.mkdir()
    first_path = first_dir / "case.probes.csv"
    second_path = second_dir / "case.probes.csv"
    application.export_probe_results(
        first_path,
        (
            ProbeSeries(
                "pressure",
                (0.5, 0.0, 0.0),
                (ProbeSample(1, 0.1, 101325.0), ProbeSample(2, 0.2, 101326.0)),
                field="p",
                unit="Pa",
            ),
        ),
    )
    application.export_probe_results(
        second_path,
        (
            ProbeSeries(
                "pressure",
                (0.5, 0.0, 0.0),
                (ProbeSample(1, 0.1, 101300.0), ProbeSample(2, 0.2, 101301.0)),
                field="p",
                unit="Pa",
            ),
        ),
    )
    monkeypatch.setattr(
        QFileDialog,
        "getOpenFileNames",
        lambda *_args: ([str(first_path), str(second_path)], "CFDX probe CSV (*.csv)"),
    )
    panel = ProbesPanel(application)

    panel._compare_results()

    assert panel._results_summary.text() == (
        "Comparison · 2 files · 2 probes · 4 samples · iterations 1–2"
    )
    assert {probe.name for probe in panel._display_probe_series} == {
        "case.probes / pressure",
        "case.probes (2) / pressure",
    }
    assert len(panel._probe_plot.plotItem.listDataItems()) == 2
    assert panel._export_results_button.isEnabled()
    panel.close()
    app.processEvents()
