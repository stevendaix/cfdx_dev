from pathlib import Path

import pytest
from cfdx.probe_validation import ProbeSample, ProbeSeries
from cfdx.results_series import ResultFrame, ResultSeries

from cfdx import CFDXSession
from cfdx.application import (
    Application,
    DisplayObject,
    ResultDataset,
    ResultsChanged,
    build_results_state,
)


def sample_series(tmp_path: Path) -> ResultSeries:
    return ResultSeries(
        (
            ResultFrame(
                tmp_path / "result_000.vtu",
                0,
                0.0,
                True,
                ("pressure", "velocity"),
                0,
                "metadata",
            ),
            ResultFrame(
                tmp_path / "result_001.vtu", 1, 0.1, True, ("pressure",), 1, "metadata"
            ),
        )
    )


def test_results_snapshot_exposes_stable_frames_and_fields(tmp_path: Path) -> None:
    state = build_results_state(sample_series(tmp_path), directory=tmp_path)

    assert state.selected_frame_id == "result:0"
    assert state.selected_frame.path.endswith("result_000.vtu")
    assert state.field_names == ("pressure", "velocity")


def test_application_routes_result_selection_as_results_changed(tmp_path: Path) -> None:
    application = Application(CFDXSession())
    events = []
    application.events.subscribe(ResultsChanged, events.append)
    application.set_results(sample_series(tmp_path), directory=tmp_path)

    state = application.select_result_frame("result:1")
    state = application.select_result_field("pressure")

    assert state.results.selected_frame_id == "result:1"
    assert state.results.selected_field == "pressure"
    assert len(events) == 3

    with pytest.raises(KeyError):
        application.select_result_field("temperature")


def test_result_dataset_is_renderer_independent(tmp_path: Path) -> None:
    dataset = ResultDataset(
        stable_id="dataset:0",
        frame_id="result:0",
        path=str(tmp_path / "result_000.vtu"),
        fields=("p", "U"),
        time=0.25,
        iteration=12,
    )
    assert dataset.field_available("p")
    assert not dataset.field_available("T")


def test_display_object_round_trips_without_gui_dependencies() -> None:
    obj = DisplayObject(
        stable_id="display:pressure",
        kind="contour",
        dataset_id="dataset:0",
        field="p",
        parameters=(("levels", "12"), ("opacity", "0.75")),
    )
    assert DisplayObject.from_dict(obj.to_dict()) == obj
    assert obj.parameter_map() == {"levels": "12", "opacity": "0.75"}


def test_display_object_rejects_unknown_renderer_kind() -> None:
    with pytest.raises(ValueError, match="unsupported display object kind"):
        DisplayObject(stable_id="display:x", kind="vtk_actor", dataset_id="dataset:0")


def test_application_probe_csv_contract_round_trips_solver_history(
    tmp_path: Path,
) -> None:
    application = Application(CFDXSession())
    original = (
        ProbeSeries(
            "pressure",
            (0.5, 0.5, 0.5),
            (ProbeSample(1, 0.1, 101325.0), ProbeSample(2, 0.2, 101300.0)),
            field="p",
            unit="Pa",
        ),
    )
    path = application.export_probe_results(tmp_path / "case.probes.csv", original)

    assert path.is_file()
    loaded = application.read_probe_results(path)
    assert loaded == original
