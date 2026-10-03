from pathlib import Path

import pytest

from cfdx import CFDXSession
from cfdx.application import Application, ResultsChanged, build_results_state
from cfdx.results_series import ResultFrame, ResultSeries


def sample_series(tmp_path: Path) -> ResultSeries:
    return ResultSeries(
        (
            ResultFrame(tmp_path / "result_000.vtu", 0, 0.0, True, ("pressure", "velocity"), 0, "metadata"),
            ResultFrame(tmp_path / "result_001.vtu", 1, 0.1, True, ("pressure",), 1, "metadata"),
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
