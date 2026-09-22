from pathlib import Path
from cfdx.results_series import discover_result_series

def test_result_series_discovers_and_orders_vtk_files(tmp_path):
    for name in ("10.vtu","2.vtu","1.vtu","scratch.txt","latest.tmp"): (tmp_path/name).write_text("")
    series=discover_result_series(tmp_path)
    assert [p.name for p in series.paths]==["1.vtu","2.vtu","10.vtu"]
    assert series.is_transient
    assert series.frame(1).time==2.0

def test_result_series_ignores_unsupported_files(tmp_path):
    (tmp_path/"solution.vtu").write_text("")
    (tmp_path/"solution.csv").write_text("")
    series=discover_result_series(tmp_path)
    assert len(series.frames)==1
    assert series.frame(0).time is None

def test_empty_result_files_are_retained_as_incomplete(tmp_path):
    from cfdx.results_series import discover_result_series
    (tmp_path/"step_0.vtu").write_text("",encoding="utf-8")
    (tmp_path/"step_1.vtu").write_text("<dummy>",encoding="utf-8")
    series=discover_result_series(tmp_path)
    assert [f.path.name for f in series.frames]==["step_0.vtu", "step_1.vtu"]
    assert series.frames[0].complete is False
    assert series.frames[1].complete is True


def test_unreadable_latest_frame_is_marked_incomplete(tmp_path):
    (tmp_path / "step_0.vtu").write_text("<not-a-vtk-file>", encoding="utf-8")
    series = discover_result_series(tmp_path, inspect_fields=True)
    assert len(series.frames) == 1
    assert series.frames[0].complete is False


def test_zero_length_result_is_reported_as_incomplete(tmp_path):
    path = tmp_path / "2.vtu"
    path.write_bytes(b"")
    series = discover_result_series(tmp_path)
    assert len(series.frames) == 1
    assert series.frames[0].complete is False
