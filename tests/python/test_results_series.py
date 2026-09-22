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