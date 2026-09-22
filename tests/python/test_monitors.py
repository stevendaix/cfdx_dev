from cfdx.monitors import MonitorSample, MonitorSeries

def test_monitor_series_is_iteration_synchronized():
    series=MonitorSeries("pressure",[])
    series.append(MonitorSample(1,0.1,{"p":2.0}))
    series.append(MonitorSample(2,0.2,{"p":2.1}))
    assert series.at(2).time==0.2
