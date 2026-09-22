import pytest
from cfdx.probe_validation import ProbeSample, ProbeSeries

def test_probe_series_validates_quantitatively():
    series=ProbeSeries("pressure",(0.5,0.5,0.5),(
        ProbeSample(1,0.1,1.0),ProbeSample(2,0.2,1.1)))
    assert series.validate_against((1.0,1.1),rtol=1e-12)==0.0

def test_probe_series_rejects_out_of_tolerance():
    series=ProbeSeries("pressure",(0.5,), (ProbeSample(1,0.1,1.2),))
    with pytest.raises(ValueError,match="exceeds tolerance"):
        series.validate_against((1.0,),rtol=1e-8)

def test_probe_series_rejects_non_monotonic_samples():
    series=ProbeSeries("pressure",(0.5,),(
        ProbeSample(2,0.2,1.0),ProbeSample(1,0.1,1.0)))
    with pytest.raises(ValueError,match="monotonic"):
        series.validate_against((1.0,1.0))
