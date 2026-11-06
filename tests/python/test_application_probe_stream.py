from cfdx.probe_stream import ProbeStreamSample

from cfdx import CFDXSession
from cfdx.application import (
    Application,
    ApplicationStateChanged,
    ProbeSampleReceived,
)


def test_probe_samples_use_a_dedicated_application_event() -> None:
    application = Application(CFDXSession())
    state_events = []
    sample_events = []
    application.events.subscribe(ApplicationStateChanged, state_events.append)
    application.events.subscribe(ProbeSampleReceived, sample_events.append)

    application.events.publish(
        ProbeSampleReceived(
            application.state,
            ProbeStreamSample("pressure", "p", "Pa", 1, 0.0, 101325.0),
        )
    )

    assert state_events == []
    assert len(sample_events) == 1
    assert sample_events[0].sample.name == "pressure"
