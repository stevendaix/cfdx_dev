from cfdx import CFDXSession
from cfdx.application import RunCenterModel
from cfdx.metrics import SolverMetrics
from cfdx.runner import ProcessResult


class FakeController:
    def __init__(self) -> None:
        self.session = CFDXSession()
        self.latest_metrics = None
        self.error = None
        self.on_output = None
        self.on_metrics = None
        self.on_complete = None


def test_run_center_tracks_output_metrics_and_completion() -> None:
    controller = FakeController()
    model = RunCenterModel(controller, max_output_lines=2)
    changes = []
    model.on_change = changes.append

    controller.on_output("Iteration 4 Time = 0.25", False)
    controller.on_metrics(SolverMetrics(iteration=4, time=0.25, cfl=0.5))
    controller.on_output("solver warning", True)
    controller.on_output("last line", False)
    controller.on_complete(ProcessResult(0, ("solver",)))

    assert model.state.output == ("STDERR | solver warning", "last line")
    assert model.state.latest_metrics.cfl == 0.5
    assert len(changes) == 5
    model.close()
    assert controller.on_output is None
