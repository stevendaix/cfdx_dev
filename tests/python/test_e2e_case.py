from pathlib import Path
import sys

from cfdx import CFDXSession, ExecutionController, SolverRunner


def test_short_headless_case(tmp_path: Path) -> None:
    solver = tmp_path / "solver.py"
    output_dir = tmp_path / "results"
    solver.write_text(
        "import json\n"
        "from pathlib import Path\n"
        f"out = Path({str(output_dir)!r})\n"
        "out.mkdir(parents=True, exist_ok=True)\n"
        "history = {'format': 'CFDX-CONVERGENCE', 'schema_version': 1, "
        "'converged': True, 'iterations': 2, 'history': ["
        "{'iteration': 1, 'momentum': 1e-3}, "
        "{'iteration': 2, 'momentum': 1e-11}]}\n"
        "execution = {'format': 'CFDX-EXECUTION', 'schema_version': 1, "
        "'converged': True, 'iterations': 2, 'process_exit_code': 0, "
        "'convergence_reason': 'synthetic certified test run', "
        "'convergence_status': 1, 'artifacts': {'convergence_json': True}}\n"
        "(out / 'convergence.json').write_text(json.dumps(history), encoding='utf-8')\n"
        "(out / 'execution.json').write_text(json.dumps(execution), encoding='utf-8')\n"
        "print('Iteration 1 Time = 0.1 CFL: 0.2')\n"
        "print('Iteration 2 Time = 0.2 CFL: 0.2')\n",
        encoding="utf-8",
    )
    session = CFDXSession()
    session.case.execution.output_dir = str(output_dir)
    controller = ExecutionController(session, SolverRunner([sys.executable, str(solver)]))
    controller.start()
    assert controller.runner._thread is not None
    controller.runner._thread.join(timeout=10)
    assert session.state.value == "CONVERGED"
    assert controller.convergence_verdict is True
    assert controller.convergence_evidence is not None
    assert session.iteration == 2
    assert session.time == 0.2
