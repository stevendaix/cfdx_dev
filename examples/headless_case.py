"""Minimal headless CFDX orchestration example."""
from __future__ import annotations

import sys
from pathlib import Path

from cfdx import CFDXSession, ExecutionController, SolverRunner


def main() -> int:
    script = Path(__file__).with_name("mock_solver.py")
    session = CFDXSession()
    controller = ExecutionController(session, SolverRunner([sys.executable, str(script)]))
    controller.start()
    assert controller.runner._thread is not None
    controller.runner._thread.join(timeout=10)
    print(controller.session.state.value)
    print(f"iteration={session.iteration} time={session.time:g}")
    return 0 if session.state.value == "CONVERGED" else 1


if __name__ == "__main__":
    raise SystemExit(main())
