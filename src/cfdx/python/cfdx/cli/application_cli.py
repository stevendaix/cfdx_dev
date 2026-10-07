"""Headless CLI using the same Application contracts as the GUI."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

from ..application.application import Application
from ..execution import ExecutionController
from ..runner import SolverRunner
from ..tui import TuiRenderer


def _application(case: Path, dat: Path | None = None) -> Application:
    app = Application()
    app.open_project(case, with_dat=dat is not None, dat_path=dat)
    return app


def cmd_status(args: argparse.Namespace) -> int:
    app = _application(args.case, args.dat)
    print(TuiRenderer.render(app.session))
    return 0


def cmd_info(args: argparse.Namespace) -> int:
    app = _application(args.case, args.dat)
    state = app.state\n    print(json.dumps({\n        "project_path": state.project.path,\n        "simulation_state": state.simulation_state.value,\n        "iteration": state.execution.iteration,\n        "time": state.execution.time,\n        "requires_restart": state.requires_restart,\n    }, indent=2, default=str))
    return 0


def cmd_results(args: argparse.Namespace) -> int:
    app = _application(args.case)
    state = app.open_results(args.directory)
    print(json.dumps({
        "directory": str(args.directory),
        "frames": [frame.to_dict() for frame in state.results.frames],
        "fields": list(state.results.field_names),
    }, indent=2, default=str))
    return 0


def _controller(app: Application, args: argparse.Namespace) -> ExecutionController:
    command = [str(args.solver), str(args.case)]
    if args.solver_args:
        command.extend(args.solver_args)
    controller = ExecutionController(
        app.session,
        SolverRunner(command, cwd=args.case.parent),
    )
    app.attach_controller(controller)
    return controller


def _wait(controller: ExecutionController, timeout: float | None) -> int:
    start = time.monotonic()
    while controller.runner.running:
        if timeout is not None and time.monotonic() - start > timeout:
            controller.stop()
            raise TimeoutError(f"solver exceeded timeout of {timeout:g} s")
        time.sleep(0.05)
    return 0 if controller.session.state.value == "CONVERGED" else 1


def cmd_run(args: argparse.Namespace) -> int:
    app = _application(args.case, args.dat)
    controller = _controller(app, args)
    if args.dat is not None:
        controller.restart(args.dat)
    else:
        app.run()
    code = _wait(controller, args.timeout)
    print(TuiRenderer.render(app.session))
    if controller.error is not None:
        print(f"CFDX ERROR: solver exited with {controller.error.returncode}")
    return code


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="cfdx", description="CFDX application CLI")
    sub = parser.add_subparsers(dest="command", required=True)

    for name, func in (("status", cmd_status), ("info", cmd_info)):
        command = sub.add_parser(name)
        command.add_argument("case", type=Path)
        command.add_argument("--dat", type=Path)
        command.set_defaults(func=func)

    results = sub.add_parser("results")
    results.add_argument("case", type=Path)
    results.add_argument("directory", type=Path)
    results.set_defaults(func=cmd_results)

    run = sub.add_parser("run")
    run.add_argument("case", type=Path)
    run.add_argument("--solver", type=Path, required=True)
    run.add_argument("--dat", type=Path)
    run.add_argument("--timeout", type=float)
    run.add_argument("solver_args", nargs=argparse.REMAINDER)
    run.set_defaults(func=cmd_run)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return int(args.func(args))
    except (OSError, RuntimeError, TimeoutError, ValueError) as exc:
        print(f"CFDX ERROR: {exc}")
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
