# CFDX GUI and TUI

The application layer is shared by headless, TUI and GUI clients. Python owns
orchestration state; numerical kernels remain in C++.

## Headless

Run the deterministic orchestration example:

    python examples/headless_case.py

## GUI

Install the optional visualization stack:

    pip install -e ".[gui]"

Then launch:

    cfdx-gui

The GUI is intentionally optional. A headless installation does not require
PySide6, PyVista or Qt.

## State and execution

CFDXSession is the single source of truth for case state, revisions,
iteration/time and checkpoints. SolverRunner consumes stdout and stderr on
background threads. SolverMetricsParser extracts iteration, physical time,
CFL and named residuals.

The GUI must never duplicate these values in a second model. Views observe the
session and execution callbacks.

## 3D results

PyVistaRenderer is an optional adapter for VTU and other formats understood by
PyVista. The renderer contract is also implemented by NullRenderer so
post-processing code can be tested without a display server.

## CI/headless use

All core orchestration and renderer-contract tests are designed to run without
a graphical display. Optional Qt/PyVista tests are skipped when those
dependencies are not installed.
