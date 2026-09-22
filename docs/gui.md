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


## Case files

The GUI exposes the same persistence contract as the headless API:

- Read Case loads the HDF5 case without consuming a DAT restart artifact.
- Read Case + DAT validates the paired DAT SHA-256 before accepting it.
- Save Case writes the HDF5 source of truth.
- Save Case + DAT writes the HDF5 case and its canonical case.dat sibling.

There is no GUI-specific checkpoint format. A modified case is marked in the
window and the GUI asks whether to save before Run or Close.

## Execution controls

Run, Pause, Resume and Stop are routed through ExecutionController.
On POSIX systems, the underlying SolverRunner controls the solver process
group, so Pause/Resume and Stop affect MPI children as well as the parent
process. Windows explicitly reports pause/resume as unsupported rather than
simulating it with a state-only transition.

A Run is refused until the case is saved and a solver executable is configured.
This prevents launching a solver against an unsaved or stale case.

## 3D workspace

The GUI keeps a persistent case tree and provides a separate 3D Results
workspace. VTU/VTK datasets can be loaded directly into the optional
PyVista/Qt view. The visualization adapter remains separate from the solver
and can therefore be disabled on headless HPC installations.
