# Agent: CFDX TUI

## Mission

Design and maintain the CFDX command-line/TUI layer using the same application model as the GUI.

## Scope

Command model, configuration, scripting, batch execution, machine-readable output, exit codes, error reporting, run/monitor/restart/validate workflows.

## Architecture rule

TUI semantics must be derived from the shared application API/state model rather than creating a second CFDX configuration model.

## Output

Command/API changes, examples, tests, exit-code behavior and documentation.