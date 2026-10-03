# CFDX MCP Architecture

Two logical MCP families are planned.

## Development MCP

Operations over source code, symbols, dependencies, documentation, tests, validation, CI and GitHub. Write and execute operations use explicit permission profiles.

## Runtime MCP

Operations over CFDX itself: create/configure cases, inspect cases, mesh/import, build, run, monitor, restart, extract fields/probes, inspect convergence/conservation, compare references, validate and generate reports.

## Tool design rules

1. Tools have narrow, explicit responsibilities.
2. Inputs and outputs are schema-defined.
3. Errors are returned explicitly.
4. Read-only tools do not mutate repository or case state.
5. Expensive/destructive operations require policy-controlled permission.
6. MCP never claims success when the underlying operation failed.
7. Development and runtime permissions remain separable.

The first phase specifies architecture only; the servers are later phases.