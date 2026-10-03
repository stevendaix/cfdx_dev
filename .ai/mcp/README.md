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

## Read-only Development MCP query contract

The first executable MCP-facing surface is deliberately read-only. It maps to the reproducible index/query layer and must preserve the evidence contract.

### Operations

| Operation | Required input | Result |
|---|---|---|
| `index.validate` | repository root, index path | freshness status and explicit stale reasons |
| `code.symbol` | pattern | matching symbols with repository paths and source locations |
| `code.references` | pattern | indexed references with source locations and target symbols |
| `code.dependencies` | pattern | indexed dependency edges |
| `evidence.test` | pattern | discovered tests plus explicit test-to-symbol relations |
| `evidence.validation` | pattern | validation cases plus explicit validation-to-symbol relations |

### Result contract

Every response must identify the index revision/freshness status before presenting indexed evidence. Stale indexes must not be presented as current. Empty query results mean “no indexed result found”; they do not prove repository absence.

Evidence relations are descriptive only. In particular, an explicit test-to-symbol or validation-to-symbol relation is not execution coverage, pass/fail evidence, or qualification evidence.

### Permission

All operations in this first surface are `READ_ONLY`. They must not modify repository, build, case, or runtime state. Development/maintainer mutation tools are separate future surfaces.

### Transport independence

This contract does not mandate a particular MCP SDK or transport. The CLI/index implementation remains the executable reference until a server adapter is introduced and tested against the same schemas and evidence semantics.
