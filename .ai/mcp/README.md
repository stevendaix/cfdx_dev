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

## Read-only Development MCP repository surface

The read-only Development MCP also exposes repository inspection without modifying Git state:

| Operation | Input | Result |
|---|---|---|
| `repository.file_structure` | optional relative path | tracked repository files below that path |
| `repository.status` | none | current Git branch/status lines |

The file-structure operation uses tracked files only, and rejects paths that escape the repository. These operations are descriptive and do not imply build, test, validation, or qualification evidence.

## Read-only Development MCP documentation surface

The read-only Development MCP exposes repository-grounded documentation lookup without modifying Git state:

| Operation | Input | Result |
|---|---|---|
| `documentation.search` | pattern + optional relative path | tracked text matches with file/line context |
| `documentation.read` | tracked relative path | file content at `HEAD` |

Both operations reject repository-escaping paths. `documentation.read` only reads tracked files through Git and therefore does not expose arbitrary filesystem paths.

## Read-only Development MCP query surface

The first implemented Development MCP surface is deliberately read-only. It maps to the reproducible index/query layer and preserves the evidence contract.

| Operation | Input | Result |
|---|---|---|
| `index.validate` | no tool arguments; repository/index are server configuration | freshness status and explicit stale reasons |
| `code.symbol` | pattern | matching symbols with paths/source locations |
| `code.references` | pattern | indexed references with source locations/targets |
| `code.dependencies` | pattern | indexed dependency edges |
| `evidence.test` | pattern | discovered tests plus explicit test-to-symbol relations |
| `evidence.validation` | pattern | validation cases plus explicit validation-to-symbol relations |

Every response identifies freshness before presenting indexed evidence. Stale indexes are never presented as current. Empty valid queries mean “no indexed result found”; they do not prove repository absence.

Evidence relations are descriptive only. An explicit test-to-symbol or validation-to-symbol relation is not execution coverage, pass/fail evidence, verification, or qualification evidence.

## Adapter

`.ai/mcp/development_server.py` implements the ten read-only operations using the official MCP Python SDK 2.3.0. `.ai/scripts/query_index.py` remains the executable reference for deterministic index semantics.

The adapter is configured with `--root` / `--index` or the `CFDX_AI_ROOT` / `CFDX_AI_INDEX` environment variables. stdio is the default transport; Streamable HTTP is available for controlled deployments.

The adapter performs freshness validation before every indexed evidence query and returns structured MCP results. It does not modify repository, build, case, or runtime state.

## Read-only Runtime MCP artifact surface

The first Runtime MCP surface is deliberately limited to inspection of existing CFDX HDF5 artifacts. It does not create cases, modify files, launch solvers, monitor processes, restart simulations, or alter runtime state.

| Operation | Input | Result |
|---|---|---|
| `case.inspect` | root-relative `.cfdx.h5` path | case attributes, parsed setup configuration, dataset metadata, runtime-state presence |
| `checkpoint.inspect` | root-relative `.dat.h5` path | checkpoint attributes, cell-ID presence, field metadata (names/shapes/dtypes/sizes) |

The Runtime MCP root is configured with `--root` or `CFDX_RUNTIME_ROOT`. Paths are required to be relative to that root and repository/filesystem escape attempts are rejected. Field arrays are not returned by the inspection surface, and no solver executable is invoked.

The adapter is `.ai/mcp/runtime_server.py`; its contract is exercised by `.ai/mcp/test_runtime_server.py`. All tools are annotated read-only and closed-world. This surface is an inspection foundation, not evidence that a simulation is valid, converged, verified, or qualified.

## Permissions

All operations in this first surface are `READ_ONLY`. Tool annotations communicate the intended behavior to MCP hosts, but they are not a security boundary. Write/execute operations belong to separately permissioned Development or Runtime MCP surfaces.

## Testing

`.ai/mcp/test_development_server.py` exercises the server in-process through the MCP client and verifies the ten-tool surface, path-safety checks, stale-index reporting, and refusal to present stale evidence.

The first phase specified architecture only; this adapter is the first executable Development MCP component. Transport-specific deployment and broader MCP integration tests remain future work.
