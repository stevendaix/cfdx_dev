# Development MCP Query Contract

## Purpose

Expose the reproducible CFDX code-intelligence index through a narrow read-only MCP surface without changing the evidence semantics of the underlying tools.

## Repository operations

- `repository.file_structure`: list tracked repository files below an optional relative path.
- `repository.status`: inspect the current Git branch/status without modifying repository state.

## Documentation operations

- `documentation.search`: search tracked repository text by pattern, optionally scoped to a repository-relative path.
- `documentation.read`: read a tracked repository file at `HEAD`.

Both operations are descriptive only and never modify repository state.

## Operations

- `index.validate`: validate repository revision and indexed file hashes.
- `code.symbol`: query indexed symbols by name/qualified-name pattern.
- `code.references`: query indexed references by target symbol pattern.
- `code.dependencies`: query indexed dependency edges by path pattern.
- `evidence.test`: query discovered tests and explicit test-to-symbol relations.
- `evidence.validation`: query validation cases and explicit validation-to-symbol relations.

## Server configuration

The server is bound to one repository/index pair at startup. Configuration can be supplied with:

- `--root` / `CFDX_AI_ROOT`: repository root.
- `--index` / `CFDX_AI_INDEX`: SQLite index path.

The MCP tool inputs therefore remain narrow: `index.validate` has no arguments and the five query operations accept only `pattern`.

## Required response fields

Each operation returns:

- operation identifier;
- index freshness status;
- repository revision when available;
- deterministic result rows;
- explicit errors, if any.

A stale index is an error condition for evidence use. Query operations return no indexed rows when the index is stale. An empty result on a valid index is not proof of repository absence.

## Permission

The first server surface is `READ_ONLY`. All six tools are annotated as read-only and do not modify repository files, Git history, CFDX cases, build trees, or runtime state. Development/maintainer mutation tools are separate future surfaces.

Tool annotations are advisory metadata, not a security boundary; host policy must enforce permissions independently.

## Transport

The reference adapter uses the official MCP Python SDK 2.3.0 and supports the standard `stdio` transport by default. Streamable HTTP is available for controlled local/server deployments. No SSE-specific implementation is added.

## Reference implementation

`.ai/mcp/development_server.py` is the MCP adapter. `.ai/scripts/query_index.py` remains the executable evidence/query reference underneath it.

The adapter validates freshness before every evidence query and converts the deterministic tabular CLI results into structured MCP output.

## Evidence limitations

Indexed relationships are not execution traces. Test-to-symbol and validation-to-symbol links do not establish test coverage, pass/fail status, verification, or qualification. Agents must confirm source, tests, documentation, and current CI/V&V evidence before drawing conclusions.

## Testing

`.ai/mcp/test_development_server.py` uses the MCP SDK's in-process client/server path to verify:

1. the complete six-tool surface is exposed;
2. tools are read-only annotated;
3. stale indexes are explicitly reported;
4. stale evidence queries do not present indexed rows.

Transport-level smoke testing remains a separate concern; the contract tests do not require a network port.
