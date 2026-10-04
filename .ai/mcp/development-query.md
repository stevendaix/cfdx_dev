# Development MCP Query Contract

## Purpose

Expose the reproducible CFDX code-intelligence index through a narrow read-only MCP surface without changing the evidence semantics of the underlying tools.

## Operations

- index.validate: validate repository revision and indexed file hashes.
- code.symbol: query indexed symbols by name/qualified-name pattern.
- code.references: query indexed references by target symbol pattern.
- code.dependencies: query indexed dependency edges by path pattern.
- evidence.test: query discovered tests and explicit test-to-symbol relations.
- evidence.validation: query validation cases and explicit validation-to-symbol relations.

## Required response fields

Each operation must return:
- operation identifier;
- index freshness status;
- repository/index revision when available;
- deterministic result rows;
- explicit errors, if any.

A stale index is an error condition for evidence use. An empty result is not proof of absence.

## Safety

The initial surface is READ_ONLY. It must not mutate repository files, Git history, CFDX cases, build trees, or runtime state. Expensive execution and all mutation belong to separately permissioned MCP surfaces.

## Evidence limitations

Indexed relationships are not execution traces. Test-to-symbol and validation-to-symbol links do not establish test coverage, pass/fail status, verification, or qualification. Agents must confirm source, tests, documentation, and current CI/V&V evidence before drawing conclusions.

## Reference implementation

.ai/scripts/query_index.py is the current executable reference. An MCP adapter must preserve its deterministic query semantics and freshness checks.
