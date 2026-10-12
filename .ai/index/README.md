# CFDX Code Intelligence Index

The code-intelligence layer provides reproducible structural queries for CFDX agents.

## Evidence contract

The index complements, but never replaces, repository source truth. Query flow is:

index result -> freshness check -> source confirmation -> relevant tests/docs -> decision

A missing result does not prove absence. Stale or partial coverage must be reported explicitly.

## Canonical model

- files — repository path, language and content identity;
- symbols — functions, methods, classes, structs, namespaces and other definitions;
- references — symbol-to-symbol references;
- dependencies — include/import/build relationships;
- tests — discovered test sources and implementation symbols when known;
- validation_cases — documented/executable V&V cases and implementation links when known.

SQLite is an implementation target for AI tooling, not part of the CFDX runtime ABI.

## Parser coverage

- Python: standard-library AST.
- C/C++ fallback: conservative declaration/include scanner, coarse navigation only.
- C/C++ semantic path: optional libclang using CMake `compile_commands.json` and the real translation-unit arguments.
- Parser provenance, coverage, parse failures and diagnostics are retained in `index_metadata`.

Semantic coverage is only valid for successfully parsed translation units. It does not establish whole-program completeness, complete macro semantics or complete virtual-dispatch analysis.

## Test and V&V inventory

`.ai/scripts/index_test_validation.py` adds a conservative inventory layer:

- test sources are discovered from `tests/` using explicit filename conventions;
- validation cases are discovered from `tests/validation/` directories and supported manifest files;
- discovery provenance is stored explicitly;
- no test-to-symbol or validation-to-symbol relationship is inferred by this layer.

This distinction is intentional: inventory is evidence of repository artifacts, not proof of coverage.

## Query tooling

`.ai/scripts/query_index.py` supports:

- `--symbol PATTERN`
- `--references PATTERN`
- `--dependency PATTERN`
- `--test PATTERN` — reports discovered test artifacts and only explicit test-to-symbol links;
- `--validation PATTERN` — reports validation cases and only explicit validation-to-symbol links;
- `--validate`

`--validate` compares indexed content hashes with current files and, where available, the repository revision. A non-zero result means the index must not be treated as current.

The query tool uses only Python's standard library and has no CFDX runtime dependency.

## Required properties

1. Reproducible.
2. Incrementally updateable.
3. Traceable to repository paths/source locations.
4. Version-aware.
5. Non-authoritative.
6. Deterministic.
7. Explicit about unsupported syntax and incomplete coverage.

## Planned evolution

1. schema and validation contract;
2. deterministic repository metadata/index builder;
3. structural and semantic parser backends;
4. references/call/dependency extraction;
5. test and V&V inventory;
6. semantic test/validation linkage;
7. incremental updates and agent-facing evidence queries;
8. reproducible AI tooling workflow and MCP query interface.

## Semantic test/V&V linkage

The linkage layer is deliberately evidence-driven. It creates links only from explicit repository references, such as a literal test source path in CMake registration or a validation manifest. Filename similarity, directory proximity, symbol-name similarity, and inferred execution coverage are not evidence and must not create links.

Each link records provenance in the `relation` column of `test_links` or `validation_links`. The query tool exposes these relations without upgrading them into execution or coverage claims. An absent link means no explicit evidence was found, not that the test is unrelated.
