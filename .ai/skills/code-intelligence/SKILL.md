# Skill: Code Intelligence

## Purpose

Provide agents with structural repository navigation that is more precise than literal text search while preserving source-grounded decisions.

## Current executable coverage

- Repository inventory is deterministic and hash-based.
- Python files are structurally parsed with the standard-library AST for classes, functions and import dependencies.
- C/C++ have two explicit paths:
  - semantic backend: libclang driven by CMake `compile_commands.json`, when optional Python bindings and compatible libclang are available;
  - fallback backend: conservative declaration/include scanning for coarse navigation.
- Parser provenance and coverage are recorded in `index_metadata`; parse failures and Clang diagnostics are retained.
- Semantic results are limited to successfully parsed translation units and their actual compilation arguments. The index does not establish whole-program completeness, macro semantics, or virtual-dispatch completeness.
- When libclang is unavailable, the tool exits explicitly and the heuristic backend remains the only C/C++ structural coverage. Heuristic results must never be presented as semantic claims.

## Operating procedure

1. Establish the repository revision and relevant build configuration.
2. Prefer the CMake-generated `compile_commands.json` for C++.
3. Use exact per-translation-unit compiler arguments with output-only flags removed.
4. Prefer libclang semantic parsing when available.
5. Record parser provenance, coverage and diagnostics.
6. Confirm important results against current source, tests and build metadata.
7. Expand to textual search when index coverage is incomplete.

## Evidence rules

- An index result is derived evidence, not source truth.
- A missing index result does not prove that a symbol or relationship is absent.
- A stale index must be regenerated or bypassed.
- Semantic references require a successful Clang parse of the relevant translation unit.
- Numerical, V&V and qualification conclusions require their corresponding evidence; code intelligence alone cannot establish them.

## Query patterns

### Find a definition

1. query indexed symbols by qualified name/name;
2. confirm declaration/definition in source;
3. inspect overloads/templates and owning type/namespace.

### Find callers

1. locate the target symbol;
2. query references targeting it;
3. inspect call sites and dispatch/virtual semantics;
4. confirm relevant build/test coverage.

### Impact analysis

Inspect direct references, dependent files/modules, tests, validation cases, public documentation, and CMake/build ownership. Do not assume complete index coverage.

## Update policy

Index generation should remain deterministic and incremental. Store parser provenance, schema/index version, repository revision and content hashes. Never commit disposable index databases.

## Completion criteria

Code intelligence is complete only when structural queries are reproducible, stale/partial coverage is observable, semantic C++ parsing is exercised where the optional backend is available, fallback behavior is explicit, tests cover compile-command handling and semantic relationships, and no runtime CFDX dependency is introduced accidentally.
