# Skill: Code Intelligence

## Purpose

Provide agents with structural repository navigation that is more precise than literal text search while preserving source-grounded decisions.

## Scope

Use this skill for symbol definitions and references, callers/callees, include/import and dependency analysis, inheritance and ownership relationships, implementation impact analysis, and locating tests and validation cases connected to implementation.

This skill does not define a specific vendor, editor, LLM, parser or MCP implementation.

## Operating procedure

1. Establish the repository revision and relevant build configuration.
2. Start with the narrowest useful query.
3. Use the structural index when the question depends on relationships rather than literal text.
4. Confirm important results against current source, tests and build metadata.
5. Expand to textual search when index coverage is incomplete.
6. Preserve uncertainty when parsing is partial, stale or ambiguous.
7. Report source locations supporting the conclusion.

## Current executable coverage

- Repository inventory is deterministic and hash-based.
- Python files are structurally parsed with the standard-library AST module for classes, functions and import dependencies.
- C/C++ files currently use a deliberately conservative declaration/include scanner. This is useful for coarse navigation only and is **not** semantic C++ parsing.
- Parser provenance is stored in index_metadata; parse failures are recorded rather than silently treated as empty files.
- C++ semantic references, overload resolution, templates, inheritance and macro-aware analysis remain future work and must not be inferred from the heuristic scanner.

## Evidence rules

- An index result is derived evidence, not source truth.
- A missing index result does not prove that a symbol or relationship is absent.
- A stale index must be regenerated or bypassed.
- Generated identifiers must not be used as permanent API identifiers.
- Structural relationships should be rechecked after source changes.
- Numerical, V&V and qualification conclusions require their corresponding evidence; code intelligence alone cannot establish them.

## Parser strategy

Prefer the strongest available parser for the language and question:

- C++: compile database plus Clang AST/libclang where available;
- Python: Python AST;
- cross-language/text fallback: Tree-sitter or repository search;
- build/dependency context: CMake metadata and compile_commands.json;
- editor/LSP data may supplement, but must not become a hidden authoritative dependency.

Unsupported syntax and parser failures must be surfaced as incomplete coverage.

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

For an API or numerical implementation change, inspect direct references, dependent files/modules, tests, validation cases, public documentation, and CMake/build ownership.

Do not assume that a complete index implies complete impact coverage.

## Update policy

Index generation should be deterministic and incremental. Store repository revision, schema/index version and content hashes. Never commit disposable generated index databases unless a future repository policy explicitly requires it.

## Completion criteria

A code-intelligence implementation is complete only when structural queries are reproducible, stale/partial coverage is observable, C++ and Python paths are covered to their documented scope, tests exercise indexing and query semantics, failure modes are explicit, and no runtime CFDX dependency is introduced accidentally.
