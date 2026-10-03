# CFDX Code Intelligence Index

The code-intelligence layer provides reproducible structural queries for CFDX agents.

## Purpose

The index complements, but never replaces, the repository source of truth. It may answer:

- where a symbol is defined;
- which files/symbols reference another symbol;
- callers and callees;
- include/import and dependency relationships;
- inheritance relationships;
- tests associated with implementation;
- validation cases exercising implementation;
- impact candidates for a proposed change.

Text search remains authoritative for literal occurrences. Structural index results are derived evidence and must be checked against current source before implementation or review decisions.

## Canonical logical model

The initial logical model contains:

- files — repository path, language and content identity;
- symbols — functions, methods, classes, structs, namespaces, modules and other definitions;
- references — symbol-to-symbol or symbol-to-file references;
- dependencies — include/import/build relationships;
- tests — test source and the implementation symbols they exercise when known;
- validation_cases — documented or executable V&V cases and the implementation symbols they exercise when known.

The model is intentionally vendor-neutral. Parser backends may use Clang AST/libclang, Python AST, Tree-sitter, LSP information, compile_commands.json, or other repository-appropriate mechanisms.

## Required properties

An implementation must be:

1. Reproducible.
2. Incrementally updateable.
3. Traceable to repository paths and source locations.
4. Version-aware.
5. Non-authoritative.
6. Deterministic.
7. Explicit about unsupported syntax and incomplete coverage.

## Query/evidence contract

Agents should treat an index answer as:

index result -> source confirmation -> relevant tests/docs -> decision

An index query alone is not sufficient evidence for an API change, numerical claim, review conclusion or qualification status.

## Initial schema

schema.sql defines the backend-neutral SQLite representation used by future indexers. SQLite is an implementation target, not part of the CFDX runtime ABI.

## Planned evolution

1. schema and validation contract;
2. deterministic repository metadata/index builder;
3. C++ structural parser backend;
4. Python structural parser backend;
5. references/call/dependency extraction;
6. test and V&V linkage;
7. incremental updates and query tooling.

The first implementation should remain dependency-light and should not add a runtime dependency to CFDX.
