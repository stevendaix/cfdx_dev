---
name: cfdx-cpp-development
description: Develop and audit CFDX C++20 code using project architecture, ownership, error handling and test conventions.
version: 1
---

# Required context

Inspect nearby CMake targets, headers, implementation files, interfaces, tests and existing style before editing.

# Procedure

1. Identify the owning target/module and public API boundary.
2. Trace lifetime, ownership, error propagation and thread/parallel assumptions.
3. Preserve existing abstractions unless the requirement explicitly changes them.
4. Implement the smallest coherent C++20 change.
5. Build the affected target and run relevant tests.
6. Review warnings, diagnostics and ABI/API impact.

# Numerical rule

For numerical code, pair implementation evidence with the applicable mathematical and verification workflow.

# Failure handling

Never hide compiler errors, warnings that are project-gated, sanitizer findings or test failures by weakening the gate.
