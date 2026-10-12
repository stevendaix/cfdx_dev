---
name: cfdx-python-development
description: Develop CFDX Python tooling, bindings and utilities consistently with the repository's API and test architecture.
version: 1
---

# Procedure

1. Inspect package layout, Python version requirements and existing public APIs.
2. Trace C++/pybind11 boundaries when bindings are involved.
3. Preserve pathlib-oriented path handling and explicit errors.
4. Add focused tests for changed behavior.
5. Run the narrowest relevant test suite, then required CI checks.

# Rules

Do not invent APIs from examples alone. Verify exported symbols and runtime behavior in the repository.
