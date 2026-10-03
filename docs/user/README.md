# CFDX User Guide

**Status: IN PROGRESS — core getting-started, concepts, workflow and reference material are now documented; detailed task-specific reference remains under migration.**

The User domain answers how to use CFDX without requiring internal C++ knowledge. It describes supported behavior and points to executable evidence rather than promoting planned or unqualified capabilities.

```{toctree}
:maxdepth: 2

getting-started/README
getting-started/installation
getting-started/quickstart
concepts/README
workflow/README
reference/README
```

## Documentation boundary

- Theory explains the mathematics.
- User documentation explains how to operate the software.
- Developer documentation explains implementation contracts.
- V&V documentation defines evidence and acceptance.
- `docs/validation/` stores executed evidence.

Existing root documents such as `docs/boundary_conditions.md`, `docs/mesh_import.md` and `docs/gui.md` remain migration sources until reconciled.

User documentation must never silently present planned or unqualified capabilities as production features.
