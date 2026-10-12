---
name: cfdx-repository-navigation
description: Navigate the CFDX repository from an issue or change request to authoritative implementation, tests and documentation.
version: 1
---

# Purpose

Establish the repository context required before changing or auditing CFDX.

# Procedure

1. Read the complete issue/PR requirement and acceptance criteria.
2. Identify relevant top-level modules and documentation.
3. Locate candidate implementation symbols/files.
4. Locate tests, validation cases and CI workflows exercising the area.
5. Trace public APIs and configuration/data flow.
6. Record authoritative paths and unresolved ambiguities.

# Evidence rules

Current repository contents and executable tests take precedence over stale issue text or model assumptions.

# Failure handling

If a path, symbol or test cannot be found, report the absence explicitly and search adjacent authoritative locations before concluding it does not exist.
