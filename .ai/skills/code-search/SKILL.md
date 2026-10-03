---
name: cfdx-code-search
description: Perform layered textual and structural searches across CFDX source, tests, documentation and configuration.
version: 1
---

# Purpose

Find relevant code and evidence efficiently without treating search results as proof of behavior.

# Procedure

1. Start with exact symbols, identifiers and issue terminology.
2. Search implementation, tests, documentation and workflows separately.
3. Follow definitions, callers, derived classes and configuration paths.
4. Cross-check search results against build/test evidence.
5. Prefer structural code intelligence when available; otherwise use verified text-search tools.

# Rules

Distinguish where text occurs from what code actually executes. Never infer an implementation solely from a matching string.
