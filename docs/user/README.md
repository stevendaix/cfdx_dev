# CFDX User Documentation

User documentation answers practical questions: how to create a case, configure a model, run a simulation, inspect results, restart a calculation and diagnose a failure.

It must remain independent of internal implementation details.

## Structure

```text
getting-started/
concepts/
workflow/
reference/
examples/
```

## Rules

- Start from the user's task, not from the C++ class hierarchy.
- Do not document an implementation detail unless it changes user-visible behaviour.
- Every supported option must have a defined status.
- Unsupported or experimental features must be labelled explicitly.
- Commands and configuration examples must be executable or marked as illustrative.
- User documentation must not create a second source of truth for V&V status.
