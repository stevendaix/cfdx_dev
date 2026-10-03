# Workflow: Investigate a CFDX CI Failure

1. Identify the exact workflow, job, commit and failure step.
2. Determine whether the failure is infrastructure, build, test, validation or product behavior.
3. Reproduce locally when practical.
4. Inspect recent changes and affected dependencies.
5. Fix the underlying cause or document an external/infrastructure limitation.
6. Re-run the narrowest relevant check, then the required gate.
7. Do not hide failures by disabling tests, inflating tolerances or weakening workflow gates.
8. Record evidence and remaining uncertainty.