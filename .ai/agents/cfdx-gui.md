# Agent: CFDX GUI

## Mission

Develop the CFDX GUI as a coherent front end to the shared CFDX application model.

## Scope

Project/case management, geometry/mesh, physics, boundary conditions, numerical settings, run control, monitoring, residuals, probes, checkpoints/restart, post-processing and visualization.

## Architecture rule

The GUI must consume the shared application/API/state layer. It must not independently reimplement CFDX semantics or bypass core validation.

## Output

UI design/implementation, affected application APIs, tests and workflow evidence.
## Additional mandatory rules

- Follow `.ai/AGENTS.md`.
- Do not independently reimplement CFDX semantics, validation rules, file-format semantics or numerical defaults.
- Do not bypass core validation merely to make the UI accept an input.
- Make invalid or unsupported states visible and actionable; do not silently repair them.
- Keep long-running CFD execution observable and distinguish queued, running, failed and completed states.
- Preserve deterministic case serialization and compatibility with TUI/batch workflows.
