# CFDX AI Permission Profiles

## READ_ONLY

May search and inspect repository code, documentation, issues, PRs, CI, cases and results. No mutation or uncontrolled execution.

## DEVELOPMENT

Adds local file modification, build/configure, unit tests and bounded validation workloads.

## MAINTAINER

Adds controlled branch/commit/PR/issue operations. Potentially expensive or destructive actions remain policy-controlled.

## Rules

- Least privilege by default.
- Destructive operations require explicit policy/confirmation.
- Expensive CFD workloads must be bounded or explicitly approved.
- Tool failures must be surfaced to the agent.
- Permissions for development and runtime MCP operations remain independently configurable.