# CFDX AI Framework

This directory is the vendor-neutral source of truth for AI-assisted CFDX development and operation.

## Purpose

The framework defines repository-grounded AI knowledge, reusable skills, specialized agents, workflows, future MCP interfaces, permission rules and vendor adapters.

The framework is deliberately independent of a specific model or vendor.

## Source of truth

The contents of `.ai/` are canonical. Vendor-specific configurations must reference or derive from these files rather than becoming independent copies.

## Evidence hierarchy

1. Current executable repository evidence.
2. Current tests and validation results.
3. Current CFDX documentation.
4. Current issues, PRs and CI evidence.
5. Explicitly identified external references.
6. Model knowledge, only as a hypothesis to verify.

Agents must distinguish implementation, testing, verification, validation and qualification.

## Numerical integrity

AI-assisted development follows the same CFDX rules as human development: no tolerance inflation, no disabling tests merely to obtain green CI, no silent numerical fallback, no qualification claim from implementation alone, and preservation of negative evidence.

## MCP split

Two logical MCP families are planned: a Development MCP for code, symbols, documentation, tests, validation metadata, CI and GitHub; and a Runtime MCP for CFDX cases, meshes, execution, monitoring, checkpoints, fields, probes and validation.

Their permissions and responsibilities remain distinct.

## Status

This foundation phase is documentation/configuration only. Code intelligence, MCP servers and vendor integrations are subsequent phases tracked by the parent issue.

## Project context

The framework complements the existing CFDX roadmap and numerical maturity work. It does not redefine CFDX qualification levels.
## Global agent rules

`.ai/AGENTS.md` defines mandatory cross-cutting operating rules for every agent, workflow and MCP. It covers repository-grounded reasoning, scope discipline, engineering/numerical integrity, evidence classification, reproducibility, permissions, security and completion criteria. More restrictive skill- or tool-specific rules may add constraints but must not weaken these rules.
