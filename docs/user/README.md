# CFDX User Guide

**Status: TOPO — migration of existing user documentation is not complete.**

The User domain answers **how to use CFDX** without requiring knowledge of internal C++ implementation.

## Planned navigation

### Getting started
- installation and prerequisites;
- build and first execution;
- first case;
- mesh creation/import;
- selecting physics and numerics;
- reading results.

### Concepts
- cases and files;
- meshes and regions;
- fields and boundary conditions;
- physics models;
- numerics and solver selection;
- convergence and stopping criteria;
- checkpoints and restart;
- post-processing.

### Workflow
- case setup;
- mesh workflow;
- physics setup;
- numerics setup;
- execution and monitoring;
- restart;
- probes/time series;
- result export;
- V&V workflow.

### Reference
- case schema;
- numerics options;
- boundary conditions;
- solver options;
- output formats;
- command-line/TUI/GUI interfaces.

## Existing material

The repository already contains useful documents such as `docs/boundary_conditions.md`, `docs/mesh_import.md` and `docs/gui.md`. These are migration sources and must be reconciled with the User domain before removal or deprecation.

## Rule

User documentation describes supported behavior. It must not silently document planned or unqualified capabilities as production features.
