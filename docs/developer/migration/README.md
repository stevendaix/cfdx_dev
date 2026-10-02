# Documentation Migration Map

**Status: ACTIVE MIGRATION — this page is the authoritative map for reconciling the legacy documentation tree with the four-domain publication model.**

The migration is deliberately incremental. Existing documents remain valid migration sources until their content has an explicit owner in **Theory**, **User**, **Developer**, or **V&V**.

## Ownership rules

| Content | Primary owner |
|---|---|
| Mathematical definition, derivation, numerical properties | Theory |
| Task-oriented procedure and supported workflow | User |
| Software architecture, contracts, implementation details | Developer |
| Requirements, verification, validation, qualification evidence | V&V |

A document can link to another domain, but a technical fact must have one authoritative owner.

## Current migration inventory

| Legacy/source family | Primary destination | Migration action | Status |
|---|---|---|---|
| `docs/application/` | User | Split task/workflow material from implementation details; preserve historical material until reconciled | IN PROGRESS |
| `docs/development/` | Developer | Reconcile architecture, build, testing, implementation and contribution material with the new Developer sections | IN PROGRESS |
| `docs/boundary_conditions.md` | Theory + Developer + User | Theory owns mathematical BC definitions; Developer owns implementation contracts; User owns configuration workflow | PLANNED |
| `docs/mesh_import.md` | User + Developer + Theory | User owns import workflow; Developer owns supported formats/contracts; Theory owns mesh concepts | PLANNED |
| `docs/gui.md` | User + Developer | User owns workflows; Developer owns GUI/core boundary and architecture | PLANNED |
| `docs/validation/` | V&V | Reconcile evidence documents into the stable V&V portal without deleting evidence | IN PROGRESS |
| `docs/theory/04-gradients-reconstruction/` | Theory | Keep as the executable scientific-documentation pilot and template | COMPLETE PILOT |
| `docs/theory/*/README.md` | Theory | Replace topology pages progressively with complete course chapters | TOPO + PILOT |
| `docs/user/*/README.md` | User | Replace topology pages with task-oriented guides based on existing workflows | TOPO |
| `docs/developer/*/README.md` | Developer | Replace topology pages with repository-grounded implementation documentation | TOPO |
| `docs/vv/*/README.md` | V&V | Build stable navigation over existing evidence; do not duplicate evidence unnecessarily | TOPO + MIGRATION |

## Migration procedure

For each source document:

1. Identify the factual claims and their current status.
2. Split mixed documents by ownership rather than copying the whole document into several domains.
3. Move mathematical explanation to **Theory**.
4. Move supported usage and procedures to **User**.
5. Move implementation contracts and architecture to **Developer**.
6. Move executable evidence and maturity claims to **V&V**.
7. Link the resulting pages instead of duplicating the same facts.
8. Keep historical material explicitly marked as historical when it is still useful.
9. Remove an obsolete source only after its replacement is complete, linked, and covered by CI.
10. Update the documentation map and navigation in the same change.

## Migration acceptance criteria

A migrated family is complete only when:

- every retained fact has a single authoritative owner;
- links from the old entry point lead to the new source of truth;
- implementation claims are backed by repository evidence;
- numerical capability claims point to V&V evidence;
- no document is orphaned from the Sphinx hierarchy;
- obsolete duplicates are removed only after replacement is proven;
- the migration does not silently change solver capability or qualification status.

## Immediate next families

The first detailed migration pass is:

1. boundary conditions;
2. mesh import;
3. GUI;
4. application/user workflows;
5. development/build/CI material;
6. validation evidence families.

The N2 gradients chapter remains the reference format for future executable Theory chapters.
