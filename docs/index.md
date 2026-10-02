# CFDX Documentation

CFDX documentation is organised as four complementary domains:

- **Theory** — CFD and numerical methods as a coherent scientific course.
- **User** — task-oriented workflows and supported usage.
- **Developer** — architecture, contracts and implementation.
- **V&V** — evidence from requirements through qualification.

## Documentation domains

```{toctree}
:maxdepth: 2
:caption: Documentation domains

theory/README
user/README
developer/README
vv/README
```

## Theory course

```{toctree}
:maxdepth: 2
:caption: Numerical methods

theory/00-foundations/README
theory/01-conservation-laws/README
theory/02-finite-volume-method/README
theory/03-meshes/README
theory/04-gradients-reconstruction/README
theory/05-fluxes/README
theory/06-time-integration/README
theory/07-pressure-velocity-coupling/README
theory/08-linear-algebra/README
theory/09-amg-mgr-schur/README
theory/10-turbulence/README
theory/11-heat-transfer/README
theory/12-radiation/README
theory/13-multiphysics/README
theory/14-numerical-analysis/README
theory/15-verification-validation/README
theory/16-cfdx-data-model-and-file-formats/README
theory/17-computational-chain-and-code-map/README
theory/18-source-tree-physics-audit/README
```

## Migration control

The migration map is the authoritative index for reconciling the existing documentation with the new publication domains.

```{toctree}
:maxdepth: 2
developer/migration/README
```

It records ownership, migration status, acceptance criteria and the next source families to migrate. It does not replace the source documents.

## Complete document hierarchy

The visible navigation is curated. A hidden global tree keeps every source document reachable during migration and prevents accidental orphaning.

```{toctree}
:hidden:
:glob:
README
application/*
development/*
theory/*/README
theory/04-gradients-reconstruction/*
user/*/README
developer/*/README
vv/*/README
validation/*
boundary_conditions
gui
mesh_import
```

## Complete document hierarchy

The visible navigation above is curated. A hidden global tree keeps every source document reachable during migration and prevents accidental orphaning. This follows the standard Sphinx toctree hierarchy mechanism.

```{toctree}
:hidden:
:glob:

README
application/*
development/*
theory/*/README
theory/04-gradients-reconstruction/*
user/*/README
developer/*/README
vv/*/README
validation/*
boundary_conditions
gui
mesh_import
```

## Existing documentation under migration

The existing material remains published until its replacement is proven. Migration is performed family-by-family; evidence is not deleted merely to make the new hierarchy cleaner.

```{toctree}
:maxdepth: 2
:caption: Existing documentation — migration source
README
application/*
development/*
boundary_conditions
gui
mesh_import
validation/*
```

These sources must eventually point to their new authoritative owner. Obsolete duplicates are removed only after the replacement is complete and covered by CI.

## Scientific publication

The documentation stack is Sphinx + MyST + PyData Sphinx Theme with MathJax and BibTeX. Executable scientific content is generated from version-controlled Python sources; generated figures and reports are outputs, not sources of truth.

## References

The project uses one authoritative BibTeX bibliography.

```{bibliography}
:all:
```
