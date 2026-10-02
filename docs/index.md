# CFDX Documentation

CFDX documentation is organised as four complementary domains:

- **Theory** — CFD and numerical methods as a coherent scientific course.
- **User** — task-oriented workflows and supported usage.
- **Developer** — architecture, contracts and implementation.
- **V&V** — evidence from requirements through qualification.

## Documentation domains

```{toctree}
:maxdepth: 2

theory/README
user/README
developer/README
vv/README
```

## Theory course

```{toctree}
:maxdepth: 2

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
```

## Complete document hierarchy

The visible navigation above is curated. A hidden global tree keeps every source document reachable during migration and prevents accidental orphaning. This follows the standard Sphinx toctree hierarchy mechanism. citeturn2search1

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

The repository already contains useful material outside the new domain structure. Until migration is complete, these documents remain published through a dedicated legacy/migration tree.

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

These are migration sources. Their content and status must be reconciled with the new Theory/User/Developer/V&V source-of-truth rules before obsolete copies are removed.

## Scientific publication

The documentation stack is Sphinx + MyST + PyData Sphinx Theme with MathJax and BibTeX. Executable scientific content is generated from version-controlled Python sources; generated figures and reports are outputs, not sources of truth.
