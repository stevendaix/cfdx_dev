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
```

## Existing documentation to migrate

The repository already contains useful material outside the new four-domain structure, including boundary-condition, mesh-import, GUI, development and validation documents. These are migration sources, not automatically authoritative duplicates.

Migration must preserve links and provenance before obsolete copies are removed.

## Scientific publication

The documentation stack is Sphinx + MyST + PyData Sphinx Theme with MathJax and BibTeX. Executable scientific content is generated from version-controlled Python sources; generated figures and reports are outputs, not sources of truth.
