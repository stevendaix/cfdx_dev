# CFDX Documentation

CFDX documentation is a single scientific publication system with four domains:

- **Theory** — executable CFD and numerical-method course.
- **User** — task-oriented usage and workflows.
- **Developer** — architecture, contracts and implementation.
- **V&V** — evidence and qualification governance.

## Source model

The authoritative source for quantitative Theory is version-controlled Python in Jupytext percent format. A Python source contains Markdown cells, equations, executable experiments and deterministic checks. Jupytext's percent format uses explicit `# %%` cell markers and is designed to remain an ordinary diffable Python file.

```text
Python .py
   │
   ├── scientific explanation
   ├── equations
   ├── executable experiment
   └── verification checks
   │
   ▼
MyST-NB + Sphinx
   │
   ├── HTML / PyData theme
   └── Typst publication pipeline
```

Committed `.ipynb` files are not the source of truth.

## Domains

```{toctree}
:maxdepth: 2

theory/README
user/README
developer/README
vv/README
```

## Theory course

Each chapter has a navigation README and an executable `chapter.py` source.

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
theory/16-data-model-file-formats/README
theory/17-computational-chain/README
theory/18-source-tree-physics-audit/README
```

## Executable Theory sources

```{toctree}
:hidden:

theory/00-foundations/chapter
theory/01-conservation-laws/chapter
theory/02-finite-volume-method/chapter
theory/03-meshes/chapter
theory/04-gradients-reconstruction/chapter
theory/05-fluxes/chapter
theory/06-time-integration/chapter
theory/07-pressure-velocity-coupling/chapter
theory/08-linear-algebra/chapter
theory/09-amg-mgr-schur/chapter
theory/10-turbulence/chapter
theory/11-heat-transfer/chapter
theory/12-radiation/chapter
theory/13-multiphysics/chapter
theory/14-numerical-analysis/chapter
theory/15-verification-validation/chapter
theory/16-data-model-file-formats/chapter
theory/17-computational-chain/chapter
theory/18-source-tree-physics-audit/chapter
```

## Existing documentation

Existing Markdown under `docs/application`, `docs/development`, `docs/validation`, `docs/user`, `docs/developer` and `docs/vv` remains available during migration. It must not become a second source of truth for Theory claims.

## Publication rule

Generated HTML, figures, tables and reports are outputs. Source equations, source data and executable experiments remain version-controlled inputs. V&V status comes from evidence, not manually copied prose.
