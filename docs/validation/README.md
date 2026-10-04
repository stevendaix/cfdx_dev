# CFDX Validation Evidence

**Status: executable evidence architecture — scientific campaign content is Python/Jupytext first.**

This tree stores **executed CFDX evidence**. The methodology and acceptance rules live in `docs/vv/`; validation does not redefine them.

## Source-of-truth model

- **`chapter.py` (Jupytext percent format)** — scientific campaign source: setup contract, equations, executable checks, data reduction, plots and quantitative gates.
- **`README.md`** — navigation and short scope only; it must not become a second source of quantitative truth.
- **`*.json`** — machine-readable registries and retained result records.
- **CI artifacts / raw logs** — immutable execution evidence for a campaign revision.

A Markdown document may remain when its role is genuinely governance, provenance, migration history or navigation. Scientific campaign results should not be maintained as hand-edited Markdown tables.

## Evidence lifecycle

```text
requirement → campaign contract → executable Jupytext source → CFDX execution
           → raw artifact → derived metrics → machine-readable result
           → review → promotion in the qualification registry
```

A reference value is never a CFDX result. A green CI job is not automatically validation or qualification.

## Campaign domains

```{toctree}
:maxdepth: 2

00-code-verification/README
01-mms/README
02-convergence/README
03-conservation/README
04-boundedness/README
05-linear-solvers/README
06-pressure-velocity/README
07-operators/README
08-turbulence/README
09-thermal/README
10-radiation/README
11-multiphysics/README
14-benchmarks/README
14-benchmarks/01-couette/chapter
14-benchmarks/02-poiseuille/chapter
14-benchmarks/03-ghia/chapter
14-benchmarks/04-vmfl036/chapter
14-benchmarks/05-bfs/chapter
14-benchmarks/06-naca0012/chapter
14-benchmarks/07-thermal/chapter
14-benchmarks/08-radiation/chapter
15-qualification/README
```

The numbered directories are the campaign taxonomy. A directory receives a `chapter.py` when executable evidence is actually implemented; we do not create empty Python notebooks merely to satisfy the layout.

## Ownership

- `docs/vv/` owns V&V methodology, terminology and acceptance-gate definitions.
- `docs/validation/` owns executed campaign evidence and reproducible data reduction.
- `14-benchmarks/` owns physical/reference benchmark campaigns.
- `15-qualification/` owns qualification evidence and promotion logic; the machine-readable registry remains authoritative for status.
- Historical plans and specialist contracts remain until their facts are absorbed and links are audited. They are not alternative status authorities.

## No-false-positive rule

Never change a tolerance, remove a failed case, relabel a diagnostic as qualification, or copy a reference result into a CFDX result merely to obtain a green report.

Promotion requires identifiable evidence for the declared scope, revision, configuration, metric and acceptance criterion.
