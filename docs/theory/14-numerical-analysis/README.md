# Numerical Analysis

This chapter is a scientific course, not an API catalogue. It follows the chain: physical motivation → variables and units → assumptions → governing equations → integral formulation → discrete formulation → numerical properties → CFDX implementation → executable experiment → V&V evidence → limitations → improvements → references.

## Detailed course structure

### 1. Consistency

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 2. Truncation error

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 3. Stability

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 4. Convergence

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 5. Accuracy

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 6. Dissipation and dispersion

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 7. Monotonicity

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 8. Boundedness

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 9. Conditioning

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 10. Error budget

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 11. Observed order

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

### 12. Mesh-quality effects

- Physical motivation and scope.
- Definitions, symbols, dimensions and SI units.
- Assumptions and domain of validity.
- Continuous mathematical formulation.
- Control-volume and finite-volume formulation where applicable.
- Discrete/algebraic formulation and sign conventions.
- Conservation, consistency, stability, boundedness and accuracy properties.
- Exact CFDX implementation path and source files.
- Executable verification experiment or test.
- Benchmark/reference evidence where applicable.
- Limitations, failure modes and improvement paths.

## Mandatory equation template

Every important equation must state its physical origin, assumptions, variables and units, coordinate/sign convention, continuous form, integral form, discrete approximation, algebraic contribution, numerical properties, implementation path, executable verification, benchmark/reference, limitations and bibliography.

## Evidence vocabulary

**Implemented** means a code path exists. **Verified** means a defined mathematical property has executable evidence. **Validated** means the computed physical result has been compared against an appropriate independent reference/experiment. **Qualified** means the declared V&V population and acceptance gates support the intended scope. These states must never be conflated.

## Scientific figures

Use generated figures for geometry, control volumes, stencils, matrix/block structure, convergence and error studies. Figures are outputs of executable sources and are never the source of truth. Interactive Plotly/Altair/PyVista material should have a static interpretation where practical.

## Repository traceability

Every implementation claim must point to the actual CFDX source and test evidence discovered during audit. Missing or partial functionality must be marked explicitly; no undocumented API or numerical result may be invented.
