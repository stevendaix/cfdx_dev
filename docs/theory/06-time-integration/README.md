# Time Integration

This chapter is a scientific course, not an API catalogue. It follows the chain: physical motivation → variables and units → assumptions → governing equations → integral formulation → discrete formulation → numerical properties → CFDX implementation → executable experiment → V&V evidence → limitations → improvements → references.

## Detailed course structure

### 1. Semi-discrete equations

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

### 2. Explicit Euler

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

### 3. Implicit Euler

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

### 4. Crank-Nicolson

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

### 5. BDF2

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

### 6. Variable time step

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

### 7. CFL and stability

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

### 8. Pseudo-transient

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

### 9. Time-step control

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

### 10. Restart history

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

### 11. Temporal verification

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
