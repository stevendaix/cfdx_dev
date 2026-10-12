# 13. CFDX Implementation Contract

Theory defines the numerical method; implementation must preserve that definition.

The CFDX contract should identify:

- gradient API and data structures;
- supported dimensions and cell types;
- stencil construction;
- weighting policy;
- conditioning diagnostics;
- boundary policy;
- face reconstruction interface;
- limiter interface;
- diagnostic outputs;
- failure semantics;
- parallel and deterministic behaviour.

The authoritative relationship is:

```
Theory
  ↓
Numerical contract
  ↓
Implementation
  ↓
Executable tests
  ↓
V&V evidence
```

Implementation details belong in the Developer domain; this chapter states the mathematical consequences and required contract.
