# 1. Introduction

Gradient reconstruction is a central numerical operation in finite-volume CFD. Cell-centred values are known from the discrete solution, while viscous terms, diffusion, higher-order convection, boundary reconstruction and several physical models require spatial derivatives.

This chapter develops the subject from the continuous gradient to practical polyhedral finite-volume algorithms. The emphasis is on mathematical assumptions, geometry, accuracy, conditioning, boundary treatment and reproducible verification.

## Learning objectives

After this chapter, the reader should be able to:

- define and interpret a scalar-field gradient;
- derive Green–Gauss and least-squares reconstruction;
- explain weighted least squares on irregular meshes;
- distinguish gradient reconstruction from face-value reconstruction;
- diagnose rank deficiency and poor conditioning;
- distinguish linear exactness from demonstrated asymptotic order;
- design reproducible gradient verification studies;
- connect theory, CFDX implementation contracts and V&V evidence.

## Roadmap

The progression is:

```
continuous field
      ↓
continuous gradient
      ↓
finite-volume discretisation
      ↓
geometric reconstruction
      ↓
gradient algorithm
      ↓
face reconstruction
      ↓
flux evaluation
      ↓
verification
```

A gradient method must therefore be assessed together with the geometric and reconstruction assumptions surrounding it.
