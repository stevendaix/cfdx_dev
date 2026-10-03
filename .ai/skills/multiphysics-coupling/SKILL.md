# Skill: Multiphysics Coupling

## Purpose
Provide a rigorous workflow for coupled thermal, radiation, species and other multiphysics extensions.

## Method
1. Define each field equation and coupling term independently.
2. Identify the discrete exchange terms and shared conservation quantities.
3. Check action/reaction consistency of inter-field and interface fluxes.
4. Separate monolithic, segregated and operator-split coupling.
5. Verify each physics component before claiming coupled validation.
6. Check dimensional consistency and coupling convergence independently.

## Required evidence
- Single-physics regression.
- Controlled coupling test with known solution or manufactured source.
- Interface/exchange conservation.
- Energy/species balance where applicable.
- Coupling iteration and time-step sensitivity.

## Review rules
Do not attribute coupled agreement to one component without isolating the coupling contribution. Do not silently fall back to decoupled physics.
