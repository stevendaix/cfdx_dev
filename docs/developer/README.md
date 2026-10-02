# CFDX Developer Documentation

Developer documentation describes how CFDX is built, structured, tested and extended.

It is the contract between the scientific design and the software implementation.

## Structure

```text
architecture/
development/
implementation/
api/
contributing/
```

## Theory boundary

Developer pages may state mathematical contracts, but the full derivation belongs to Theory.

A developer page should answer:

- what interface exists;
- what invariants it guarantees;
- what data it consumes and produces;
- what ownership/lifetime rules apply;
- what errors are reported;
- how it is tested.

It should link to Theory for the numerical method and to V&V for the evidence.

## Implementation contract

A numerical implementation document should identify the mathematical method, software component, inputs/outputs, invariants, supported and unsupported configurations, verification hooks, relevant tests and related issue/PR.

Implementation presence is never equivalent to qualification.
