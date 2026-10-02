# CFDX Verification, Validation & Qualification

The V&V documentation defines how CFDX correctness and maturity are demonstrated.

The detailed evidence store currently lives under [docs/validation/](../validation/README.md). This `docs/vv/` tree is the stable publication/navigation layer. Migration of evidence documents is incremental and governed by the V&V governance rules.

## V&V chain

```text
requirement
    ↓
implementation
    ↓
code verification
    ↓
solution verification
    ↓
validation (when physically applicable)
    ↓
qualification
```

## Authoritative current material

- [V&V governance](../validation/CFDX_VV_GOVERNANCE.md)
- [Qualification campaign](../validation/CFDX_QUALIFICATION_CAMPAIGN.md)
- [Qualification case template](../validation/CFDX_QUALIFICATION_CASE_TEMPLATE.md)
- [#461 requirements audit](../validation/NUMERICAL_METHOD_REQUIREMENTS_AUDIT.json)
- [Numerical maturity audit](../validation/NUMERICAL_METHOD_MATURITY_AUDIT.md)
- [Validation reporting](../validation/VALIDATION_REPORT.md)
- [Documentation map](../validation/DOCUMENTATION_MAP.md)

## Rules

- Evidence is executable whenever practical.
- Acceptance criteria are declared before observing the result.
- Failures are retained.
- Human-readable documents summarize evidence; machine-readable artifacts remain the regression source of truth.
- Status promotion requires the corresponding evidence.
- Cross-code agreement is supporting evidence, not an automatic oracle.

See [CFDX V&V governance](../validation/CFDX_VV_GOVERNANCE.md) for the complete rules.
