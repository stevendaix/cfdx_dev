# Boundary condition schema

The boundary editor is patch-aware: the available form is selected from the real mesh patch type, not from a static GUI list. Scalar boundary fields mirror the existing C++ `ScalarBoundaryCondition` contract (`FIXED_VALUE`, `ZERO_GRADIENT`, `FIXED_GRADIENT`, with value and gradient). Unsupported patch types are rejected by the schema instead of silently mapped to an unrelated condition.
