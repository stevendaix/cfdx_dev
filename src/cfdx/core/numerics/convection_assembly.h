#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include <vector>

namespace cfdx::core::numerics {

//! Assemble convection term: ∇·(U ⊗ u) for Navier-Stokes momentum equation.
//!
//! Uses FVM interpolation (upwind, linear, or limited schemes) to compute
//! face fluxes F_f = U_f · Sf and interpolate velocity to faces.
//!
//! Args:
//!   mesh: CFDX mesh
//!   velocity: velocity field (U, V, W components)
//!   scheme: "upwind", "linear", or "limited"
//!   component_idx: 0=U, 1=V, 2=W (component of velocity being assembled)
//!   A: output SparseMatrix (CSR)
//! Returns: true on success
bool assembleConvectionCSR(const Mesh& mesh,
                           const std::vector<ScalarCellField>& velocity,
                           const std::string& scheme,
                           int component_idx,
                           SparseMatrix& A);

//! Convenience wrapper for momentum equation assembly (convection + diffusion + pressure gradient).
bool assembleMomentumCSR(const Mesh& mesh,
                        const std::vector<ScalarCellField>& velocity,
                        const ScalarCellField& pressure,
                        const std::string& convection_scheme,
                        SparseMatrix& A,
                        Vector& b);

} // namespace cfdx::core::numerics
