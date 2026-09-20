#pragma once

#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include <vector>

namespace cfdx::core::numerics {

// ============================================
// Vertical Slice 1 (Poisson) — Laplacian Assembly
// ============================================

//! Assemble the Laplacian operator in CSR format for the Poisson equation.
//!
//! Given the mesh topology (faces, owner, neighbour) and geometry
//! (cell volumes, face areas, face normals), this constructs the sparse
//! matrix A for -∇² p = b (discretized via FVM Gauss divergence).
//!
//! The assembly follows:
//!   A_ii = Σ_f (|Sf| / |d_of|) * (owner+neighbour contribution)
//!   A_ij = - Σ_f (|Sf| / |d_of|)  for face f connecting owner i and neighbour j
//!
//! This produces a symmetric positive-definite matrix suitable for CG
//! with Hypre BoomerAMG preconditioner.
//!
//! Args:
//!   mesh: CFDX mesh with topology and geometry
//!   A: output SparseMatrix (CSR)
//! Returns: true on success
bool assembleLaplacianCSR(const Mesh& mesh, SparseMatrix& A);

//! Assemble the right-hand side vector b from a source term field.
//! Args:
//!   mesh: CFDX mesh
//!   source: Source term scalar field (e.g. divergence of a velocity field)
//!   b: output Vector
//! Returns: true on success
bool assemblePoissonRHS(const Mesh& mesh, const ScalarCellField& source, Vector& b);

//! Convenience wrapper: build full Poisson system A*p = b from source field.
//! Returns LinearSystem with assembled A and b.
LinearSystem buildPoissonSystem(const Mesh& mesh, const ScalarCellField& source);

} // namespace cfdx::core::numerics
