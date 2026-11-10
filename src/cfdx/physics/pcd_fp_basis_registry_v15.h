#pragma once

// V15 exact-CSR certificate selection. Every nonlinear outer iteration that
// rebuilds Fp selects a basis for THAT matrix: no stale basis cache and no
// silent fallback. The directory is prepared by an independent SVD/sparse
// certifier; on-the-fly SVD inside the CFD solve would be unbounded and is
// deliberately not performed. A missing new fingerprint is a hard failure.
#include "cfdx/physics/pcd_fp_bisided_basis_v13.h"
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

inline std::string pcd_fp_basis_registry_key_v15(
    const cfdx::core::SparseMatrix& matrix) {
    std::ostringstream out;
    out << std::hex << std::setfill('0') << std::setw(16)
        << pcd_fp_v13_hash_matrix(matrix);
    return out.str();
}

inline PcdFpBiNullspaceV13 load_pcd_fp_basis_registry_v15(
    const std::string& directory, const cfdx::core::SparseMatrix& matrix) {
    if (directory.empty() || !matrix.is_consistent() ||
        matrix.n_rows() == 0 || matrix.n_rows() != matrix.n_cols())
        throw std::invalid_argument("V15 Fp registry requires a valid square CSR and directory");
    const auto basis_file = std::filesystem::path(directory) /
        (pcd_fp_basis_registry_key_v15(matrix) + ".cfpb");
    if (!std::filesystem::is_regular_file(basis_file))
        throw std::runtime_error(
            "V15 Fp current numeric matrix has no certified basis: " +
            basis_file.string());
    return load_pcd_fp_bisided_basis_v13(basis_file.string(), matrix);
}

} // namespace cfdx::physics
