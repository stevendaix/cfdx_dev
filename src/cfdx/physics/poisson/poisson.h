#pragma once

#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

/// Canonical scalar Poisson problem:
///
///     -div(Gamma * grad(phi)) = S
///
/// where Gamma is a positive diffusivity and S is a volumetric source.
/// Cell-integrated source contributions therefore use S_c * V_c.
///
/// Boundary fluxes use the outward physical flux convention
///
///     q_n = Gamma * grad(phi) . n
///
/// so the integrated compatibility condition for a pure-Neumann problem is
///     integral(S dV) + integral(q_n dA) = 0.
enum class PoissonBoundaryType {
    DIRICHLET,
    NEUMANN,
};

/// Boundary data for the canonical Poisson contract.
///
/// Values are face-based and must have one entry per mesh face. A NaN entry
/// means that the face is not prescribed by this boundary object.
///
/// For DIRICHLET, face_values contains phi_b.
/// For NEUMANN, face_values contains q_n = Gamma * grad(phi) . n.
///
/// The explicit type prevents the old convention where an unspecified
/// Dirichlet value silently became a zero-gradient boundary.
struct PoissonBoundaryCondition {
    PoissonBoundaryType type = PoissonBoundaryType::DIRICHLET;
    std::vector<double> face_values;

    static PoissonBoundaryCondition dirichlet(std::size_t n_faces)
    {
        return {
            PoissonBoundaryType::DIRICHLET,
            std::vector<double>(n_faces, std::numeric_limits<double>::quiet_NaN())
        };
    }

    static PoissonBoundaryCondition neumann(std::size_t n_faces)
    {
        return {
            PoissonBoundaryType::NEUMANN,
            std::vector<double>(n_faces, std::numeric_limits<double>::quiet_NaN())
        };
    }

    void validate(std::size_t n_faces) const
    {
        if (face_values.size() != n_faces)
            throw std::invalid_argument(
                "PoissonBoundaryCondition: face_values size must equal mesh.n_faces()");
    }
};

/// Canonical input contract for scalar Poisson/Laplace solves.
///
/// A zero source represents the Laplace equation. The source is a volumetric
/// density and is integrated by the assembly as source[c] * cell_volume[c].
struct PoissonProblem {
    double diffusivity = 1.0;
    std::vector<double> source;
    PoissonBoundaryCondition boundary;

    void validate(std::size_t n_cells, std::size_t n_faces) const
    {
        if (!(diffusivity > 0.0))
            throw std::invalid_argument(
                "PoissonProblem: diffusivity must be positive");
        if (source.size() != n_cells)
            throw std::invalid_argument(
                "PoissonProblem: source size must equal mesh.n_cells()");
        boundary.validate(n_faces);
    }
};

} // namespace cfdx::core
