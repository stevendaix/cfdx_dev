#pragma once

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace cfdx::core {

enum class PoissonBoundaryType {
    DIRICHLET,
    NEUMANN,
};

/// Boundary data for the canonical Poisson contract.
///
/// type remains the default type for all faces. face_types, when populated,
/// provides an explicit type per boundary face and enables mixed conditions.
///
/// For DIRICHLET, face_values contains phi_b.
/// For NEUMANN, face_values contains q_n = Gamma * grad(phi) . n.
///
/// A NaN value means that the corresponding boundary face is not prescribed.
struct PoissonBoundaryCondition {
    PoissonBoundaryType type = PoissonBoundaryType::DIRICHLET;
    std::vector<double> face_values;
    std::vector<PoissonBoundaryType> face_types;

    static PoissonBoundaryCondition dirichlet(std::size_t n_faces)
    {
        return {
            PoissonBoundaryType::DIRICHLET,
            std::vector<double>(
                n_faces, std::numeric_limits<double>::quiet_NaN()),
            {}
        };
    }

    static PoissonBoundaryCondition neumann(std::size_t n_faces)
    {
        return {
            PoissonBoundaryType::NEUMANN,
            std::vector<double>(
                n_faces, std::numeric_limits<double>::quiet_NaN()),
            {}
        };
    }

    static PoissonBoundaryCondition mixed(
        std::size_t n_faces,
        PoissonBoundaryType default_type = PoissonBoundaryType::DIRICHLET)
    {
        return {
            default_type,
            std::vector<double>(
                n_faces, std::numeric_limits<double>::quiet_NaN()),
            std::vector<PoissonBoundaryType>(n_faces, default_type)
        };
    }

    PoissonBoundaryType type_for_face(std::size_t face) const
    {
        return face_types.empty() ? type : face_types.at(face);
    }

    void validate(std::size_t n_faces) const
    {
        if (!face_values.empty() && face_values.size() != n_faces)
            throw std::invalid_argument(
                "PoissonBoundaryCondition: face_values size must equal "
                "mesh.n_faces() when provided");
        if (!face_types.empty() && face_types.size() != n_faces)
            throw std::invalid_argument(
                "PoissonBoundaryCondition: face_types size must equal "
                "mesh.n_faces() when provided");
        for (double value : face_values) {
            if (!std::isfinite(value) && !std::isnan(value))
                throw std::invalid_argument(
                    "PoissonBoundaryCondition: face values must be finite "
                    "or NaN (unprescribed)");
        }
    }
};

enum class PoissonGaugeType {
    REFERENCE_CELL,
};

struct PoissonGauge {
    PoissonGaugeType type = PoissonGaugeType::REFERENCE_CELL;
    std::size_t reference_cell = 0;
    double reference_value = 0.0;
};

struct PoissonProblem {
    double diffusivity = 1.0;
    std::vector<double> source;
    PoissonBoundaryCondition boundary;

    void validate(std::size_t n_cells, std::size_t n_faces) const
    {
        if (!(diffusivity > 0.0) || !std::isfinite(diffusivity))
            throw std::invalid_argument(
                "PoissonProblem: diffusivity must be finite and positive");
        if (source.size() != n_cells)
            throw std::invalid_argument(
                "PoissonProblem: source size must equal mesh.n_cells()");
        for (double value : source) {
            if (!std::isfinite(value))
                throw std::invalid_argument(
                    "PoissonProblem: source values must be finite");
        }
        boundary.validate(n_faces);
    }
};

} // namespace cfdx::core
