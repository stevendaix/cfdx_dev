#pragma once

#include "cfdx/physics/steady_incompressible_solver.h"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>

namespace cfdx::io {

// Persist the authoritative solver convergence result in the canonical runtime
// artifact. This is execution evidence only; it is not a validation verdict.
inline void write_convergence_history(
    const std::string& path,
    const cfdx::physics::IncompressibleSolveResult& result)
{
    if (!std::isfinite(result.reference_momentum_residual))
        throw std::invalid_argument("write_convergence_history: non-finite reference residual");

    const std::string tmp_path = path + ".tmp";
    std::ofstream out(tmp_path, std::ios::trunc);
    if (!out)
        throw std::runtime_error("write_convergence_history: cannot open " + path);
    out << std::setprecision(std::numeric_limits<double>::max_digits10);

    try {
        out << "{\n"
            << "  \"format\": \"CFDX-CONVERGENCE\",\n"
            << "  \"schema_version\": 1,\n"
            << "  \"converged\": " << (result.converged ? "true" : "false") << ",\n"
            << "  \"iterations\": " << result.iterations << ",\n"
            << "  \"reference_momentum_residual\": " << result.reference_momentum_residual << ",\n"
            << "  \"convergence_status\": " << static_cast<int>(result.convergence_status) << ",\n"
            << "  \"convergence_reason\": \"" << result.convergence_reason << "\",\n"
            << "  \"history_count\": " << result.history.size() << "\n"
            << "}\n";
    } catch (...) {
        out.close();
        std::remove(tmp_path.c_str());
        throw;
    }

    if (!out) {
        out.close();
        std::remove(tmp_path.c_str());
        throw std::runtime_error("write_convergence_history: write failed for " + path);
    }
    out.close();
    if (std::rename(tmp_path.c_str(), path.c_str()) != 0) {
        std::remove(tmp_path.c_str());
        throw std::runtime_error("write_convergence_history: cannot replace " + path);
    }
}

} // namespace cfdx::io
