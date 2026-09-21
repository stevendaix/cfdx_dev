#pragma once

#include <cstddef>

namespace cfdx::runtime::gpu {

// Backend-neutral kernels. The CPU reference path is deterministic and is
// also used by tests when CUDA is unavailable.
void gradient_gauss_reference(
    const double* phi,
    const double* face_sx,
    const double* face_sy,
    const double* face_sz,
    const std::size_t* owner,
    const std::int64_t* neighbour,
    const double* cell_volume,
    std::size_t n_faces,
    std::size_t n_cells,
    double* grad_x,
    double* grad_y,
    double* grad_z);

void divergence_reference(
    const double* phi_face,
    const double* face_sx,
    const double* face_sy,
    const double* face_sz,
    const std::size_t* owner,
    const std::int64_t* neighbour,
    std::size_t n_faces,
    std::size_t n_cells,
    double* div);

} // namespace cfdx::runtime::gpu
