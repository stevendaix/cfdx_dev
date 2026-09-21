#include "gpu_kernels.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace cfdx::runtime::gpu {

void gradient_gauss_reference(
    const double* phi,
    const double* sx, const double* sy, const double* sz,
    const std::size_t* owner, const std::int64_t* neighbour,
    const double* volume, std::size_t n_faces, std::size_t n_cells,
    double* gx, double* gy, double* gz)
{
    std::fill(gx, gx + n_cells, 0.0);
    std::fill(gy, gy + n_cells, 0.0);
    std::fill(gz, gz + n_cells, 0.0);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t o = owner[f];
        const std::int64_t n = neighbour[f];
        const double vf = (n >= 0) ? 0.5 * (phi[o] + phi[static_cast<std::size_t>(n)]) : phi[o];
        gx[o] += vf * sx[f]; gy[o] += vf * sy[f]; gz[o] += vf * sz[f];
        if (n >= 0) {
            const std::size_t j = static_cast<std::size_t>(n);
            gx[j] -= vf * sx[f]; gy[j] -= vf * sy[f]; gz[j] -= vf * sz[f];
        }
    }
    for (std::size_t c = 0; c < n_cells; ++c) {
        const double inv = volume[c] > 0.0 ? 1.0 / volume[c] : 0.0;
        gx[c] *= inv; gy[c] *= inv; gz[c] *= inv;
    }
}

void divergence_reference(
    const double* phi_face,
    const double* sx, const double* sy, const double* sz,
    const std::size_t* owner, const std::int64_t* neighbour,
    std::size_t n_faces, std::size_t n_cells, double* div)
{
    std::fill(div, div + n_cells, 0.0);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t o = owner[f];
        div[o] += phi_face[f];
        if (neighbour[f] >= 0)
            div[static_cast<std::size_t>(neighbour[f])] -= phi_face[f];
    }
    for (std::size_t c = 0; c < n_cells; ++c) div[c] /= 1.0; // caller supplies integrated flux
}

} // namespace cfdx::runtime::gpu
