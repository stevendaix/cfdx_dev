#pragma once

#include "cfdx/physics/finite_volume_transport.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {

// Positive diagonal approximation to the Jacobian of the fully reconstructed
// deferred TVD face flux, estimated from TWO accepted Picard assembly inputs.
// The extra diagonal J and RHS J*U_old cancel EXACTLY on the current state.
// Consequently the true accepted-state FV operator and its physical residual
// remain the original second-order TVD operator, not this stabilization.
// A negative secant derivative is never injected as an antidiffusive diagonal.
struct TvdSecantImplicitControls {
    bool enabled = false;
    double spectral_radius_cap = 0.5;
    double minimum_relative_displacement = 1e-7;
};
inline void validate_tvd_secant_implicit(const TvdSecantImplicitControls& c) {
    if (!std::isfinite(c.spectral_radius_cap) ||
        !(c.spectral_radius_cap > 0.0 && c.spectral_radius_cap <= 2.0) ||
        !std::isfinite(c.minimum_relative_displacement) ||
        !(c.minimum_relative_displacement > 0.0 &&
          c.minimum_relative_displacement < 1.0))
        throw std::invalid_argument("TVD secant implicit: invalid stabilization controls");
}
struct TvdSecantImplicitResult {
    double max_diagonal_added = 0.0;
    double max_ratio_to_face_flux = 0.0;
    std::size_t regularized_cells = 0;
};

inline std::vector<double> tvd_deferred_source_by_cell(
    const cfdx::core::Mesh& mesh, const std::vector<double>& face_correction)
{
    if (face_correction.size() != mesh.n_faces())
        throw std::invalid_argument("TVD secant: incorrect flux size");
    std::vector<double> rhs(mesh.n_cells(),0.0);
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const auto n=mesh.ownership().neighbour(f);
        if (n < 0) continue;
        const auto o=mesh.ownership().owner(f);
        if (o>=rhs.size() || static_cast<std::size_t>(n)>=rhs.size() ||
            !std::isfinite(face_correction[f]))
            throw std::invalid_argument("TVD secant: invalid face correction");
        rhs[o]-=face_correction[f];
        rhs[static_cast<std::size_t>(n)]+=face_correction[f];
    }
    return rhs;
}

inline TvdSecantImplicitResult add_tvd_secant_implicit_diagonal(
    const cfdx::core::Mesh& mesh,
    ScalarEquation& eq,
    const std::vector<double>& raw_current,
    const std::vector<double>& raw_previous,
    const std::vector<double>& u_current,
    const std::vector<double>& u_previous,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& mass_flux,
    const TvdSecantImplicitControls& controls)
{
    validate_tvd_secant_implicit(controls);
    const auto nc=mesh.n_cells();
    if (u_current.size()!=nc || u_previous.size()!=nc ||
        eq.diagonal.size()!=nc || eq.rhs.size()!=nc ||
        eq.matrix.n_rows()!=nc || eq.matrix.n_cols()!=nc ||
        mass_flux.size()!=mesh.n_faces())
        throw std::invalid_argument("TVD secant: invalid input dimensions");
    const auto b=tvd_deferred_source_by_cell(mesh,raw_current);
    const auto previous=tvd_deferred_source_by_cell(mesh,raw_previous);
    std::vector<double> spectral(nc,0.0);
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const double flux=mass_flux(f);
        if (!std::isfinite(flux))
            throw std::invalid_argument("TVD secant: nonfinite face mass flux");
        const auto o=mesh.ownership().owner(f);
        if (o>=nc) throw std::invalid_argument("TVD secant: bad face owner");
        spectral[o]+=std::abs(flux);
        const auto n=mesh.ownership().neighbour(f);
        if (n>=0) {
            if (static_cast<std::size_t>(n)>=nc)
                throw std::invalid_argument("TVD secant: bad face neighbour");
            spectral[static_cast<std::size_t>(n)]+=std::abs(flux);
        }
    }
    auto* values=eq.matrix.values_data();
    const auto* cols=eq.matrix.columns_data();
    const auto* offsets=eq.matrix.row_offsets_data();
    TvdSecantImplicitResult result;
    for (std::size_t c=0;c<nc;++c) {
        const double velocity=u_current[c], old=u_previous[c];
        if (!std::isfinite(velocity) || !std::isfinite(old))
            throw std::invalid_argument("TVD secant: nonfinite velocity");
        const double delta=velocity-old;
        const double eps=controls.minimum_relative_displacement *
            std::max({1.0,std::abs(velocity),std::abs(old)});
        if (std::abs(delta)<=eps || spectral[c]<=0.0) continue;
        const double slope=-(b[c]-previous[c])/delta;
        if (!(slope>0.0) || !std::isfinite(slope)) continue;
        const double coeff=std::min(slope,controls.spectral_radius_cap*spectral[c]);
        bool diagonal_found=false;
        for (auto k=offsets[c];k<offsets[c+1];++k)
            if (cols[k]==c) {
                values[k]+=coeff;
                diagonal_found=true;
                break;
            }
        if (!diagonal_found) throw std::invalid_argument("TVD secant: missing CSR diagonal");
        eq.diagonal[c]+=coeff;
        eq.rhs(c)+=coeff*velocity;
        if (!std::isfinite(eq.diagonal[c]) || !std::isfinite(eq.rhs(c)))
            throw std::runtime_error("TVD secant: nonfinite regularized momentum matrix");
        result.regularized_cells++;
        result.max_diagonal_added=std::max(result.max_diagonal_added,coeff);
        result.max_ratio_to_face_flux=std::max(
            result.max_ratio_to_face_flux,coeff/spectral[c]);
    }
    return result;
}
} // namespace cfdx::physics

