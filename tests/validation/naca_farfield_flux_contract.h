#pragma once
// Independent kinematic flux consistency audit for prescribed-velocity faces.
// This is NOT a nonlinear, conservation, force or GCI certificate.
#include "cfdx/physics/finite_volume_transport.h"
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <algorithm>

namespace cfdx::validation {

struct NacaFarfieldFluxContract {
    std::size_t inspected_faces = 0;
    std::size_t worst_face = 0;
    double mismatch_signed = 0.0;
    double mismatch_l1 = 0.0;
    double mismatch_linf = 0.0;
    double prescribed_l1 = 0.0;
    double face_scale = 0.0;

    [[nodiscard]] double relative_l1() const {
        return mismatch_l1 / std::max(prescribed_l1, 1e-300);
    }
};

inline NacaFarfieldFluxContract audit_naca_fixed_velocity_farfield(
    const cfdx::core::Mesh& mesh,
    const cfdx::physics::FvGeometry& geometry,
    const cfdx::core::Field<double,cfdx::core::Location::FACE>& phi,
    const std::string& patch_name,
    const cfdx::core::Vec3& fixed_velocity, double density) {
    using namespace cfdx::core;
    if (!std::isfinite(density) || density <= 0.0 ||
        !std::isfinite(fixed_velocity.x) || !std::isfinite(fixed_velocity.y) ||
        !std::isfinite(fixed_velocity.z))
        throw std::invalid_argument("NACA fixed-velocity flux audit: invalid physical data");
    if (phi.size() != mesh.n_faces() ||
        geometry.face_area_vectors.size() != mesh.n_faces())
        throw std::invalid_argument("NACA fixed-velocity flux audit: face count mismatch");
    const auto id=mesh.boundary().find(patch_name);
    if (id>=mesh.boundary().n_patches())
        throw std::invalid_argument("NACA fixed-velocity flux audit: missing patch");
    const auto& patch=mesh.boundary().patch(id);
    if (patch.face_ids.empty())
        throw std::invalid_argument("NACA fixed-velocity flux audit: empty patch");
    NacaFarfieldFluxContract audit;
    for (const auto f:patch.face_ids) {
        if (f>=mesh.n_faces() || mesh.ownership().neighbour(f)>=0)
            throw std::invalid_argument("NACA fixed-velocity flux audit: nonboundary patch face");
        const Vec3 sf=geometry.face_area_vectors[f];
        const double expected=density*(sf.x*fixed_velocity.x+
                                       sf.y*fixed_velocity.y+sf.z*fixed_velocity.z);
        const double actual=phi(f);
        const double scale=density*sf.mag()*fixed_velocity.mag();
        if (!std::isfinite(expected) || !std::isfinite(actual) ||
            !std::isfinite(scale) || !(sf.mag()>0.0))
            throw std::invalid_argument("NACA fixed-velocity flux audit: nonfinite/degenerate face");
        const double error=actual-expected;
        const double abs_error=std::abs(error);
        audit.mismatch_signed+=error;
        audit.mismatch_l1+=abs_error;
        audit.prescribed_l1+=std::abs(expected);
        audit.face_scale=std::max(audit.face_scale,scale);
        if (abs_error>audit.mismatch_linf) {
            audit.mismatch_linf=abs_error;
            audit.worst_face=f;
        }
        ++audit.inspected_faces;
    }
    return audit;
}

inline bool naca_kinematic_face_flux_matches(
    const NacaFarfieldFluxContract& audit, double relative_tolerance=1e-10) {
    if (!std::isfinite(relative_tolerance) || relative_tolerance<0.0 ||
        audit.inspected_faces==0 || !std::isfinite(audit.mismatch_linf))
        return false;
    // Scale by the physical velocity and face area. A global signed balance
    // can be almost zero even if individual boundary faces are inconsistent.
    return audit.mismatch_linf <= relative_tolerance*std::max(audit.face_scale,1e-300);
}
} // namespace cfdx::validation
