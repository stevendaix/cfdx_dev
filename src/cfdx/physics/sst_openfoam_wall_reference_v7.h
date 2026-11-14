#pragma once

// V7 reference-only evaluator of OpenFOAM Foundation omegaWallFunction's
// pointwise viscous/log switch or laminar-weighted blend. This intentionally
// DOES NOT replace CFDX's SST boundary discretisation: OpenFOAM constrains
// owner-cell omega AND production G, and a comparable nut wall function is
// necessary. Using the pointwise formula alone cannot qualify BFS/SST.

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace cfdx::physics {

enum class ReferenceOmegaWallBlend { Stepwise, LaminarWeighted };

struct OpenFoamOmegaWallPoint {
    double reynolds_y = 0.0;
    double y_plus_from_k = 0.0;
    double laminar_fraction = 0.0;
    double viscous_omega = 0.0;
    double logarithmic_omega = 0.0;
    double reference_omega = 0.0;
};

// Formula corresponds to the face-local omega quantity (not the full
// patch-cell manipulation) from OpenFOAM Foundation omegaWallFunction v13.
// wall_normal_velocity_gradient is mag(snGrad(U_wall)) [1/s], and is not
// interchangeable with the first-cell viscous y+ proxy in CFDX.
inline OpenFoamOmegaWallPoint openfoam_omega_wall_reference(
    double k, double molecular_kinematic_viscosity,
    double face_owner_normal_distance,
    double wall_normal_velocity_gradient,
    ReferenceOmegaWallBlend blend=ReferenceOmegaWallBlend::Stepwise,
    double cmu=0.09, double beta1=0.075, double kappa=0.41,
    double laminar_yplus=11.06)
{
    if (!(k>=0.0) || !std::isfinite(k) ||
        !(molecular_kinematic_viscosity>0.0) ||
        !std::isfinite(molecular_kinematic_viscosity) ||
        !(face_owner_normal_distance>0.0) ||
        !std::isfinite(face_owner_normal_distance) ||
        !(wall_normal_velocity_gradient>=0.0) ||
        !std::isfinite(wall_normal_velocity_gradient) ||
        !(cmu>0.0) || !std::isfinite(cmu) ||
        !(beta1>0.0) || !std::isfinite(beta1) ||
        !(kappa>0.0) || !std::isfinite(kappa) ||
        !(laminar_yplus>0.0) || !std::isfinite(laminar_yplus))
        throw std::invalid_argument("OpenFOAM omega wall reference: invalid inputs");
    OpenFoamOmegaWallPoint result;
    const double y=face_owner_normal_distance;
    const double nu=molecular_kinematic_viscosity;
    const double sqrt_cmu=std::sqrt(cmu);
    result.reynolds_y=y*std::sqrt(k)/nu;
    result.y_plus_from_k=std::sqrt(sqrt_cmu)*result.reynolds_y;
    result.viscous_omega=6.0*nu/(beta1*y*y);
    if (blend==ReferenceOmegaWallBlend::Stepwise) {
        result.laminar_fraction=result.y_plus_from_k < laminar_yplus ? 1.0 : 0.0;
        // OpenFOAM's stepwise turbulent branch uses uStar=sqrt(Cmu5*k).
        // This differs from the blended branch below at finite wall shear.
        const double u_star=std::sqrt(sqrt_cmu*k);
        result.logarithmic_omega=u_star/(sqrt_cmu*kappa*y);
    } else {
        // Foundation v13 blended branch uses Rey, not yPlus, in exp(-Rey/11).
        result.laminar_fraction=std::exp(-result.reynolds_y/11.0);
        const double u_star=std::sqrt(
            result.laminar_fraction*nu*wall_normal_velocity_gradient +
            (1.0-result.laminar_fraction)*sqrt_cmu*k);
        result.logarithmic_omega=u_star/(sqrt_cmu*kappa*y);
    }
    result.reference_omega=result.laminar_fraction*result.viscous_omega+
        (1.0-result.laminar_fraction)*result.logarithmic_omega;
    if (!std::isfinite(result.reference_omega) ||
        !(result.reference_omega>0.0) ||
        !std::isfinite(result.y_plus_from_k) ||
        !std::isfinite(result.viscous_omega) ||
        !std::isfinite(result.logarithmic_omega))
        throw std::domain_error("OpenFOAM omega wall reference: overflow");
    return result;
}

} // namespace cfdx::physics
