// Audit : reconstruction du cisaillement pariétal NACA0012 (GG vs WLS vs 1D).
//
// Test d'audit de POST-TRAITEMENT : aucune méthode numérique du solveur n'est
// modifiée. Trois reconstructions du cisaillement pariétal visqueux sont
// comparées sur les faces du patch airfoil :
//
//   gg    : gradient de Gauss avec conditions aux limites (statu quo du test
//           de qualification NACA0012), traction visqueuse complète
//           tau = mu (grad u + grad u^T) évaluée aux centres des cellules
//           adjacentes à la paroi.
//   wls   : moindres carrés pondérés (compute_gradient_weighted_least_squares),
//           ghosts DIRICHLET par face sur airfoil/farfield, même traction.
//   one_d : reconstruction tangentielle 1D mu |u_t| / d_n par face.
//
// Partie M (quick, sans solveur) : MMS pariétal sur le canal channel_n16.h5.
//   M1 : u = sin(pi y), oracle |du/dy| = pi sur les deux parois.
//   M2 : u = y (1 - y), oracle |du/dy| = 1 sur les deux parois.
//   Attendus analytiques (ny = 16, dy = 1/16, centre de la 1re cellule à
//   dy/2 de la paroi) :
//     M1 : gg/wls = 1.117 %, one_d = 0.156 % -> gate 5 % (trois méthodes)
//     M2 : gg/wls = 7.812 % (biais centre-cellule sur profil quadratique),
//          one_d = 3.125 % -> gate 12 % (gg/wls), 5 % (one_d)
//
// Partie N (long) : NACA0012 Re = 1000, alpha = 0, n64/n128/n256, config de
// solveur IDENTIQUE au test de qualification (copie verbatim). L'écart relatif
// entre les trois reconstructions (rel_spread = (max-min)/mean des Cd_visc)
// est LA MESURE : il n'est PAS gated. Les gates sont uniquement défensifs :
// finitude, Cd_visc > 0, continuité < 1e-7, convergence, d_n > 0.

#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/core/numerics/gradient.h"
#include "cfdx/core/numerics/gradient_stencil.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

constexpr double kPi = 3.14159265358979323846;

struct WallShearAudit {
    double cd_visc_gg = 0.0;
    double cd_visc_wls = 0.0;
    double cd_visc_1d = 0.0;
    double rel_spread = 0.0;
    double dn_min = 0.0;
    double dn_mean = 0.0;
    double continuity = 0.0;
    std::size_t iterations = 0;
};

// Taille du vecteur de conditions par face : max(face id)+1 sur tous les
// patches (le reconstructeur n'indexe que des faces de frontière).
std::size_t boundary_bc_size(const Mesh& mesh)
{
    std::size_t n = 0;
    for (std::size_t p = 0; p < mesh.boundary().n_patches(); ++p)
        for (const auto f : mesh.boundary().patch(p).face_ids)
            n = std::max(n, f + 1);
    return n;
}

// ---------------------------------------------------------------------------
// Partie M : MMS pariétal sur le canal (quick, sans solveur)
// ---------------------------------------------------------------------------

struct MmsWallShear {
    double gg = 0.0;
    double wls = 0.0;
    double one_d = 0.0;
};

MmsWallShear wall_shear_mms(const std::string& path, int field_id)
{
    if (field_id != 1 && field_id != 2)
        throw std::invalid_argument("wall_shear_mms: unknown MMS field id");
    Mesh mesh;
    if (!cfdx::io::read_mesh_hdf5(path, mesh))
        throw std::runtime_error("cannot read channel mesh: " + path);
    if (!mesh.topo_validate().ok)
        throw std::runtime_error("invalid channel topology: " + path);
    const auto geometry = build_fv_geometry(mesh);
    const std::size_t n_cells = mesh.n_cells();

    const auto bottom_id = mesh.boundary().find("bottom");
    const auto top_id = mesh.boundary().find("top");
    if (bottom_id >= mesh.boundary().n_patches() ||
        top_id >= mesh.boundary().n_patches())
        throw std::runtime_error("channel wall patches (bottom/top) are missing");

    Field<double, Location::CELL> u(n_cells, "mms_u", "m/s", 1);
    for (std::size_t c = 0; c < n_cells; ++c) {
        const double y = geometry.cell_centres[c].y;
        u(c) = (field_id == 1) ? std::sin(kPi * y) : y * (1.0 - y);
    }

    // gg : Gauss + CL par patch. Paroi = DIRICHLET u = 0 (exact pour les deux
    // champs MMS) ; les autres faces en zéro-gradient n'affectent pas la
    // composante normale du gradient aux cellules de paroi.
    ScalarBoundaryConditions gbc;
    gbc["bottom"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    gbc["top"] = {ScalarBoundaryType::FIXED_VALUE, 0.0, 0.0};
    gbc["inlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    gbc["outlet"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    gbc["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    gbc["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    const auto g_gg = gauss_gradient_with_boundary(u, mesh, geometry, gbc);

    // wls : ghosts DIRICHLET u = 0 sur les faces de paroi. Les faces
    // front/back (EMPTY) gardent le zéro-gradient par défaut, qui fournit le
    // rang 3 en z sur ce maillage à couche z unique.
    std::vector<BoundaryGradientCondition> bc(boundary_bc_size(mesh));
    for (const auto f : mesh.boundary().patch(bottom_id).face_ids)
        bc[f] = BoundaryGradientCondition{BoundaryGradientConditionType::DIRICHLET, 0.0};
    for (const auto f : mesh.boundary().patch(top_id).face_ids)
        bc[f] = BoundaryGradientCondition{BoundaryGradientConditionType::DIRICHLET, 0.0};
    const auto g_wls = compute_gradient_weighted_least_squares(
        u, mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
        std::numeric_limits<double>::infinity(),
        BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, &bc);

    // est = -(grad u . n), n sortant du fluide ; one_d = |u| / d_n (la
    // vitesse MMS est purement tangentielle aux parois horizontales).
    const double oracle = (field_id == 1) ? kPi : 1.0;
    MmsWallShear out;
    const auto consider = [&](const auto& wall_id) {
        for (const auto f : mesh.boundary().patch(wall_id).face_ids) {
            const std::size_t o = mesh.ownership().owner(f);
            const Vec3 Sf = geometry.face_area_vectors[f];
            const Vec3 n = Sf.normalized();
            const double d_n = (geometry.face_centres[f] - geometry.cell_centres[o]).dot(n);
            if (!(d_n > 0.0) || !std::isfinite(d_n))
                throw std::runtime_error("wall_shear_mms: non-positive wall distance");
            const auto grad_err = [&](const Field<double, Location::CELL>& g) {
                const double est = -(g.component_data(0)[o] * n.x
                                    + g.component_data(1)[o] * n.y
                                    + g.component_data(2)[o] * n.z);
                return std::abs(est - oracle) / oracle;
            };
            out.gg = std::max(out.gg, grad_err(g_gg));
            out.wls = std::max(out.wls, grad_err(g_wls));
            out.one_d = std::max(out.one_d,
                                std::abs(std::abs(u(o)) / d_n - oracle) / oracle);
        }
    };
    consider(bottom_id);
    consider(top_id);
    return out;
}

void run_channel_mms(const std::string& path)
{
    struct Spec { const char* field; int id; double gate_grad; double gate_one_d; };
    const Spec specs[] = {
        {"sin_pi_y", 1, 5e-2, 5e-2},
        {"y_1_minus_y", 2, 12e-2, 5e-2},
    };
    for (const auto& spec : specs) {
        const MmsWallShear r = wall_shear_mms(path, spec.id);
        const double oracle = (spec.id == 1) ? kPi : 1.0;
        std::cout << "WALL_SHEAR_MMS field=" << spec.field << " method=gg"
                  << " max_rel_err=" << r.gg << " oracle=" << oracle << "\n";
        std::cout << "WALL_SHEAR_MMS field=" << spec.field << " method=wls"
                  << " max_rel_err=" << r.wls << " oracle=" << oracle << "\n";
        std::cout << "WALL_SHEAR_MMS field=" << spec.field << " method=one_d"
                  << " max_rel_err=" << r.one_d << " oracle=" << oracle << "\n";
        if (!std::isfinite(r.gg) || !std::isfinite(r.wls) || !std::isfinite(r.one_d))
            throw std::runtime_error(std::string("wall-shear MMS non-finite error, field=") + spec.field);
        if (r.gg > spec.gate_grad || r.wls > spec.gate_grad)
            throw std::runtime_error(std::string("wall-shear MMS gg/wls gate failed, field=") + spec.field);
        if (r.one_d > spec.gate_one_d)
            throw std::runtime_error(std::string("wall-shear MMS one_d gate failed, field=") + spec.field);
    }
}

// ---------------------------------------------------------------------------
// Partie N : NACA0012 Re = 1000, alpha = 0, trois niveaux de maillage
// ---------------------------------------------------------------------------

WallShearAudit run_naca(const std::string& path, std::size_t level)
{
    Mesh mesh;
    if (!cfdx::io::read_mesh_hdf5(path, mesh))
        throw std::runtime_error("cannot read NACA0012 mesh: " + path);
    if (!mesh.topo_validate().ok)
        throw std::runtime_error("invalid NACA0012 topology: " + path);

    constexpr double rho = 1.0, Uinf = 1.0, Re = 1000.0, chord = 1.0;
    constexpr double nu = Uinf * chord / Re;
    Field<double, Location::CELL> U(mesh.n_cells(), "U", "m/s", 3);
    Field<double, Location::CELL> p(mesh.n_cells(), "p", "Pa", 1);
    U.fill(0.0);
    p.fill(0.0);
    VelocityBoundaryConditions ubc;
    ubc["farfield"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {Uinf, 0, 0}};
    ubc["airfoil"] = {VelocityBoundaryCondition::Type::FIXED_VALUE, {0, 0, 0}};
    ubc["front"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ubc["back"] = {VelocityBoundaryCondition::Type::ZERO_GRADIENT, {0, 0, 0}};
    ScalarBoundaryConditions pbc;
    pbc["farfield"] = {ScalarBoundaryType::FIXED_VALUE, 0, 0};
    pbc["airfoil"] = {ScalarBoundaryType::ZERO_GRADIENT, 0, 0};
    pbc["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0, 0};
    pbc["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0, 0};

    // Config de solveur verbatim du test de qualification NACA0012.
    IncompressibleSolverControls c;
    c.algorithm = PressureVelocityAlgorithm::SIMPLE;
    c.density = rho;
    c.kinematic_viscosity = nu;
    c.linear_max_iterations = 4000;
    c.linear_tolerance = 1e-10;
    c.convergence.max_iterations = 8000;
    c.convergence.relative_tolerance = 1e-8;
    c.convergence.continuity_tolerance = 1e-8;
    c.coupling.alpha_u = 0.6;
    c.coupling.alpha_p = 0.25;
    c.use_bounded_convection = true;
    c.convection_scheme = ConvectionScheme::UPWIND;
    c.pressure_reference_cell = 0;
    c.pressure_reference_value = 0.0;
    c.diagnostics.iteration_trace = true;
    c.diagnostics.iteration_trace_frequency = 200;

    const auto geometry = build_fv_geometry(mesh);
    for (std::size_t i = 0; i < mesh.n_cells(); ++i)
        U.component_data(0)[i] = Uinf;

    const auto solve = solve_steady_incompressible(mesh, U, p, ubc, pbc, c);
    if (!solve.converged || solve.history.empty())
        throw std::runtime_error("NACA0012 audit solve did not converge");
    WallShearAudit out;
    out.continuity = solve.history.back().continuity_linf;
    out.iterations = solve.iterations;

    const std::size_t n_cells = mesh.n_cells();
    Field<double, Location::CELL> ux(n_cells, "ux", "m/s", 1), uy(n_cells, "uy", "m/s", 1);
    for (std::size_t i = 0; i < n_cells; ++i) {
        ux(i) = U.component_data(0)[i];
        uy(i) = U.component_data(1)[i];
    }

    // gg : conditions gradient verbatim du test de qualification.
    ScalarBoundaryConditions gx, gy;
    gx["farfield"] = {ScalarBoundaryType::FIXED_VALUE, Uinf, 0};
    gx["airfoil"] = {ScalarBoundaryType::FIXED_VALUE, 0, 0};
    gx["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0, 0};
    gx["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0, 0};
    gy["farfield"] = {ScalarBoundaryType::FIXED_VALUE, 0, 0};
    gy["airfoil"] = {ScalarBoundaryType::FIXED_VALUE, 0, 0};
    gy["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0, 0};
    gy["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0, 0};
    const auto gux_gg = gauss_gradient_with_boundary(ux, mesh, geometry, gx);
    const auto guy_gg = gauss_gradient_with_boundary(uy, mesh, geometry, gy);

    // wls : ghosts DIRICHLET par face — airfoil 0, farfield Uinf (ux) / 0
    // (uy), front/back zéro-gradient par défaut.
    std::vector<BoundaryGradientCondition> bcx(boundary_bc_size(mesh));
    std::vector<BoundaryGradientCondition> bcy(boundary_bc_size(mesh));
    const auto airfoil_id = mesh.boundary().find("airfoil");
    const auto farfield_id = mesh.boundary().find("farfield");
    if (airfoil_id >= mesh.boundary().n_patches() ||
        farfield_id >= mesh.boundary().n_patches())
        throw std::runtime_error("NACA0012 required boundary patches are missing");
    for (const auto f : mesh.boundary().patch(airfoil_id).face_ids) {
        bcx[f] = BoundaryGradientCondition{BoundaryGradientConditionType::DIRICHLET, 0.0};
        bcy[f] = BoundaryGradientCondition{BoundaryGradientConditionType::DIRICHLET, 0.0};
    }
    for (const auto f : mesh.boundary().patch(farfield_id).face_ids) {
        bcx[f] = BoundaryGradientCondition{BoundaryGradientConditionType::DIRICHLET, Uinf};
        bcy[f] = BoundaryGradientCondition{BoundaryGradientConditionType::DIRICHLET, 0.0};
    }
    const auto gux_wls = compute_gradient_weighted_least_squares(
        ux, mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
        std::numeric_limits<double>::infinity(),
        BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, &bcx);
    const auto guy_wls = compute_gradient_weighted_least_squares(
        uy, mesh, GradientWeighting::INVERSE_DISTANCE_SQUARED,
        std::numeric_limits<double>::infinity(),
        BoundaryGradientPolicy::ZERO_GRADIENT_GHOST, &bcy);

    // Intégration de la traction visqueuse (composante x) sur les faces du
    // patch airfoil. Formule verbatim du test de qualification pour gg.
    const double q_dyn = 0.5 * rho * Uinf * Uinf;
    double fx_gg = 0.0, fx_wls = 0.0, fx_1d = 0.0;
    double dn_min = std::numeric_limits<double>::infinity();
    double dn_sum = 0.0;
    std::size_t n_wall = 0;
    const auto& patch = mesh.boundary().patch(airfoil_id);
    for (const auto f : patch.face_ids) {
        const std::size_t o = mesh.ownership().owner(f);
        const auto Sf = geometry.face_area_vectors[f];
        const auto n = Sf.normalized();
        const double d_n = (geometry.face_centres[f] - geometry.cell_centres[o]).dot(n);
        if (!(d_n > 0.0) || !std::isfinite(d_n))
            throw std::runtime_error("NACA0012 audit: non-positive wall distance");
        const double txx_gg = 2.0 * rho * nu * gux_gg.component_data(0)[o];
        const double txy_gg = rho * nu * (gux_gg.component_data(1)[o] + guy_gg.component_data(0)[o]);
        fx_gg += -(txx_gg * n.x + txy_gg * n.y) * Sf.mag();
        const double txx_wls = 2.0 * rho * nu * gux_wls.component_data(0)[o];
        const double txy_wls = rho * nu * (gux_wls.component_data(1)[o] + guy_wls.component_data(0)[o]);
        fx_wls += -(txx_wls * n.x + txy_wls * n.y) * Sf.mag();
        const Vec3 uc{U.component_data(0)[o], U.component_data(1)[o], U.component_data(2)[o]};
        const Vec3 ut = uc - n * uc.dot(n);
        fx_1d += rho * nu * (ut.x / d_n) * Sf.mag();
        dn_min = std::min(dn_min, d_n);
        dn_sum += d_n;
        ++n_wall;
    }
    if (n_wall == 0)
        throw std::runtime_error("NACA0012 audit: airfoil patch has no faces");

    out.cd_visc_gg = fx_gg / q_dyn;
    out.cd_visc_wls = fx_wls / q_dyn;
    out.cd_visc_1d = fx_1d / q_dyn;
    out.dn_min = dn_min;
    out.dn_mean = dn_sum / static_cast<double>(n_wall);
    const double lo = std::min(out.cd_visc_gg, std::min(out.cd_visc_wls, out.cd_visc_1d));
    const double hi = std::max(out.cd_visc_gg, std::max(out.cd_visc_wls, out.cd_visc_1d));
    const double mean = (out.cd_visc_gg + out.cd_visc_wls + out.cd_visc_1d) / 3.0;
    out.rel_spread = (hi - lo) / mean;

    // Gates défensifs uniquement : la spread est la mesure, jamais un gate.
    if (!std::isfinite(out.cd_visc_gg) || !std::isfinite(out.cd_visc_wls) ||
        !std::isfinite(out.cd_visc_1d) || !std::isfinite(out.rel_spread))
        throw std::runtime_error("NACA0012 audit: non-finite Cd_visc");
    if (!(out.cd_visc_gg > 0.0) || !(out.cd_visc_wls > 0.0) || !(out.cd_visc_1d > 0.0))
        throw std::runtime_error("NACA0012 audit: Cd_visc not strictly positive");
    if (!(out.continuity < 1e-7))
        throw std::runtime_error("NACA0012 audit: continuity gate failed");

    std::cout << "NACA0012_WALL_SHEAR_AUDIT level=" << level
              << " Cd_visc_gg=" << out.cd_visc_gg
              << " Cd_visc_wls=" << out.cd_visc_wls
              << " Cd_visc_1d=" << out.cd_visc_1d
              << " rel_spread=" << out.rel_spread
              << " dn_min=" << out.dn_min
              << " dn_mean=" << out.dn_mean
              << " continuity=" << out.continuity
              << " iterations=" << out.iterations << "\n";
    return out;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const bool quick = argc >= 2 && std::string(argv[1]) == "--quick";
        if (quick) {
            if (argc != 3)
                throw std::invalid_argument(
                    "usage: test_naca0012_wall_shear_audit --quick channel_n16.h5");
            run_channel_mms(argv[2]);
            std::cout << "NACA0012_WALL_SHEAR_AUDIT_QUICK: PASS mms_wall_shear_gates_ok\n";
            return 0;
        }
        if (argc != 5)
            throw std::invalid_argument(
                "usage: test_naca0012_wall_shear_audit channel_n16.h5"
                " naca0012_laminar_n64.h5 naca0012_laminar_n128.h5"
                " naca0012_laminar_n256.h5");
        run_channel_mms(argv[1]);
        for (std::size_t i = 0; i < 3; ++i)
            run_naca(argv[i + 2], 64u << i);
        std::cout << "NACA0012_WALL_SHEAR_AUDIT: DIAGNOSTIC spread_is_measurement_not_gated\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "NACA0012_WALL_SHEAR_AUDIT: FAIL: " << e.what() << "\n";
        return 1;
    }
}
