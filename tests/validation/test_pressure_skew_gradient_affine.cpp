#include "cfdx/physics/steady_incompressible_solver.h"
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
using namespace cfdx::core;
using namespace cfdx::physics;

namespace {
Mesh make_skewed_two_cell_mesh() {
    Mesh m;
    m.points().resize(12);
    const double xyz[12][3] = {
        {0,0,0},{1,0,0},{1.5,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1.5,1,1},{0,1,1},
        {2,0,0},{2,1,0},{2,0,1},{2,1,1}
    };
    for (std::size_t i=0;i<12;++i)
        m.points().set(i, xyz[i][0], xyz[i][1], xyz[i][2]);
    m.faces().push_face({0,3,7,4});
    m.faces().push_face({1,2,6,5});
    m.faces().push_face({8,9,11,10});
    m.faces().push_face({0,1,5,4});
    m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,3,2,1});
    m.faces().push_face({4,5,6,7});
    m.faces().push_face({1,8,10,5});
    m.faces().push_face({2,6,11,9});
    m.faces().push_face({1,2,9,8});
    m.faces().push_face({5,10,11,6});
    m.ownership().resize(11);
    for (std::size_t f=0; f<11; ++f) {
        m.ownership().set_owner(f, (f == 2 || f >= 7) ? 1 : 0);
        m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
    }
    m.ownership().set_owner(1,0);
    m.ownership().set_neighbour(1,1);
    m.cells().push_cell({0,1,3,4,5,6});
    m.cells().push_cell({1,2,7,8,9,10});
    for (std::size_t f=0; f<11; ++f) {
        if (f==1) continue;
        Patch b;
        b.name = "p" + std::to_string(f);
        b.type = PatchType::WALL;
        b.face_ids = {f};
        m.boundary().add_patch(b);
    }
    return m;
}
double analytic(const Vec3& v) {return 0.3+2.0*v.x-3.0*v.y;}
}
int main() {
    try {
        const auto mesh=make_skewed_two_cell_mesh();
        const auto g=build_fv_geometry(mesh);
        Field<double, Location::CELL> p(mesh.n_cells(),"p","Pa",1);
        for(std::size_t c=0;c<mesh.n_cells();++c) p(c)=analytic(g.cell_centres[c]);
        ScalarBoundaryConditions pbc;
        for (std::size_t f=0; f<mesh.n_faces();++f) {
            if(mesh.ownership().neighbour(f)>=0) continue;
            const auto patch=g.face_patch[f];
            pbc[mesh.boundary().patch(patch).name] = {
                ScalarBoundaryType::FIXED_VALUE,analytic(g.face_centres[f]),0};
        }
        const auto basic=gauss_gradient_with_boundary(p,mesh,g,pbc);
        const auto improved=gauss_gradient_with_boundary(p,mesh,g,pbc,true);
        double before=0, after=0;
        for(std::size_t c=0;c<mesh.n_cells();++c) {
            const Vec3 raw{basic.component_data(0)[c]-2.0,basic.component_data(1)[c]+3.0,basic.component_data(2)[c]};
            const Vec3 corrected{improved.component_data(0)[c]-2.0,improved.component_data(1)[c]+3.0,improved.component_data(2)[c]};
            before=std::max(before,raw.mag());
            after=std::max(after,corrected.mag());
        }
        std::cout << "SKEW_AFFINE_GRADIENT before=" << before << " after=" << after << "\n";
        if (!(before>1e-4 && after<before && after<1e-6))
            throw std::runtime_error("skew correction failed affine-exactness gate");
        // Differentiate the exact eight-sweep skew-corrected pressure
        // gradient. Old and new use identical fixed-value BCs, whereas the
        // Jacobian action must use their homogeneous variations.
        Field<double, Location::CELL> perturbation(mesh.n_cells(),"dp","Pa",1);
        for (std::size_t c=0; c<mesh.n_cells(); ++c)
            perturbation(c) = 0.3 - static_cast<double>(c);
        const auto jacobian = gauss_pressure_gradient_jacobian_action(
            perturbation,mesh,g,pbc,true);
        constexpr double step = 1e-5;
        auto perturbed = p;
        for (std::size_t c=0; c<mesh.n_cells(); ++c)
            perturbed(c) += step * perturbation(c);
        const auto updated = gauss_gradient_with_boundary(perturbed,mesh,g,pbc,true);
        for(std::size_t c=0; c<mesh.n_cells(); ++c)
            for(std::size_t component=0; component<3; ++component) {
                const double reference =
                    (updated.component_data(component)[c]-improved.component_data(component)[c])/step;
                const double analytic_jv = jacobian.component_data(component)[c];
                if (std::abs(reference-analytic_jv) >
                    1e-7 * std::max(1.0,std::abs(analytic_jv)))
                    throw std::runtime_error("skew gradient Jacobian action failed finite difference");
            }
        std::cout << "SKEW_GAUSS_JACOBIAN_DIRECTIONAL: PASS\n";
        // Use the production Rhie-Chow face reconstruction: this exercises
        // the same deferred face operator as the pressure correction.
        Field<double,Location::CELL> h(mesh.n_cells(),"h","m/s",3);
        h.fill(0.0);
        std::array<std::vector<double>,3> inverse_diagonal;
        for (auto& component : inverse_diagonal)
            component.assign(mesh.n_cells(),1.0);
        const VelocityBoundaryConditions velocity_bc;
        const auto flux_old=make_rhie_chow_flux_from_hbya(
            mesh,g,h,p,inverse_diagonal,1.0,velocity_bc,pbc);
        const auto flux_new=make_rhie_chow_flux_from_hbya(
            mesh,g,h,p,inverse_diagonal,1.0,velocity_bc,pbc,true);
        double internal_before=0.0, internal_after=0.0;
        for(std::size_t f=0; f<mesh.n_faces(); ++f) {
            if(mesh.ownership().neighbour(f)<0) continue;
            const auto o=mesh.ownership().owner(f);
            const auto n=static_cast<std::size_t>(mesh.ownership().neighbour(f));
            const double coeff=0.5*(g.cell_volumes[o]+g.cell_volumes[n]);
            const Vec3 gradient{2.0,-3.0,0.0};
            const double exact_flux=-coeff*gradient.dot(g.face_area_vectors[f]);
            internal_before=std::max(internal_before,std::abs(flux_old(f)-exact_flux));
            internal_after=std::max(internal_after,std::abs(flux_new(f)-exact_flux));
        }
        double boundary_error=0.0;
        for(std::size_t f=0; f<mesh.n_faces(); ++f) {
            if(mesh.ownership().neighbour(f)>=0) continue;
            const auto o=mesh.ownership().owner(f);
            const Vec3 gradient{2.0,-3.0,0.0};
            const double exact_flux=-g.cell_volumes[o]*gradient.dot(g.face_area_vectors[f]);
            boundary_error=std::max(boundary_error,std::abs(flux_new(f)-exact_flux));
        }
        if (!(boundary_error<1e-10))
            throw std::runtime_error("Rhie-Chow fixed-pressure boundary mismatch");
        std::cout << "SKEW_RHIE_CHOW boundary_error=" << boundary_error << "\n";
        std::cout << "SKEW_RHIE_CHOW internal_before=" << internal_before
                  << " internal_after=" << internal_after << "\n";
        if (!(internal_after<1e-10 && internal_after<=internal_before))
            throw std::runtime_error("Rhie-Chow affine pressure flux mismatch");
        IncompressibleSolverControls controls;
        controls.algorithm=PressureVelocityAlgorithm::COUPLED;
        controls.pressure_skew_corrected_gradient=true;
        bool rejected=false;
        try { validate_incompressible_controls(controls,mesh.n_cells()); }
        catch (const std::invalid_argument&) { rejected=true; }
        if (!rejected) throw std::runtime_error(
            "COUPLED with unassembled skew Jacobian must be rejected");
        std::cout << "PRESSURE_SKEW_GRADIENT_AFFINE: PASS\n";
    } catch (const std::exception& e) {
        std::cerr << "PRESSURE_SKEW_GRADIENT_AFFINE: FAIL " << e.what() << "\n";
        return 1;
    }
    return 0;
}
