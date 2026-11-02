// Independent physical benchmark: de Vahl Davis, square natural-convection cavity.
// Uses the actual production Boussinesq flow-energy fixed-point iteration.
// The benchmark must not be marked qualified unless Nu and velocity QoIs pass.
// Reference: G. de Vahl Davis, IJNMF 3 (1983), DOI 10.1002/fld.1650030305.
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/physics/production_energy_stage.h"
#include "cfdx/physics/production_radiation_contract.h"
#include "cfdx/physics/production_radiation_s2s_stage.h"
#include "cfdx/physics/production_buoyancy_stage.h"
#include "common/test_harness.h"

#include <cmath>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <sstream>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::io;
using namespace cfdx::physics;
using namespace cfdx::testing;

// ---------------------------------------------------------------------------
// Structured box mesh: nx x ny x 1 hexahedra of size dx x dy x 1 along x,
// y, z. Two point columns per (i, j) node (l = 0, 1 for z = 0, z = 1).
// Patches in order: left, right, bottom, top, front (z = 0), back (z = 1).
// Same construction as the stage-2 / stage-3 VVs so all tests exercise
// the identical production geometry path. Cell (i, j) has index
// i * ny + j.
// ---------------------------------------------------------------------------
static Mesh box_mesh(std::size_t nx, std::size_t ny, double dx, double dy)
{
    const auto pid = [&](std::size_t i, std::size_t j, std::size_t l) {
        return (i * (ny + 1) + j) * 2 + l;
    };
    const std::size_t n_points = 2 * (nx + 1) * (ny + 1);
    const std::size_t base_fy = (nx + 1) * ny;
    const std::size_t base_fz = base_fy + nx * (ny + 1);
    const std::size_t n_faces = base_fz + 2 * nx * ny;
    const auto fx = [&](std::size_t i, std::size_t j) { return i * ny + j; };
    const auto fy = [&](std::size_t i, std::size_t j) {
        return base_fy + i * (ny + 1) + j;
    };
    const auto fz = [&](std::size_t i, std::size_t j, std::size_t l) {
        return base_fz + 2 * (i * ny + j) + l;
    };
    const auto cell = [&](std::size_t i, std::size_t j) {
        return i * ny + j;
    };

    Mesh m;
    m.points().resize(n_points);
    for (std::size_t i = 0; i <= nx; ++i)
        for (std::size_t j = 0; j <= ny; ++j)
            for (std::size_t l = 0; l < 2; ++l)
                m.points().set(pid(i, j, l),
                               static_cast<double>(i) * dx,
                               static_cast<double>(j) * dy,
                               static_cast<double>(l));

    for (std::size_t i = 0; i <= nx; ++i)
        for (std::size_t j = 0; j < ny; ++j) {
            // Face normals point owner-to-neighbour for internal faces and
            // outward from the owner for boundary faces.
            if (i == 0)
                m.faces().push_face({pid(i, j, 0), pid(i, j, 1),
                                     pid(i, j + 1, 1), pid(i, j + 1, 0)});
            else
                m.faces().push_face({pid(i, j, 0), pid(i, j + 1, 0),
                                     pid(i, j + 1, 1), pid(i, j, 1)});
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j <= ny; ++j) {
            if (j == 0)
                m.faces().push_face({pid(i, j, 0), pid(i + 1, j, 0),
                                     pid(i + 1, j, 1), pid(i, j, 1)});
            else
                m.faces().push_face({pid(i, j, 0), pid(i, j, 1),
                                     pid(i + 1, j, 1), pid(i + 1, j, 0)});
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            for (std::size_t l = 0; l < 2; ++l) {
                if (l == 0)
                    m.faces().push_face({pid(i, j, l), pid(i, j + 1, l),
                                         pid(i + 1, j + 1, l), pid(i + 1, j, l)});
                else
                    m.faces().push_face({pid(i, j, l), pid(i + 1, j, l),
                                         pid(i + 1, j + 1, l), pid(i, j + 1, l)});
            }

    m.ownership().resize(n_faces);
    for (std::size_t j = 0; j < ny; ++j) {
        m.ownership().set_owner(fx(0, j), cell(0, j));
        m.ownership().set_neighbour(fx(0, j), FaceOwnership::BOUNDARY);
        m.ownership().set_owner(fx(nx, j), cell(nx - 1, j));
        m.ownership().set_neighbour(fx(nx, j), FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 1; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j) {
            m.ownership().set_owner(fx(i, j), cell(i - 1, j));
            m.ownership().set_neighbour(fx(i, j),
                                        static_cast<int>(cell(i, j)));
        }
    for (std::size_t i = 0; i < nx; ++i) {
        m.ownership().set_owner(fy(i, 0), cell(i, 0));
        m.ownership().set_neighbour(fy(i, 0), FaceOwnership::BOUNDARY);
        m.ownership().set_owner(fy(i, ny), cell(i, ny - 1));
        m.ownership().set_neighbour(fy(i, ny), FaceOwnership::BOUNDARY);
    }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 1; j < ny; ++j) {
            m.ownership().set_owner(fy(i, j), cell(i, j - 1));
            m.ownership().set_neighbour(fy(i, j),
                                        static_cast<int>(cell(i, j)));
        }
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
        for (std::size_t l = 0; l < 2; ++l) {
            m.ownership().set_owner(fz(i, j, l), cell(i, j));
            m.ownership().set_neighbour(fz(i, j, l),
                                        FaceOwnership::BOUNDARY);
        }

    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            m.cells().push_cell({fx(i, j), fx(i + 1, j),
                                 fy(i, j), fy(i, j + 1),
                                 fz(i, j, 0), fz(i, j, 1)});

    Patch left;
    left.name = "left";
    left.type = PatchType::WALL;
    for (std::size_t j = 0; j < ny; ++j) left.face_ids.push_back(fx(0, j));
    Patch right;
    right.name = "right";
    right.type = PatchType::WALL;
    for (std::size_t j = 0; j < ny; ++j) right.face_ids.push_back(fx(nx, j));
    Patch bottom;
    bottom.name = "bottom";
    bottom.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i) bottom.face_ids.push_back(fy(i, 0));
    Patch top;
    top.name = "top";
    top.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i) top.face_ids.push_back(fy(i, ny));
    Patch front;
    front.name = "front";
    front.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            front.face_ids.push_back(fz(i, j, 0));
    Patch back;
    back.name = "back";
    back.type = PatchType::WALL;
    for (std::size_t i = 0; i < nx; ++i)
        for (std::size_t j = 0; j < ny; ++j)
            back.face_ids.push_back(fz(i, j, 1));
    m.boundary().add_patch(left);
    m.boundary().add_patch(right);
    m.boundary().add_patch(bottom);
    m.boundary().add_patch(top);
    m.boundary().add_patch(front);
    m.boundary().add_patch(back);
    return m;
}

// ---------------------------------------------------------------------------
// Helpers. All flow solves use the production solver controls of the
// production driver (SIMPLE + upwind, strict continuity 1e-12, linear
// 1e-11 / 1000, pressure reference cell 0): the VV exercises the same
// numerical path, never a loosened test-only variant.
// ---------------------------------------------------------------------------

// Single complete material carrying the Boussinesq beta through the
// established extra_properties convention (strict case-double parse).
static BoundarySpec wall_thermal(const std::string& patch, double temperature)
{
    BoundarySpec bc;
    bc.patch_name = patch;
    bc.type = BCType::WALL;
    bc.value_type = BCValueType::WALL_THERMAL;
    bc.wall_temperature = temperature;
    return bc;
}

static VelocityBoundaryConditions no_slip_velocity_bcs(const Mesh& mesh)
{
    VelocityBoundaryConditions ubc;
    for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i)
        ubc[mesh.boundary().patch(i).name] = {
            VelocityBoundaryCondition::Type::FIXED_VALUE,
            {0.0, 0.0, 0.0}};
    return ubc;
}

static ScalarBoundaryConditions zero_gradient_pressure_bcs(const Mesh& mesh)
{
    ScalarBoundaryConditions pbc;
    for (std::size_t i = 0; i < mesh.boundary().n_patches(); ++i)
        pbc[mesh.boundary().patch(i).name] = {
            ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    return pbc;
}

static IncompressibleSolverControls vv_controls(double density,
                                                double kinematic_viscosity)
{
    IncompressibleSolverControls controls;
    controls.algorithm = PressureVelocityAlgorithm::SIMPLE;
    controls.convection_scheme = ConvectionScheme::UPWIND;
    controls.convergence.max_iterations = 2000;
    controls.convergence.continuity_tolerance = 1e-12;
    controls.linear_max_iterations = 1000;
    controls.linear_tolerance = 1e-11;
    controls.pressure_reference_cell = 0;
    controls.pressure_reference_value = 0.0;
    controls.density = density;
    controls.kinematic_viscosity = kinematic_viscosity;
    return controls;
}


namespace {
constexpr double kLength = 1.0;
constexpr double kPr = 0.71;
struct Benchmark {
    double ra;
    double nu_ref;
    double ux_ref;
    double uy_ref;
    std::vector<std::size_t> grids;
};

static Benchmark benchmark_for(double ra)
{
    if (ra == 1000.0) return {1000.0, 1.118, 3.649, 3.697, {8, 12, 16}};
    if (ra == 10000.0) return {10000.0, 2.243, 16.178, 19.617, {16, 24, 32}};
    throw std::invalid_argument("unsupported Rayleigh number: expected 1000 or 10000");
}
constexpr double kGravity = 9.81;
constexpr double kHot = 305.0;
constexpr double kCold = 295.0;

struct CavityMeasures {
    bool converged;
    double nu_hot;
    double nu_cold;
    double ux_max;
    double uy_max;
    double max_speed;
    double energy_balance;
    std::size_t outer_iterations;
};

CavityMeasures solve_case(std::size_t n, const Benchmark& benchmark)
{
    if (n < 4 || n > 128) throw std::invalid_argument("invalid cavity resolution");
    auto mesh = box_mesh(n, n, kLength / static_cast<double>(n),
                         kLength / static_cast<double>(n));
    if (!mesh.topo_validate().ok) throw std::runtime_error("invalid cavity mesh topology");
    const auto geometry = build_fv_geometry(mesh);
    // Reduce dimensional viscosity at Ra=1e4 to keep beta*DeltaT < 0.1.
    // The nondimensional Rayleigh/Prandtl benchmark is unchanged.
    const double nu = benchmark.ra >= 10000.0 ? 0.005 : 0.01;
    const double alpha = nu / kPr;
    const double beta = benchmark.ra * nu * alpha / (kGravity *
                         (kHot-kCold) * kLength*kLength*kLength);
    // beta*DeltaT << 0.1: Boussinesq approximation remains admissible.
    if (!(beta*(kHot-kCold) < 0.1)) throw std::runtime_error("Boussinesq invalid");
    MaterialSpec fluid;
    fluid.name="fluid";
    fluid.density=1.0;
    fluid.dynamic_viscosity=nu*fluid.density;
    fluid.specific_heat=1.0;
    fluid.thermal_conductivity=alpha*fluid.density*fluid.specific_heat;
    // Serialize the physical coefficient with full numerical precision.
    std::ostringstream beta_text;
    beta_text << std::setprecision(17) << beta;
    fluid.extra_properties["thermal_expansion_coefficient"]=beta_text.str();
    CaseSetup setup;
    setup.energy_model="buoyant";
    setup.ref_temperature=(kHot+kCold)/2.0;
    setup.gravity_vector="0 -9.81 0";
    setup.materials.push_back(fluid);
    setup.boundary_conditions.push_back(wall_thermal("left", kHot));
    setup.boundary_conditions.push_back(wall_thermal("right", kCold));
    setup.has_initial_condition=true;
    setup.initial_condition.temperature=setup.ref_temperature;
    Field<double,Location::CELL> velocity(mesh.n_cells(),"U","m/s",3);
    Field<double,Location::CELL> pressure(mesh.n_cells(),"p","Pa",1);
    auto controls=vv_controls(fluid.density, nu);
    const auto result=run_production_buoyant_outer_loop(
        mesh, no_slip_velocity_bcs(mesh),
        zero_gradient_pressure_bcs(mesh), controls, setup, velocity, pressure);
    CavityMeasures q{};
    q.converged = result.converged && result.energy.converged &&
        result.fixed_point_temperature_passed &&
        result.fixed_point_velocity_passed &&
        result.fixed_point_pressure_passed;
    q.outer_iterations=result.outer_iterations;
    q.energy_balance=result.energy.global_balance_normalized;

    // Independent hot/cold wall heat-flux reconstruction, face-to-cell distance,
    // not the nonlinear convergence report (Nu = q*L/[k*DeltaT*A]).
    // z-thickness is 1 m; hot/cold wall total area is L*1.
    auto nu_wall=[&](const std::string& name, double twall, bool hot) {
        double q_wall=0.0;
        for(std::size_t patch_id=0; patch_id<mesh.boundary().n_patches(); ++patch_id){
            const auto& patch=mesh.boundary().patch(patch_id);
            if(patch.name!=name) continue;
            for(const auto face: patch.face_ids){
                const auto owner=static_cast<std::size_t>(mesh.ownership().owner(face));
                const double area=geometry.face_area_vectors[face].mag();
                const double distance=(geometry.face_centres[face]-geometry.cell_centres[owner]).mag();
                if(!(area>0.0 && distance>0.0)) throw std::runtime_error("invalid wall geometry");
                const double dt= hot ? twall-result.temperature(owner) :
                                         result.temperature(owner)-twall;
                q_wall += fluid.thermal_conductivity*dt*area/distance;
            }
        }
        return q_wall*kLength /
            (fluid.thermal_conductivity*(kHot-kCold)*kLength*1.0);
    };
    q.nu_hot=nu_wall("left",kHot,true);
    q.nu_cold=nu_wall("right",kCold,false);
    // Cell-centred approximations along x=0.5 (Ux) and y=0.5 (Uy).
    // These are deliberately lower-accuracy velocity QoIs on coarse grids.
    for(std::size_t i=0;i<n;++i) for(std::size_t j=0;j<n;++j) {
        const std::size_t c=i*n+j;
        const auto& cc=geometry.cell_centres[c];
        const double ux=velocity.component_data(0)[c];
        const double uy=velocity.component_data(1)[c];
        const double uz=velocity.component_data(2)[c];
        if(!std::isfinite(ux) || !std::isfinite(uy) || !std::isfinite(result.temperature(c)))
            throw std::runtime_error("nonfinite Boussinesq field");
        q.max_speed=std::max(q.max_speed,std::sqrt(ux*ux+uy*uy+uz*uz));
        if (std::abs(cc.x-0.5)<=0.5/static_cast<double>(n)+1e-12)
            q.ux_max=std::max(q.ux_max,std::abs(ux)*kLength/alpha);
        if (std::abs(cc.y-0.5)<=0.5/static_cast<double>(n)+1e-12)
            q.uy_max=std::max(q.uy_max,std::abs(uy)*kLength/alpha);
    }
    std::cout << std::setprecision(12)
              << "CAVITY_DVD Ra=" << benchmark.ra << " Pr=" << kPr
              << " n=" << n << " converged=" << q.converged
              << " outer=" << q.outer_iterations
              << " Nu_hot=" << q.nu_hot << " Nu_cold=" << q.nu_cold
              << " ustar_max=" << q.ux_max << " vstar_max=" << q.uy_max
              << " speed_max=" << q.max_speed
              << " energy_balance=" << q.energy_balance
              << " reference_Nu=" << benchmark.nu_ref
              << " reference_u=" << benchmark.ux_ref
              << " reference_v=" << benchmark.uy_ref << '\n';
    return q;
}
} // namespace

int main(int argc, char** argv)
{
    // Default: three-grid prequalification. --grid selects a diagnostic
    // single-grid replay; it cannot by itself prove spatial convergence.
    bool single=false;
    std::size_t n=0;
    double ra=1000.0;
    for(int i=1;i<argc;++i){
        const std::string arg=argv[i];
        if(arg=="--grid" && i+1<argc){
            const std::string value=argv[++i];
            const auto number=std::stoul(value);
            if(number<4 || number>128 || std::to_string(number)!=value)
                throw std::invalid_argument("invalid --grid");
            n=number;
            single=true;
        } else if (arg=="--ra" && i+1<argc) {
            const std::string value=argv[++i];
            if(value=="1000") ra=1000.0;
            else if(value=="10000") ra=10000.0;
            else throw std::invalid_argument("--ra requires 1000 or 10000");
        } else throw std::invalid_argument("unknown cavity option: " + arg);
    }
    const auto benchmark=benchmark_for(ra);
    run_case("de_vahl_davis_buoyancy_three_grid_prequalification",[&] {
        std::vector<std::size_t> grids = single ? std::vector<std::size_t>{n} :
            benchmark.grids;
        std::vector<double> errors;
        for(const auto grid:grids) {
            const auto q=solve_case(grid, benchmark);
            if(!q.converged)
                throw std::runtime_error("full flow-energy fixed point did not converge");
            if(!(std::isfinite(q.nu_hot) && std::isfinite(q.nu_cold) &&
                 std::isfinite(q.ux_max) && std::isfinite(q.uy_max) &&
                 std::isfinite(q.energy_balance)))
                throw std::runtime_error("nonfinite physical qualification metrics");
            if(!(q.nu_hot>0.0 && q.nu_cold>0.0 &&
                 q.max_speed>1e-12 && q.energy_balance<=1e-8))
                throw std::runtime_error("circulation, wall flux or energy conservation failed");
            if(std::abs(q.nu_hot-q.nu_cold)/std::max(1.0,std::abs(q.nu_hot))>0.02)
                throw std::runtime_error("hot/cold wall heat-flux imbalance >2 percent");
            errors.push_back(std::abs(q.nu_hot-benchmark.nu_ref));
            // The finest grid is intentionally still a *prequalification*:
            // tighter, grid-converged velocity QoIs must be demonstrated on
            // additional 32/64 grids before this is promoted to full PASS.
            if(grid>=benchmark.grids.back()) {
                if(errors.back()/benchmark.nu_ref>0.05)
                    throw std::runtime_error("Nusselt benchmark error >5 percent");
                if(std::abs(q.ux_max-benchmark.ux_ref)/benchmark.ux_ref>0.12 ||
                   std::abs(q.uy_max-benchmark.uy_ref)/benchmark.uy_ref>0.12)
                    throw std::runtime_error("coarse-grid velocity benchmark error >12 percent");
            }
        }
        if(!single && !(errors[2]<errors[1] && errors[1]<errors[0]))
            throw std::runtime_error("three-grid physical Nu error not decreasing");
        std::cout << "CAVITY_DVD_STATUS: PREQUALIFICATION_PASS"
                  << " Ra=" << benchmark.ra << " grids=" << grids.size() << '\n';
    });
    return run_all();
}
