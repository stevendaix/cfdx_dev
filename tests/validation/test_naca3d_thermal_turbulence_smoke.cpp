// NACA 3-D extrusion operator integration smoke, NOT an aerodynamic RANS/CHT
// qualification. No velocity is solved; mass flux is identically zero.
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/energy_solver.h"
#include "cfdx/physics/turbulence_solver.h"
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;

int main(int argc, char** argv) {
    if (argc != 2) throw std::invalid_argument("usage: test_naca3d_thermal_turbulence_smoke mesh.h5");
    Mesh mesh;
    if (!cfdx::io::read_mesh_hdf5(argv[1], mesh) || !mesh.topo_validate().ok)
        throw std::runtime_error("NACA3D: invalid HDF5 mesh");
    const auto geometry = build_fv_geometry(mesh);
    const auto n = mesh.n_cells();
    Field<double, Location::FACE> mass_flux(mesh.n_faces(), "flux", "kg/s", 1);
    mass_flux.fill(0.0);

    // Energy: exact isothermal equilibrium T=300 K with no source or flux,
    // regardless of the 3-D geometry or thermal conductivity.
    Field<double, Location::CELL> temperature(n, "T", "K", 1);
    Field<double, Location::CELL> source(n, "Q", "W/m3", 1);
    temperature.fill(300.0); source.fill(0.0);
    ScalarBoundaryConditions tbc;
    tbc["airfoil"] = {ScalarBoundaryType::FIXED_VALUE, 300.0, 0.0};
    tbc["farfield"] = {ScalarBoundaryType::FIXED_VALUE, 300.0, 0.0};
    tbc["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    tbc["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    EnergySolverControls ec;
    ec.conductivity = 2.0; ec.relaxation = 1.0;
    ec.max_iterations = 6; ec.tolerance = 1e-9;
    const auto energy = solve_energy(mesh, geometry, mass_flux,
                                      temperature, source, ec, tbc);
    double tmax = 0.0;
    for (std::size_t c = 0; c < n; ++c) {
        if (!std::isfinite(temperature(c))) throw std::runtime_error("NACA3D: nonfinite temperature");
        tmax = std::max(tmax, std::abs(temperature(c) - 300.0));
    }
    std::cout << "NACA3D_THERMAL cells=" << n << " converged=" << energy.converged
              << " Linf_T=" << tmax << '\n';
    if (!energy.converged || tmax > 1e-7)
        throw std::runtime_error("NACA3D: isothermal thermal MMS failed");

    // Turbulence: one transport step on the same volume mesh; bounded fields
    // and finite residuals are a smoke condition, NOT validation against TMR.
    Field<double, Location::CELL> k(n, "k", "m2/s2", 1);
    Field<double, Location::CELL> omega(n, "omega", "1/s", 1);
    Field<double, Location::CELL> strain(n, "strain", "1/s", 1);
    k.fill(0.01); omega.fill(1.0); strain.fill(0.5);
    TurbulenceTransportControls tc;
    tc.model = TurbulenceModel::KOMEGA;
    tc.molecular_viscosity = 0.001;
    ScalarBoundaryConditions kb, wb;
    kb["airfoil"] = {ScalarBoundaryType::FIXED_VALUE, 1e-6, 0.0};
    kb["farfield"] = {ScalarBoundaryType::FIXED_VALUE, 0.01, 0.0};
    wb["airfoil"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
    wb["farfield"] = {ScalarBoundaryType::FIXED_VALUE, 1.0, 0.0};
    for (auto* bc : {&kb, &wb}) {
        (*bc)["front"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
        (*bc)["back"] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};
    }
    const auto turbulence = solve_komega_transport(mesh, geometry, mass_flux,
        k, omega, strain, tc, kb, wb, 1, 1e-7);
    for (std::size_t c = 0; c < n; ++c) {
        if (!std::isfinite(k(c)) || !std::isfinite(omega(c)) ||
            k(c) < tc.k_min || omega(c) < tc.omega_min)
            throw std::runtime_error("NACA3D: invalid turbulence field");
    }
    if (!std::isfinite(turbulence.k_residual) ||
        !std::isfinite(turbulence.second_residual) ||
        turbulence.nonfinite_bound_events != 0)
        throw std::runtime_error("NACA3D: invalid turbulence residual");
    std::cout << "NACA3D_TURBULENCE cells=" << n
              << " steps=" << turbulence.iterations
              << " k_residual=" << turbulence.k_residual
              << " omega_residual=" << turbulence.second_residual
              << " nonfinite=" << turbulence.nonfinite_bound_events << '\n';
    return 0;
}
