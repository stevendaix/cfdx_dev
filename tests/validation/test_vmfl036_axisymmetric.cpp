#include "cfdx/physics/axisymmetric_vmfl036_solver.h"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#ifndef CFDX_SOURCE_DIR
#error "CFDX_SOURCE_DIR must be defined"
#endif
#ifndef CFDX_PYTHON_EXECUTABLE
#define CFDX_PYTHON_EXECUTABLE "python3"
#endif

using cfdx::physics::axisymmetric::AxisymControls;
using cfdx::physics::axisymmetric::AxisymMesh;
using cfdx::physics::axisymmetric::AxisymResult;
using cfdx::physics::axisymmetric::VMFL036AxisymmetricSolver;

namespace {
struct Level { int nr; int nt; const char* name; };

void generate(const std::filesystem::path& out, const Level& l)
{
    const auto script=std::filesystem::path(CFDX_SOURCE_DIR)/"scripts"/"generate_vmfl036_axisymmetric_mesh.py";
    std::ostringstream cmd;
    cmd << "'" << CFDX_PYTHON_EXECUTABLE << "' '" << script.string()
        << "' --nr " << l.nr << " --nt " << l.nt << " --output '" << out.string() << "'";
    if(std::system(cmd.str().c_str())!=0) throw std::runtime_error("axisymmetric mesh generation failed");
}

AxisymResult run_level(const Level& l, double mu, double cd_ref, const char* tag)
{
    const auto path=std::filesystem::temp_directory_path() /
        (std::string("cfdx_vmfl036_")+tag+"_"+l.name+".axmesh");
    generate(path,l);
    const auto mesh=AxisymMesh::read(path.string());
    if(mesh.outer_radius < 49.999 || mesh.nr != static_cast<std::size_t>(l.nr) ||
       mesh.nt != static_cast<std::size_t>(l.nt))
        throw std::runtime_error("VMFL036 mesh geometry contract failed");

    AxisymControls c;
    c.max_outer_iterations=700;
    c.momentum_sweeps=25;
    c.pressure_sweeps=100;
    c.momentum_tolerance=2e-6;
    c.pressure_tolerance=2e-7;
    c.continuity_tolerance=2e-6;
    c.alpha_u=0.65;
    c.alpha_p=0.25;
    VMFL036AxisymmetricSolver solver(mesh,1.0,mu);
    const auto result=solver.solve(c);
    std::error_code ec; std::filesystem::remove(path,ec);
    if(!result.converged) throw std::runtime_error(std::string("VMFL036 axisymmetric solver did not converge on ")+l.name);
    if(!std::isfinite(result.cd_total) || !std::isfinite(result.drag_total))
        throw std::runtime_error("VMFL036 drag is non-finite");
    std::cout << "VMFL036_AXISYM RESULT level=" << l.name
              << " iterations=" << result.iterations
              << " continuity=" << result.continuity
              << " momentum=" << result.momentum_residual
              << " pcorr=" << result.pressure_correction
              << " Cd_pressure=" << result.cd_pressure
              << " Cd_viscous=" << result.cd_viscous
              << " Cd_total=" << result.cd_total
              << " Cd_reference=" << cd_ref
              << " Cd_relative_error=" << std::abs(result.cd_total-cd_ref)/cd_ref << "\n";
    return result;
}
}

int main(int argc,char** argv)
{
    try {
        const bool quick=argc>1 && std::string(argv[1])=="--quick";
        const std::vector<Level> levels = quick
            ? std::vector<Level>{{24,36,"coarse"}}
            : std::vector<Level>{{24,36,"coarse"},{36,54,"medium"},{48,72,"fine"}};

        std::vector<AxisymResult> r;
        if(!quick) {
            const auto fluent_exact=run_level(Level{36,54,"medium"},0.02,1.0875,"re50");
            const double re50_err=std::abs(fluent_exact.cd_total-1.0875)/1.0875;
            std::cout << "VMFL036_AXISYM FLUENT_EXACT Re=50 Cd=" << fluent_exact.cd_total
                      << " reference=1.0875 relative_error=" << re50_err << "\n";
            if(re50_err>0.10) throw std::runtime_error("VMFL036 Fluent-exact Re50 comparison exceeds 10%");
        }
        for(const auto& l:levels) r.push_back(run_level(l,0.01,1.0895,"re100"));

        if(!quick){
            const double ref=1.0895;
            const double fine_err=std::abs(r.back().cd_total-ref)/ref;
            const double refine=std::abs(r.back().cd_total-r[r.size()-2].cd_total) /
                std::max(1e-12,std::abs(r.back().cd_total));
            std::cout << "VMFL036_AXISYM REFINEMENT fine_relative_error=" << fine_err
                      << " fine_to_medium_change=" << refine << "\n";
            if(fine_err>0.10) throw std::runtime_error("VMFL036 fine-grid Cd error exceeds 10%");
            if(refine>0.08) throw std::runtime_error("VMFL036 fine/medium Cd change exceeds 8%");
        }
        std::cout << "VMFL036_AXISYMMETRIC_VALIDATION: PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "VMFL036_AXISYMMETRIC_VALIDATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}
