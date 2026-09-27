#include "cfdx/physics/axisymmetric_vmfl036_solver.h"
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <limits>

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

void audit_mesh(const AxisymMesh& mesh, const Level& l)
{
    if(mesh.nr != static_cast<std::size_t>(l.nr) || mesh.nt != static_cast<std::size_t>(l.nt))
        throw std::runtime_error("VMFL036 mesh resolution mismatch");
    if(!std::isfinite(mesh.diameter) || !std::isfinite(mesh.outer_radius) ||
       !(mesh.diameter>0.0) || !(mesh.outer_radius>mesh.diameter))
        throw std::runtime_error("VMFL036 invalid mesh dimensions");

    std::size_t sphere=0, outer=0, axis=0, internal=0;
    for(std::size_t i=0;i<mesh.nodes.size();++i){
        const auto& p=mesh.nodes[i];
        if(!std::isfinite(p.x) || !std::isfinite(p.r) || p.r<0.0)
            throw std::runtime_error("VMFL036 invalid node geometry at node "+std::to_string(i));
    }
    for(std::size_t i=0;i<mesh.cells.size();++i){
        const auto& c=mesh.cells[i];
        if(!std::isfinite(c.area) || !std::isfinite(c.volume) ||
           !(c.area>0.0) || !(c.volume>0.0))
            throw std::runtime_error("VMFL036 invalid cell metric at cell "+std::to_string(i));
    }
    for(std::size_t i=0;i<mesh.faces.size();++i){
        const auto& f=mesh.faces[i];
        if(!std::isfinite(f.ds) || !std::isfinite(f.area) || !(f.ds>0.0) || f.area<0.0)
            throw std::runtime_error("VMFL036 invalid face metric at face "+std::to_string(i));
        if(f.neighbour==static_cast<std::size_t>(-1)){
            if(f.sphere) ++sphere;
            if(f.outer) ++outer;
            if(f.axis) ++axis;
            if((f.sphere?1:0)+(f.outer?1:0)+(f.axis?1:0)!=1)
                throw std::runtime_error("VMFL036 boundary face classification is not exclusive");
            const auto& p0=mesh.nodes[f.n0];
            const auto& p1=mesh.nodes[f.n1];
            if(f.axis && (std::abs(p0.r)>1e-14 || std::abs(p1.r)>1e-14))
                throw std::runtime_error("VMFL036 axis face is not geometrically on r=0");
            if(f.axis && std::abs(f.area)>1e-14)
                throw std::runtime_error("VMFL036 axis face has non-zero revolution area");
            if(!f.axis && !(f.area>0.0))
                throw std::runtime_error("VMFL036 non-axis boundary face has zero area");
        } else {
            ++internal;
            if(f.neighbour>=mesh.cells.size() || f.neighbour==f.owner)
                throw std::runtime_error("VMFL036 invalid internal face connectivity");
            const auto& o=mesh.cells[f.owner];
            const auto& n=mesh.cells[f.neighbour];
            const double orient=f.nx*(n.cx-o.cx)+f.nr*(n.cr-o.cr);
            if(!(orient>0.0))
                throw std::runtime_error("VMFL036 internal face normal is not owner-to-neighbour");
        }
    }
    if(sphere!=mesh.nt || outer!=mesh.nt || axis!=2*mesh.nr)
        throw std::runtime_error("VMFL036 boundary topology count mismatch");
    if(internal+ sphere+outer+axis != mesh.faces.size())
        throw std::runtime_error("VMFL036 face accounting mismatch");
}

AxisymResult run_level(const Level& l, double mu, double cd_ref, const char* tag)
{
    const auto path=std::filesystem::temp_directory_path() /
        (std::string("cfdx_vmfl036_")+tag+"_"+l.name+".axmesh");
    generate(path,l);
    const auto mesh=AxisymMesh::read(path.string());
    if(mesh.outer_radius < 49.999)
        throw std::runtime_error("VMFL036 mesh geometry contract failed");
    audit_mesh(mesh,l);

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
    AxisymResult result;
    try {
        result=solver.solve(c);
    } catch (...) {
        std::error_code ec; std::filesystem::remove(path,ec);
        throw;
    }
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
            if(re50_err>0.03) throw std::runtime_error("VMFL036 Fluent-exact Re50 comparison exceeds 3%");
        }
        for(const auto& l:levels) r.push_back(run_level(l,0.01,1.0895,"re100"));

        if(!quick){
            const double ref=1.0895;
            const double fine_err=std::abs(r.back().cd_total-ref)/ref;
            const double refine=std::abs(r.back().cd_total-r[r.size()-2].cd_total) /
                std::max(1e-12,std::abs(r.back().cd_total));
            std::cout << "VMFL036_AXISYM REFINEMENT fine_relative_error=" << fine_err
                      << " fine_to_medium_change=" << refine << "\n";
            if(fine_err>0.03) throw std::runtime_error("VMFL036 fine-grid Cd error exceeds 3%");
            if(refine>0.02) throw std::runtime_error("VMFL036 fine/medium Cd change exceeds 2%");
        }
        std::cout << "VMFL036_AXISYMMETRIC_VALIDATION: PASS\n";
        return 0;
    } catch(const std::exception& e) {
        std::cerr << "VMFL036_AXISYMMETRIC_VALIDATION: FAIL: " << e.what() << "\n";
        return 1;
    }
}
