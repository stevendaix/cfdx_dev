// Independent NACA trailing-edge gradient verification, not a flow validation.
// Affine manufactured field supplies its exact values at every physical face.
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/core/fvm/least_squares_gradient.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>
using namespace cfdx::core;
using namespace cfdx::physics;
static double solution(const Vec3& x) { return 2.0*x.x-3.0*x.y+0.13; }
int main(int argc,char** argv) {
 try {
    if(argc!=4) throw std::invalid_argument("expected NACA n64 n128 n256 mesh paths");
    for(int m=1;m<4;++m) {
        Mesh mesh;
        if(!cfdx::io::read_mesh_hdf5(argv[m],mesh) || !mesh.topo_validate().ok)
            throw std::runtime_error("invalid or unreadable NACA mesh");
        const auto g=build_fv_geometry(mesh);
        double max_error=0.0, max_te_error=0.0;
        std::size_t te_cells=0, extended_cells=0;
        const auto& own=mesh.ownership();
        for(std::size_t c=0;c<mesh.n_cells();++c) {
            std::vector<Vec3> coords;
            std::vector<double> values;
            const auto off=mesh.cells().offsets_data()[c];
            const auto end=mesh.cells().offsets_data()[c+1];
            for(auto k=off;k<end;++k) {
                const auto f=mesh.cells().faces_data()[k];
                const auto nb=own.neighbour(f);
                if(nb>=0) {
                    const auto other=own.owner(f)==c
                        ? static_cast<std::size_t>(nb) : own.owner(f);
                    const auto pt=g.cell_centres[other];
                    coords.push_back(pt);
                    values.push_back(solution(pt));
                } else {
                    const auto patch=g.face_patch[f];
                    if(patch<mesh.boundary().n_patches() &&
                       (mesh.boundary().patch(patch).type==PatchType::EMPTY ||
                        mesh.boundary().patch(patch).name=="front" ||
                        mesh.boundary().patch(patch).name=="back"))
                        continue;
                    const auto pt=g.face_centres[f];
                    coords.push_back(pt);
                    values.push_back(solution(pt));
                }
            }
            StencilQuality q;
            Vec3 grad;
            try {
                grad=weighted_least_squares_gradient_robust(
                    g.cell_centres[c],solution(g.cell_centres[c]),coords,values,
                    GradientWeighting::INVERSE_DISTANCE_SQUARED,2,&q);
            } catch (const std::runtime_error&) {
                // Stretched N128/N256 topology can leave the immediate
                // neighbour stencil rank deficient. Add the second ring,
                // retaining exact MMS boundary values and rank checks.
                ++extended_cells;
                for(auto k=off;k<end;++k) {
                    const auto f=mesh.cells().faces_data()[k];
                    const auto nb=own.neighbour(f);
                    if(nb<0) continue;
                    const auto other=own.owner(f)==c
                        ? static_cast<std::size_t>(nb) : own.owner(f);
                    for(auto k2=mesh.cells().offsets_data()[other];
                        k2<mesh.cells().offsets_data()[other+1]; ++k2) {
                        const auto f2=mesh.cells().faces_data()[k2];
                        const auto nr2=own.neighbour(f2);
                        Vec3 pt;
                        if(nr2>=0) {
                            const auto other2=own.owner(f2)==other
                                ? static_cast<std::size_t>(nr2) : own.owner(f2);
                            if(other2==c) continue;
                            pt=g.cell_centres[other2];
                        } else {
                            const auto patch=g.face_patch[f2];
                            if(patch<mesh.boundary().n_patches() &&
                               (mesh.boundary().patch(patch).type==PatchType::EMPTY ||
                                mesh.boundary().patch(patch).name=="front" ||
                                mesh.boundary().patch(patch).name=="back"))
                                continue;
                            pt=g.face_centres[f2];
                        }
                        coords.push_back(pt);
                        values.push_back(solution(pt));
                    }
                }
                try {
                    grad=weighted_least_squares_gradient_robust(
                        g.cell_centres[c],solution(g.cell_centres[c]),coords,values,
                        GradientWeighting::INVERSE_DISTANCE_SQUARED,2,&q);
                } catch (const std::runtime_error& exc) {
                    const auto x=g.cell_centres[c];
                    throw std::runtime_error("WLS 2-ring cell="+std::to_string(c)+
                        " xyz="+std::to_string(x.x)+":"+std::to_string(x.y)+
                        " samples="+std::to_string(coords.size())+
                        " rank="+std::to_string(q.rank)+
                        " scale="+std::to_string(q.scale)+
                        " zmin="+std::to_string(std::min_element(coords.begin(),coords.end(),
                             [](const Vec3& a,const Vec3& b){return a.z<b.z;})->z)+
                        " zmax="+std::to_string(std::max_element(coords.begin(),coords.end(),
                             [](const Vec3& a,const Vec3& b){return a.z<b.z;})->z)+
                        " "+exc.what());
                }
            }
            const double error=std::sqrt((grad.x-2)*(grad.x-2)+
                                         (grad.y+3)*(grad.y+3)+grad.z*grad.z);
            max_error=std::max(max_error,error);
            const auto x=g.cell_centres[c];
            if(x.x>=0.90 && x.x<=1.1 && std::abs(x.y)<=0.1) {
                ++te_cells;
                max_te_error=std::max(max_te_error,error);
            }
        }
        std::cout<<"NACA_TE_GRADIENT_MMS level="<< (64<<(m-1))
                 <<" cells="<<mesh.n_cells()<<" te_cells="<<te_cells
                 <<" extended_cells="<<extended_cells<<" global_linf="<<max_error<<" te_linf="<<max_te_error<<'\n';
        if(te_cells==0 || !(max_te_error<1e-8 && max_error<1e-8))
            throw std::runtime_error("trailing-edge affine gradient gate failed");
    }
    std::cout<<"NACA_TE_GRADIENT_MMS: PASS (manufactured gradients only)\n";
 } catch(const std::exception& e) {
    std::cerr<<"NACA_TE_GRADIENT_MMS: FAIL "<<e.what()<<'\n';
    return 1;
 }
 return 0;
}
