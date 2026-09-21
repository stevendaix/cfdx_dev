#include "convection_assembly.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::core::numerics;

namespace {
struct FaceData {
    Vec3 Sf;
    Vec3 centre;
};

std::vector<FaceData> geometry(const Mesh& mesh)
{
    std::vector<FaceData> g(mesh.n_faces());
    for(std::size_t f=0; f<mesh.n_faces(); ++f) {
        const auto off=mesh.faces().offsets_data()[f];
        const auto n=mesh.faces().offsets_data()[f+1]-off;
        const auto fg=compute_face_geometry(
            mesh.points().x_data(),mesh.points().y_data(),mesh.points().z_data(),
            mesh.faces().vertices_data(),off,n);
        g[f]={fg.Sf,fg.centre};
    }
    return g;
}

double interpolation_weight(const std::string& scheme, double F)
{
    if(scheme=="upwind") return F>=0.0 ? 1.0 : 0.0;
    if(scheme=="linear") return 0.5;
    if(scheme=="limited") return F>=0.0 ? 0.75 : 0.25;
    throw std::invalid_argument("assembleConvectionCSR: unsupported scheme: "+scheme);
}
}

bool assembleConvectionCSR(const Mesh& mesh,
                           const std::vector<ScalarCellField>& velocity,
                           const std::string& scheme,
                           int component_idx,
                           SparseMatrix& A)
{
    if(velocity.size()!=3 || component_idx<0 || component_idx>=3)
        return false;
    for(const auto& v:velocity)
        if(v.size()!=mesh.n_cells() || v.dimension()!=1) return false;

    const auto g=geometry(mesh);
    const std::size_t nc=mesh.n_cells();
    std::vector<std::map<std::size_t,double>> rows(nc);

    for(std::size_t f=0; f<mesh.n_faces(); ++f) {
        const std::size_t o=mesh.ownership().owner(f);
        if(o>=nc) return false;
        Vec3 Sf=g[f].Sf;
        if(mesh.ownership().neighbour(f)>=0) {
            const std::size_t n=static_cast<std::size_t>(mesh.ownership().neighbour(f));
            const Vec3 d=mesh.cells().n_cells()>n
                ? Vec3(g[f].centre.x,g[f].centre.y,g[f].centre.z)-Vec3(0.0)
                : Vec3(0.0);
            const auto co=compute_cell_geometry(
                nullptr,nullptr,nullptr,nullptr,0);
            (void)d; (void)co;
            // Re-orient the face area vector from owner to neighbour.
            const auto oo=mesh.cells().offsets_data()[o];
            const auto on=mesh.cells().offsets_data()[o+1]-oo;
            const auto no=mesh.cells().offsets_data()[n];
            const auto nn=mesh.cells().offsets_data()[n+1]-no;
            std::vector<Vec3> fcentres(mesh.n_faces());
            for(std::size_t ff=0;ff<mesh.n_faces();++ff) fcentres[ff]=g[ff].centre;
            const auto co_geom=compute_cell_geometry(
                fcentres.data(), nullptr, mesh.cells().faces_data()+oo, on);
            const auto cn_geom=compute_cell_geometry(
                fcentres.data(), nullptr, mesh.cells().faces_data()+no, nn);
            const Vec3 delta=cn_geom.centre-co_geom.centre;
            if(Sf.dot(delta)<0.0) Sf=Sf*(-1.0);

            const double F=0.5*(velocity[0](o)*Sf.x+velocity[1](o)*Sf.y+velocity[2](o)*Sf.z+
                                velocity[0](n)*Sf.x+velocity[1](n)*Sf.y+velocity[2](n)*Sf.z);
            const double w=interpolation_weight(scheme,F);
            const double a_o=F*w;
            const double a_n=F*(1.0-w);
            rows[o][o]+=a_o; rows[o][n]+=a_n;
            rows[n][o]-=a_o; rows[n][n]-=a_n;
        } else {
            const double F=velocity[0](o)*Sf.x+velocity[1](o)*Sf.y+velocity[2](o)*Sf.z;
            rows[o][o]+=F;
        }
    }

    A=SparseMatrix(nc,nc);
    for(std::size_t r=0;r<nc;++r)
        for(const auto& [col,val]:rows[r])
            if(std::abs(val)>0.0) A.push_back(r,col,val);
    A.finalize();
    return true;
}

bool assembleMomentumCSR(const Mesh& mesh,
                         const std::vector<ScalarCellField>& velocity,
                         const ScalarCellField& pressure,
                         const std::string& convection_scheme,
                         SparseMatrix& A,
                         Vector& b,
                         double diffusion_coefficient)
{
    if(velocity.size()!=3 || pressure.size()!=mesh.n_cells() || pressure.dimension()!=1 ||
       diffusion_coefficient<0.0 || !std::isfinite(diffusion_coefficient))
        return false;
    if(!assembleConvectionCSR(mesh,velocity,convection_scheme,0,A))
        return false;

    const auto g=geometry(mesh);
    const std::size_t nc=mesh.n_cells();
    std::vector<std::map<std::size_t,double>> rows(nc);
    for(std::size_t r=0;r<A.n_rows();++r)
        for(std::size_t k=A.row_offsets_data()[r];k<A.row_offsets_data()[r+1];++k)
            rows[r][A.columns_data()[k]]+=A.values_data()[k];

    b=Vector(nc,0.0);
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const std::size_t o=mesh.ownership().owner(f);
        Vec3 Sf=g[f].Sf;
        double pf=pressure(o);
        if(mesh.ownership().neighbour(f)>=0) {
            const std::size_t n=static_cast<std::size_t>(mesh.ownership().neighbour(f));
            pf=0.5*(pressure(o)+pressure(n));
            const auto oo=mesh.cells().offsets_data()[o];
            const auto on=mesh.cells().offsets_data()[o+1]-oo;
            const auto no=mesh.cells().offsets_data()[n];
            const auto nn=mesh.cells().offsets_data()[n+1]-no;
            std::vector<Vec3> fcentres(mesh.n_faces());
            for(std::size_t ff=0;ff<mesh.n_faces();++ff) fcentres[ff]=g[ff].centre;
            const auto co=compute_cell_geometry(fcentres.data(),nullptr,mesh.cells().faces_data()+oo,on);
            const auto cn=compute_cell_geometry(fcentres.data(),nullptr,mesh.cells().faces_data()+no,nn);
            if(Sf.dot(cn.centre-co.centre)<0.0) Sf=Sf*(-1.0);
            b(o)-=pf*Sf.x; b(n)+=pf*Sf.x;
        } else {
            b(o)-=pf*Sf.x;
        }
        if(diffusion_coefficient>0.0 && mesh.ownership().neighbour(f)>=0) {
            const std::size_t n=static_cast<std::size_t>(mesh.ownership().neighbour(f));
            const double d=(g[f].centre-g[f].centre).mag();
            (void)d;
            const auto oo=mesh.cells().offsets_data()[o];
            const auto on=mesh.cells().offsets_data()[o+1]-oo;
            const auto no=mesh.cells().offsets_data()[n];
            const auto nn=mesh.cells().offsets_data()[n+1]-no;
            std::vector<Vec3> fcentres(mesh.n_faces());
            for(std::size_t ff=0;ff<mesh.n_faces();++ff) fcentres[ff]=g[ff].centre;
            const auto co=compute_cell_geometry(fcentres.data(),nullptr,mesh.cells().faces_data()+oo,on);
            const auto cn=compute_cell_geometry(fcentres.data(),nullptr,mesh.cells().faces_data()+no,nn);
            const double dist=(cn.centre-co.centre).mag();
            const double D=diffusion_coefficient*g[f].Sf.mag()/std::max(dist,1e-14);
            rows[o][o]+=D; rows[o][n]-=D;
            rows[n][n]+=D; rows[n][o]-=D;
        }
    }
    A=SparseMatrix(nc,nc);
    for(std::size_t r=0;r<nc;++r)
        for(const auto& [col,val]:rows[r]) A.push_back(r,col,val);
    A.finalize();
    return true;
}
