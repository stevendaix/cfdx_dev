#include "convection_assembly.h"
#include "cfdx/core/geometry/face_geometry.h"
#include <algorithm>
#include <cmath>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::core::numerics;

namespace {

struct GeometryData {
    std::vector<Vec3> face_centres;
    std::vector<Vec3> face_area_vectors;
    std::vector<Vec3> cell_centres;
};

GeometryData build_geometry(const Mesh& mesh)
{
    GeometryData g;
    g.face_centres.resize(mesh.n_faces());
    g.face_area_vectors.resize(mesh.n_faces());
    g.cell_centres.resize(mesh.n_cells());

    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        const auto off=mesh.faces().offsets_data()[f];
        const auto n=mesh.faces().offsets_data()[f+1]-off;
        const auto fg=compute_face_geometry(
            mesh.points().x_data(),mesh.points().y_data(),mesh.points().z_data(),
            mesh.faces().vertices_data(),off,n);
        g.face_centres[f]=fg.centre;
        g.face_area_vectors[f]=fg.Sf;
    }

    for (std::size_t c=0; c<mesh.n_cells(); ++c) {
        const auto off=mesh.cells().offsets_data()[c];
        const auto n=mesh.cells().offsets_data()[c+1]-off;
        Vec3 centre;
        double weight=0.0;
        for (std::size_t k=0;k<n;++k) {
            const auto f=mesh.cells().faces_data()[off+k];
            const double a=g.face_area_vectors[f].mag();
            centre=centre+g.face_centres[f]*a;
            weight+=a;
        }
        if (!(weight>0.0)) throw std::runtime_error("convection assembly: invalid cell geometry");
        g.cell_centres[c]=centre*(1.0/weight);
    }
    return g;
}

double scheme_weight(const std::string& scheme,double F)
{
    if(scheme=="upwind") return F>=0.0 ? 1.0 : 0.0;
    if(scheme=="linear") return 0.5;
    if(scheme=="limited") {
        // Bounded blend between upwind and central interpolation.
        return F>=0.0 ? 0.75 : 0.25;
    }
    throw std::invalid_argument("assembleConvectionCSR: unsupported scheme: "+scheme);
}

double face_flux(const std::vector<ScalarCellField>& velocity,
                 std::size_t o,std::size_t n,const Vec3& Sf)
{
    const double uo=velocity[0](o)*Sf.x+velocity[1](o)*Sf.y+velocity[2](o)*Sf.z;
    if(n>=velocity[0].size()) return uo;
    const double un=velocity[0](n)*Sf.x+velocity[1](n)*Sf.y+velocity[2](n)*Sf.z;
    return 0.5*(uo+un);
}

}

bool assembleConvectionCSR(const Mesh& mesh,
                           const std::vector<ScalarCellField>& velocity,
                           const std::string& scheme,
                           int component_idx,
                           SparseMatrix& A)
{
    if(velocity.size()!=3 || component_idx<0 || component_idx>=3) return false;
    for(const auto& v:velocity)
        if(v.size()!=mesh.n_cells() || v.dimension()!=1) return false;

    const auto g=build_geometry(mesh);
    const std::size_t nc=mesh.n_cells();
    std::vector<std::map<std::size_t,double>> rows(nc);

    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const std::size_t o=mesh.ownership().owner(f);
        if(o>=nc) return false;
        const int nraw=mesh.ownership().neighbour(f);
        Vec3 Sf=g.face_area_vectors[f];

        if(nraw>=0) {
            const std::size_t n=static_cast<std::size_t>(nraw);
            if(n>=nc) return false;
            const Vec3 delta=g.cell_centres[n]-g.cell_centres[o];
            if(Sf.dot(delta)<0.0) Sf=Sf*(-1.0);
            const double F=face_flux(velocity,o,n,Sf);
            const double w=scheme_weight(scheme,F);
            const double a_o=F*w;
            const double a_n=F*(1.0-w);

            rows[o][o]+=a_o;
            rows[o][n]+=a_n;
            rows[n][o]-=a_o;
            rows[n][n]-=a_n;
        } else {
            const double F=velocity[0](o)*Sf.x+
                           velocity[1](o)*Sf.y+
                           velocity[2](o)*Sf.z;
            rows[o][o]+=F;
        }
    }

    A=SparseMatrix(nc,nc);
    for(std::size_t r=0;r<nc;++r)
        for(const auto& [col,val]:rows[r])
            if(std::isfinite(val) && std::abs(val)>0.0) A.push_back(r,col,val);
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
    if(velocity.size()!=3 || pressure.size()!=mesh.n_cells() ||
       pressure.dimension()!=1 || diffusion_coefficient<0.0 ||
       !std::isfinite(diffusion_coefficient))
        return false;

    if(!assembleConvectionCSR(mesh,velocity,convection_scheme,0,A))
        return false;

    const auto g=build_geometry(mesh);
    const std::size_t nc=mesh.n_cells();
    std::vector<std::map<std::size_t,double>> rows(nc);
    for(std::size_t r=0;r<nc;++r)
        for(std::size_t k=A.row_offsets_data()[r];k<A.row_offsets_data()[r+1];++k)
            rows[r][A.columns_data()[k]]+=A.values_data()[k];

    b=Vector(nc,0.0);
    for(std::size_t f=0;f<mesh.n_faces();++f) {
        const std::size_t o=mesh.ownership().owner(f);
        Vec3 Sf=g.face_area_vectors[f];
        const int nraw=mesh.ownership().neighbour(f);

        if(nraw>=0) {
            const std::size_t n=static_cast<std::size_t>(nraw);
            if(n>=nc) return false;
            const Vec3 delta=g.cell_centres[n]-g.cell_centres[o];
            if(Sf.dot(delta)<0.0) Sf=Sf*(-1.0);
            const double pf=0.5*(pressure(o)+pressure(n));
            b(o)-=pf*Sf.x;
            b(n)+=pf*Sf.x;

            if(diffusion_coefficient>0.0) {
                const double dist=delta.mag();
                const double D=diffusion_coefficient*Sf.mag()/std::max(dist,1e-14);
                rows[o][o]+=D; rows[o][n]-=D;
                rows[n][n]+=D; rows[n][o]-=D;
            }
        } else {
            // The legacy API has no velocity boundary values; use the
            // owner value as a zero-gradient extrapolation.
            b(o)-=pressure(o)*Sf.x;
        }
    }

    A=SparseMatrix(nc,nc);
    for(std::size_t r=0;r<nc;++r)
        for(const auto& [col,val]:rows[r])
            if(std::isfinite(val)) A.push_back(r,col,val);
    A.finalize();
    return true;
}
