#include "cfdx/core/numerics/interpolation.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/mesh/mesh.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

using namespace cfdx::core;

namespace {

struct Grid { Mesh mesh; GeometryCache geometry; std::vector<std::size_t> internal_faces; };

Grid make_hex_frustum_chain(std::size_t n)
{
    if (n<2) throw std::invalid_argument("n must be >= 2");
    Mesh m;
    m.points().resize(6*(n+1));
    const double pi=std::acos(-1.0);
    auto pid=[](std::size_t i,std::size_t k){return 6*i+k;};
    for (std::size_t i=0;i<=n;++i) {
        const double x=static_cast<double>(i)/n;
        const double r=0.22+0.08*x; // non-affine cross-section: frustum family
        for (std::size_t k=0;k<6;++k) {
            const double a=2.0*pi*static_cast<double>(k)/6.0;
            m.points().set(pid(i,k),x,0.5+r*std::cos(a),0.5+r*std::sin(a));
        }
    }

    std::vector<std::size_t> left(n),right(n);
    std::vector<std::array<std::size_t,6>> side(n);
    auto face=[&](const std::vector<std::size_t>& v){
        const auto id=m.n_faces(); m.faces().push_face(v); return id;
    };

    for (std::size_t i=0;i<n;++i) {
        left[i]=face({pid(i,0),pid(i,1),pid(i,2),pid(i,3),pid(i,4),pid(i,5)});
        right[i]=face({pid(i+1,5),pid(i+1,4),pid(i+1,3),pid(i+1,2),pid(i+1,1),pid(i+1,0)});
        for (std::size_t k=0;k<6;++k) {
            const std::size_t kn=(k+1)%6;
            side[i][k]=face({pid(i,k),pid(i+1,k),pid(i+1,kn),pid(i,kn)});
        }
    }

    m.ownership().resize(m.n_faces());
    std::vector<std::vector<std::size_t>> cf(n);
    for (std::size_t i=0;i<n;++i) {
        if (i+1<n) {
            left[i+1]=right[i];
            m.ownership().set_owner(right[i],i);
            m.ownership().set_neighbour(right[i],static_cast<std::int64_t>(i+1));
        }
        if (i==0) m.ownership().set_owner(left[i],i);
        for (const auto f:side[i]) {
            m.ownership().set_owner(f,i);
            m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
        }
        cf[i]={left[i],right[i],side[i][0],side[i][1],side[i][2],
               side[i][3],side[i][4],side[i][5]};
    }
    m.ownership().set_neighbour(left[0],FaceOwnership::BOUNDARY);
    m.ownership().set_neighbour(right[n-1],FaceOwnership::BOUNDARY);
    for (const auto& faces:cf) m.cells().push_cell(faces);

    Patch inlet; inlet.name="inlet"; inlet.type=PatchType::WALL; inlet.face_ids={left[0]};
    Patch outlet; outlet.name="outlet"; outlet.type=PatchType::WALL; outlet.face_ids={right[n-1]};
    Patch walls; walls.name="walls"; walls.type=PatchType::WALL;
    for (std::size_t i=0;i<n;++i) for (const auto f:side[i]) walls.face_ids.push_back(f);
    m.boundary().add_patch(inlet); m.boundary().add_patch(outlet); m.boundary().add_patch(walls);

    const auto geometry=make_geometry_cache(m);
    Grid g{std::move(m),geometry,{}};
    for (std::size_t f=0;f<g.mesh.n_faces();++f)
        if (g.mesh.ownership().neighbour(f)>=0) g.internal_faces.push_back(f);
    return g;
}

Field<double,Location::FACE> axial_flux(const Grid& g)
{
    Field<double,Location::FACE> f(g.mesh.n_faces(),"F","m3/s",1);
    for (std::size_t i=0;i<g.mesh.n_faces();++i) f(i)=g.geometry.face_Sf[i].x;
    return f;
}

double smooth(const Vec3& p)
{
    const double pi=std::acos(-1.0);
    return std::sin(pi*p.x)+0.15*std::cos(2.0*pi*p.x);
}

double order(double a,double b){return std::log(a/b)/std::log(2.0);}

void conservation()
{
    const auto g=make_hex_frustum_chain(32);
    const auto flux=axial_flux(g);
    Field<double,Location::CELL> phi(g.mesh.n_cells(),"phi","1",1);
    phi.fill(3.0);
    for (const auto scheme:{InterpScheme::UPWIND,InterpScheme::LINEAR,
                            InterpScheme::BLENDED,InterpScheme::QUICK,
                            InterpScheme::QUICK_BOUNDED}) {
        const auto face=interpolate_cell_to_face(phi,g.mesh,g.geometry,scheme,&flux,
                                                  LimiterType::NONE,nullptr,1.0);
        double integral=0.0;
        for (std::size_t c=0;c<g.mesh.n_cells();++c) {
            double div=0.0;
            const auto off=g.mesh.cells().offsets_data()[c];
            const auto end=g.mesh.cells().offsets_data()[c+1];
            for (auto k=off;k<end;++k) {
                const auto f=g.mesh.cells().faces_data()[k];
                const double s=g.mesh.ownership().owner(f)==c ? 1.0 : -1.0;
                div+=s*flux(f)*face(f);
            }
            integral+=g.geometry.cell_volumes[c]*div/g.geometry.cell_volumes[c];
        }
        std::cout<<"N4_POLY_CONSERVATION scheme="<<to_string(scheme)
                 <<" integral="<<integral<<"\n";
        if (std::abs(integral)>1e-11) throw std::runtime_error("polyhedral conservation failure");
    }
}

void smooth_order()
{
    const std::vector<std::size_t> ns={8u,16u,32u};
    for (const auto scheme:{InterpScheme::UPWIND,InterpScheme::LINEAR,InterpScheme::QUICK}) {
        std::vector<double> e;
        for (const auto n:ns) {
            const auto g=make_hex_frustum_chain(n);
            const auto flux=axial_flux(g);
            Field<double,Location::CELL> phi(g.mesh.n_cells(),"phi","1",1);
            for (std::size_t c=0;c<g.mesh.n_cells();++c) phi(c)=smooth(g.geometry.cell_centres[c]);
            const auto face=interpolate_cell_to_face(phi,g.mesh,g.geometry,scheme,&flux);
            double e2=0.0, w=0.0;
            for (const auto f:g.internal_faces) {
                const double d=face(f)-smooth(g.geometry.face_centres[f]);
                e2+=d*d; w+=1.0;
            }
            e.push_back(std::sqrt(e2/w));
        }
        const double p=order(e[1],e[2]);
        std::cout<<"N4_POLY_ORDER scheme="<<to_string(scheme)
                 <<" n=8 "<<e[0]<<" n=16 "<<e[1]<<" n=32 "<<e[2]
                 <<" order="<<p<<"\n";
        if (scheme==InterpScheme::UPWIND) {
            if (!(p>0.70)) throw std::runtime_error("polyhedral upwind lost first order");
        } else if (!(p>1.50)) {
            throw std::runtime_error("polyhedral second-order reconstruction gate failed");
        }
    }
}

void bounded_step()
{
    const auto g=make_hex_frustum_chain(32);
    const auto flux=axial_flux(g);
    Field<double,Location::CELL> phi(g.mesh.n_cells(),"phi","1",1);
    for (std::size_t c=0;c<g.mesh.n_cells();++c)
        phi(c)=g.geometry.cell_centres[c].x<0.5 ? 0.0 : 1.0;
    const auto grad=compute_gradient_least_squares_quadratic(phi,g.mesh);
    for (const auto lim:{LimiterType::MINMOD,LimiterType::VANLEER,
                         LimiterType::SUPERBEE,LimiterType::VAN_ALBADA,LimiterType::MC}) {
        const auto face=interpolate_cell_to_face(phi,g.mesh,g.geometry,InterpScheme::LIMITED,
                                                 &flux,lim,&grad);
        double worst=0.0;
        for (const auto f:g.internal_faces) {
            const auto o=g.mesh.ownership().owner(f);
            const auto n=static_cast<std::size_t>(g.mesh.ownership().neighbour(f));
            worst=std::max(worst,std::max(std::min(phi(o),phi(n))-face(f),
                                          face(f)-std::max(phi(o),phi(n))));
        }
        std::cout<<"N4_POLY_BOUNDED limiter="<<to_string(lim)<<" worst="<<worst<<"\n";
        if (worst>1e-12) throw std::runtime_error("polyhedral TVD boundedness failure");
    }
}

}

int main()
{
    try {
        conservation();
        smooth_order();
        bounded_step();
        std::cout<<"CONVECTION_POLYHEDRAL_QUALIFICATION: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr<<"CONVECTION_POLYHEDRAL_QUALIFICATION: FAIL: "<<e.what()<<"\n";
        return 1;
    }
}
