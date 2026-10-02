// Issue #461 — N4 complete convection qualification on skewed tetrahedral meshes.
//
// Qualification scope:
//   * conservative face-flux assembly on a genuinely 3-D tetrahedral mesh;
//   * affine skew/non-orthogonal tetrahedral refinement;
//   * constant-field conservation;
//   * linear exactness of the quadratic-upwind reconstruction;
//   * measured smooth-field order for upwind and QUICK;
//   * boundedness of all five TVD limiters and bounded QUICK;
//   * production compute_convection diagnostics.
//
// The campaign intentionally does not call central/BLENDED "second order" on
// tetrahedra: their two-point face interpolation is not generally
// face-centroid exact on non-orthogonal polyhedra. That is a documented
// discretisation property, not a reason to weaken the gates.

#include "cfdx/core/numerics/convection.h"
#include "cfdx/core/numerics/interpolation.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include "common/test_harness.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {

struct Grid {
    Mesh mesh;
    GeometryCache geometry;
    std::vector<std::size_t> interior;
};

static const std::array<std::array<unsigned, 3>, 8> kCorner = {{
    {{0,0,0}}, {{1,0,0}}, {{1,1,0}}, {{0,1,0}},
    {{0,0,1}}, {{1,0,1}}, {{1,1,1}}, {{0,1,1}}
}};

static const std::array<std::array<unsigned, 4>, 6> kTets = {{
    {{0,1,2,6}}, {{0,1,5,6}}, {{0,3,2,6}},
    {{0,3,7,6}}, {{0,4,5,6}}, {{0,4,7,6}}
}};

std::string face_key(std::size_t a, std::size_t b, std::size_t c)
{
    std::array<std::size_t,3> v{a,b,c};
    std::sort(v.begin(), v.end());
    return std::to_string(v[0]) + "/" + std::to_string(v[1]) + "/" + std::to_string(v[2]);
}

Grid make_tet_grid(std::size_t n, double shear_xy, double shear_xz)
{
    const std::size_t nv = n + 1;
    const auto vid = [nv](std::size_t x, std::size_t y, std::size_t z) {
        return x + nv * (y + nv * z);
    };

    Mesh m;
    m.points().resize(nv * nv * nv);
    for (std::size_t z=0; z<=n; ++z)
        for (std::size_t y=0; y<=n; ++y)
            for (std::size_t x=0; x<=n; ++x) {
                const double qx=static_cast<double>(x)/static_cast<double>(n);
                const double qy=static_cast<double>(y)/static_cast<double>(n);
                const double qz=static_cast<double>(z)/static_cast<double>(n);
                m.points().set(vid(x,y,z),
                    qx + shear_xy*qy + shear_xz*qz, qy, qz);
            }

    std::unordered_map<std::string,std::size_t> face_ids;
    std::vector<std::vector<std::size_t>> face_cells;
    std::vector<std::vector<std::size_t>> cell_faces;

    auto get_face = [&](std::size_t a, std::size_t b, std::size_t c) {
        const std::string key=face_key(a,b,c);
        const auto it=face_ids.find(key);
        if (it != face_ids.end()) return it->second;
        const std::size_t id=face_cells.size();
        m.faces().push_face({a,b,c});
        face_ids.emplace(key,id);
        face_cells.emplace_back();
        return id;
    };

    for (std::size_t k=0;k<n;++k)
        for (std::size_t j=0;j<n;++j)
            for (std::size_t i=0;i<n;++i) {
                std::array<std::size_t,8> c{};
                for (unsigned q=0;q<8;++q)
                    c[q]=vid(i+kCorner[q][0],j+kCorner[q][1],k+kCorner[q][2]);

                for (const auto& t : kTets) {
                    const std::array<std::size_t,4> v{
                        c[t[0]],c[t[1]],c[t[2]],c[t[3]]};
                    const std::array<std::array<std::size_t,3>,4> fverts{{
                        {{v[0],v[1],v[2]}},
                        {{v[0],v[1],v[3]}},
                        {{v[0],v[2],v[3]}},
                        {{v[1],v[2],v[3]}}
                    }};
                    std::vector<std::size_t> faces;
                    faces.reserve(4);
                    for (const auto& fv : fverts) {
                        const std::size_t f=get_face(fv[0],fv[1],fv[2]);
                        faces.push_back(f);
                        face_cells[f].push_back(cell_faces.size());
                    }
                    cell_faces.push_back(std::move(faces));
                }
            }

    m.ownership().resize(m.n_faces());
    for (std::size_t f=0;f<m.n_faces();++f) {
        const auto& cells=face_cells[f];
        if (cells.size()==1) {
            m.ownership().set_owner(f,cells[0]);
            m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
        } else if (cells.size()==2) {
            m.ownership().set_owner(f,cells[0]);
            m.ownership().set_neighbour(f,static_cast<std::int64_t>(cells[1]));
        } else {
            throw std::runtime_error("tet mesh face has invalid adjacency");
        }
    }
    for (const auto& faces : cell_faces) m.cells().push_cell(faces);

    Grid g;
    g.mesh=std::move(m);
    g.geometry=make_geometry_cache(g.mesh);

    // Deep interior cells have no boundary faces, so the manufactured
    // divergence is not contaminated by a boundary closure policy.
    for (std::size_t k=1;k+1<n;++k)
        for (std::size_t j=1;j+1<n;++j)
            for (std::size_t i=1;i+1<n;++i) {
                const std::size_t base=6*(i+n*(j+n*k));
                for (std::size_t t=0;t<6;++t) g.interior.push_back(base+t);
            }
    return g;
}

struct AnalyticField {
    static double value(const Vec3& p) {
        const double pi=std::acos(-1.0);
        return std::sin(pi*p.x)*std::cos(pi*p.y)*std::sin(pi*p.z);
    }
    static Vec3 gradient(const Vec3& p) {
        const double pi=std::acos(-1.0);
        return Vec3{
            pi*std::cos(pi*p.x)*std::cos(pi*p.y)*std::sin(pi*p.z),
            -pi*std::sin(pi*p.x)*std::sin(pi*p.y)*std::sin(pi*p.z),
            pi*std::sin(pi*p.x)*std::cos(pi*p.y)*std::cos(pi*p.z)};
    }
};

Field<double,Location::CELL> sample_smooth(const Grid& g)
{
    Field<double,Location::CELL> f(g.mesh.n_cells(),"phi","1",1);
    for (std::size_t c=0;c<g.mesh.n_cells();++c)
        f(c)=AnalyticField::value(g.geometry.cell_centres[c]);
    return f;
}

Field<double,Location::CELL> sample_linear(const Grid& g)
{
    Field<double,Location::CELL> f(g.mesh.n_cells(),"phi","1",1);
    for (std::size_t c=0;c<g.mesh.n_cells();++c) {
        const Vec3 p=g.geometry.cell_centres[c];
        f(c)=2.0*p.x-1.25*p.y+0.75*p.z+0.2;
    }
    return f;
}

Field<double,Location::CELL> sample_step(const Grid& g)
{
    Field<double,Location::CELL> f(g.mesh.n_cells(),"phi","1",1);
    for (std::size_t c=0;c<g.mesh.n_cells();++c)
        f(c)=g.geometry.cell_centres[c].x < 0.5 ? 0.0 : 1.0;
    return f;
}

Field<double,Location::FACE> make_flux(const Grid& g, const Vec3& U)
{
    Field<double,Location::FACE> flux(g.mesh.n_faces(),"flux","m3/s",1);
    for (std::size_t f=0;f<g.mesh.n_faces();++f)
        flux(f)=U.dot(g.geometry.face_Sf[f]);
    return flux;
}

Field<double,Location::CELL> divergence_from_face(
    const Grid& g, const Field<double,Location::FACE>& flux,
    const Field<double,Location::FACE>& face_value)
{
    Field<double,Location::CELL> div(g.mesh.n_cells(),"div","1/s",1);
    const auto& own=g.mesh.ownership();
    const auto* faces=g.mesh.cells().faces_data();
    const auto* offsets=g.mesh.cells().offsets_data();
    for (std::size_t c=0;c<g.mesh.n_cells();++c) {
        double sum=0.0;
        for (Offset k=offsets[c]; k<offsets[c+1]; ++k) {
            const std::size_t f=faces[k];
            const double sign=own.owner(f)==c ? 1.0 : -1.0;
            sum += sign*flux(f)*face_value(f);
        }
        div(c)=sum/g.geometry.cell_volumes[c];
    }
    return div;
}

double interior_l2(const Grid& g, const Field<double,Location::CELL>& got,
                   const Field<double,Location::CELL>& exact)
{
    double e2=0.0, w=0.0;
    for (const std::size_t c : g.interior) {
        const double e=got(c)-exact(c);
        e2 += g.geometry.cell_volumes[c]*e*e;
        w += g.geometry.cell_volumes[c];
    }
    return std::sqrt(e2/w);
}

double interior_linf(const Grid& g, const Field<double,Location::CELL>& got,
                     const Field<double,Location::CELL>& exact)
{
    double e=0.0;
    for (const std::size_t c : g.interior)
        e=std::max(e,std::abs(got(c)-exact(c)));
    return e;
}

Field<double,Location::CELL> exact_advection(const Grid& g, const Vec3& U)
{
    Field<double,Location::CELL> e(g.mesh.n_cells(),"exact","1/s",1);
    for (std::size_t c=0;c<g.mesh.n_cells();++c)
        e(c)=U.dot(AnalyticField::gradient(g.geometry.cell_centres[c]));
    return e;
}

double observed_order(double e_coarse,double e_fine)
{
    return std::log(e_coarse/e_fine)/std::log(2.0);
}

void check_conservation_constant()
{
    const Grid g=make_tet_grid(8,0.35,0.20);
    const Vec3 U{1.0,0.35,-0.20};
    const auto flux=make_flux(g,U);
    Field<double,Location::CELL> constant(g.mesh.n_cells(),"phi","1",1);
    for (std::size_t c=0;c<g.mesh.n_cells();++c) constant(c)=3.0;

    for (const InterpScheme scheme : {
            InterpScheme::UPWIND, InterpScheme::LINEAR,
            InterpScheme::BLENDED, InterpScheme::QUICK,
            InterpScheme::QUICK_BOUNDED}) {
        const auto fv=interpolate_cell_to_face(
            constant,g.mesh,g.geometry,scheme,&flux,LimiterType::NONE,nullptr,1.0);
        const auto div=divergence_from_face(g,flux,fv);
        double total=0.0;
        double max_abs=0.0;
        for (std::size_t c=0;c<g.mesh.n_cells();++c) {
            total += g.geometry.cell_volumes[c]*div(c);
            max_abs=std::max(max_abs,std::abs(div(c)));
        }
        std::cout<<"N4_CONSERVATION scheme="<<to_string(scheme)
                 <<" total="<<total<<" max_abs="<<max_abs<<"\n";
        EXPECT_NEAR(total,0.0,1e-12);
        EXPECT_TRUE(max_abs<1e-11);
    }

    const auto grad=compute_gradient_least_squares_quadratic(constant,g.mesh);
    for (const LimiterType limiter : {
            LimiterType::MINMOD, LimiterType::VANLEER,
            LimiterType::SUPERBEE, LimiterType::VAN_ALBADA,
            LimiterType::MC}) {
        const auto fv=interpolate_cell_to_face(
            constant,g.mesh,g.geometry,InterpScheme::LIMITED,&flux,limiter,&grad);
        const auto div=divergence_from_face(g,flux,fv);
        double total=0.0;
        double max_abs=0.0;
        for (std::size_t c=0;c<g.mesh.n_cells();++c) {
            total += g.geometry.cell_volumes[c]*div(c);
            max_abs=std::max(max_abs,std::abs(div(c)));
        }
        std::cout<<"N4_CONSERVATION limiter="<<to_string(limiter)
                 <<" total="<<total<<" max_abs="<<max_abs<<"\n";
        EXPECT_NEAR(total,0.0,1e-12);
        EXPECT_TRUE(max_abs<1e-11);
    }
}

void check_quick_linear_exactness()
{
    const Grid g=make_tet_grid(12,0.35,0.20);
    const Vec3 U{0.9,-0.3,0.7};
    const auto flux=make_flux(g,U);
    const auto phi=sample_linear(g);
    const auto exact=Field<double,Location::CELL>(g.mesh.n_cells(),"exact","1/s",1);
    auto expected=exact;
    for (std::size_t c=0;c<g.mesh.n_cells();++c) expected(c)=U.dot(Vec3{2.0,-1.25,0.75});

    const auto quick=interpolate_cell_to_face(
        phi,g.mesh,g.geometry,InterpScheme::QUICK,&flux);
    const auto div=divergence_from_face(g,flux,quick);
    const double linf=interior_linf(g,div,expected);
    const double l2=interior_l2(g,div,expected);
    std::cout<<"N4_QUICK_LINEAR L2="<<l2<<" Linf="<<linf<<"\n";
    EXPECT_TRUE(linf<1e-10);
}

void check_smooth_order()
{
    const Vec3 U{0.8,-0.35,0.6};
    std::vector<double> upwind_error;
    std::vector<double> quick_error;

    for (const std::size_t n : {4u,8u,16u}) {
        const Grid g=make_tet_grid(n,0.35,0.20);
        const auto flux=make_flux(g,U);
        const auto phi=sample_smooth(g);
        const auto exact=exact_advection(g,U);

        const auto up_face=interpolate_cell_to_face(
            phi,g.mesh,g.geometry,InterpScheme::UPWIND,&flux);
        const auto q_face=interpolate_cell_to_face(
            phi,g.mesh,g.geometry,InterpScheme::QUICK,&flux);
        const auto up_div=divergence_from_face(g,flux,up_face);
        const auto q_div=divergence_from_face(g,flux,q_face);

        const double eu=interior_l2(g,up_div,exact);
        const double eq=interior_l2(g,q_div,exact);
        upwind_error.push_back(eu);
        quick_error.push_back(eq);
        std::cout<<std::setprecision(12)
                 <<"N4_SMOOTH n="<<n<<" upwind_L2="<<eu
                 <<" quick_L2="<<eq<<"\n";
    }

    const double p_up=observed_order(upwind_error[1],upwind_error[2]);
    const double p_q=observed_order(quick_error[1],quick_error[2]);
    std::cout<<"N4_ORDER upwind="<<p_up<<" quick="<<p_q<<"\n";

    EXPECT_TRUE(std::isfinite(p_up));
    EXPECT_TRUE(std::isfinite(p_q));
    EXPECT_TRUE(p_up>0.70);
    EXPECT_TRUE(p_q>1.50);
    EXPECT_TRUE(quick_error[2] < quick_error[1]);
}

void check_tvd_smooth_order()
{
    const Vec3 U{0.8,-0.35,0.6};
    const std::vector<LimiterType> limiters = {
        LimiterType::MINMOD, LimiterType::VANLEER,
        LimiterType::SUPERBEE, LimiterType::VAN_ALBADA,
        LimiterType::MC};

    for (const LimiterType limiter : limiters) {
        std::vector<double> errors;
        for (const std::size_t n : {4u,8u,16u}) {
            const Grid g=make_tet_grid(n,0.35,0.20);
            const auto flux=make_flux(g,U);
            const auto phi=sample_smooth(g);
            const auto exact=exact_advection(g,U);
            const auto grad=compute_gradient_least_squares_quadratic(phi,g.mesh);
            const auto face=interpolate_cell_to_face(
                phi,g.mesh,g.geometry,InterpScheme::LIMITED,&flux,limiter,&grad);
            const auto div=divergence_from_face(g,flux,face);
            const double e=interior_l2(g,div,exact);
            errors.push_back(e);
            std::cout<<"N4_TVD_ORDER limiter="<<to_string(limiter)
                     <<" n="<<n<<" L2="<<e<<"\n";
        }
        const double p=observed_order(errors[1],errors[2]);
        std::cout<<"N4_TVD_ORDER limiter="<<to_string(limiter)
                 <<" order="<<p<<"\n";
        EXPECT_TRUE(std::isfinite(p));
        EXPECT_TRUE(p>1.20);
    }
}

void check_bounded_reconstruction()
{
    const Grid g=make_tet_grid(8,0.35,0.20);
    const Vec3 U{1.0,0.2,0.4};
    const auto flux=make_flux(g,U);
    const auto step=sample_step(g);
    const auto grad=compute_gradient_least_squares_quadratic(step,g.mesh);

    for (const LimiterType limiter : {
            LimiterType::MINMOD, LimiterType::VANLEER,
            LimiterType::SUPERBEE, LimiterType::VAN_ALBADA,
            LimiterType::MC}) {
        const auto fv=interpolate_cell_to_face(
            step,g.mesh,g.geometry,InterpScheme::LIMITED,&flux,limiter,&grad);
        double worst=0.0;
        const auto& own=g.mesh.ownership();
        for (std::size_t f=0;f<g.mesh.n_faces();++f) {
            const auto nb=own.neighbour(f);
            if (nb<0) continue;
            const std::size_t n=static_cast<std::size_t>(nb);
            const double lo=std::min(step(own.owner(f)),step(n));
            const double hi=std::max(step(own.owner(f)),step(n));
            worst=std::max(worst,std::max(lo-fv(f),fv(f)-hi));
        }
        std::cout<<"N4_BOUNDED limiter="<<to_string(limiter)<<" worst="<<worst<<"\n";
        EXPECT_TRUE(worst<=1e-12);
    }

    const auto qbv=interpolate_cell_to_face(
        step,g.mesh,g.geometry,InterpScheme::QUICK_BOUNDED,&flux);
    double worst=0.0;
    for (std::size_t f=0;f<g.mesh.n_faces();++f) {
        const auto nb=g.mesh.ownership().neighbour(f);
        if (nb<0) continue;
        const std::size_t n=static_cast<std::size_t>(nb);
        const double lo=std::min(step(g.mesh.ownership().owner(f)),step(n));
        const double hi=std::max(step(g.mesh.ownership().owner(f)),step(n));
        worst=std::max(worst,std::max(lo-qbv(f),qbv(f)-hi));
    }
    std::cout<<"N4_BOUNDED scheme=quick_bounded worst="<<worst<<"\n";
    EXPECT_TRUE(worst<=1e-12);
}

void check_production_diagnostics()
{
    const Grid g=make_tet_grid(8,0.35,0.20);
    const Vec3 U{1.0,0.2,0.4};
    const auto flux=make_flux(g,U);
    const auto phi=sample_step(g);

    ConvectionDiagnostics d{};
    (void)compute_convection(phi,flux,g.mesh,InterpScheme::QUICK_BOUNDED,
                             LimiterType::NONE,&g.geometry,
                             GradientScheme::LEAST_SQUARES_QUADRATIC,1.0,&d);
    std::cout<<"N4_DIAGNOSTICS limited_faces="<<d.limited_faces
             <<" unbounded_faces="<<d.unbounded_faces
             <<" fallback_faces="<<d.fallback_faces
             <<" reconstruction_failures="<<d.reconstruction_failures<<"\n";
    EXPECT_TRUE(d.unbounded_faces==0);
    EXPECT_TRUE(d.limited_faces>0);
}

} // namespace

int main()
{
    run_case("n4_conservation_constant",check_conservation_constant);
    run_case("n4_quick_linear_exactness",check_quick_linear_exactness);
    run_case("n4_smooth_order",check_smooth_order);
    run_case("n4_tvd_smooth_order",check_tvd_smooth_order);
    run_case("n4_bounded_reconstruction",check_bounded_reconstruction);
    run_case("n4_production_diagnostics",check_production_diagnostics);
    return run_all();
}
