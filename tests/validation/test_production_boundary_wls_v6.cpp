#include "cfdx/physics/steady_incompressible_solver.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "common/test_harness.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

namespace {
Mesh grid(std::size_t nx, std::size_t ny, std::size_t nz, double shear = 0.0,
          bool two_dimensional = false, bool normalized = false) {
    Mesh mesh;
    const auto v = [=](std::size_t i, std::size_t j, std::size_t k) {
        return (k * (ny + 1) + j) * (nx + 1) + i;
    };
    mesh.points().resize((nx+1)*(ny+1)*(nz+1));
    for (std::size_t k=0;k<=nz;++k)
        for (std::size_t j=0;j<=ny;++j)
            for (std::size_t i=0;i<=nx;++i)
                mesh.points().set(v(i,j,k),
                    static_cast<double>(i) / (normalized ? static_cast<double>(nx) : 1.0) +
                        shear * static_cast<double>(j) / (normalized ? static_cast<double>(ny) : 1.0),
                    static_cast<double>(j) / (normalized ? static_cast<double>(ny) : 1.0),
                    static_cast<double>(k) / (normalized ? static_cast<double>(nz) : 1.0));
    std::map<std::array<std::size_t,4>,std::size_t> seen;
    std::vector<std::vector<std::size_t>> cf(nx*ny*nz);
    const auto cindex = [=](std::size_t i,std::size_t j,std::size_t k) {
        return (k*ny+j)*nx+i;
    };
    for (std::size_t k=0;k<nz;++k)
        for (std::size_t j=0;j<ny;++j)
            for (std::size_t i=0;i<nx;++i) {
                const std::size_t cell=cindex(i,j,k);
                const auto add=[&](std::array<std::size_t,4> vertices) {
                    auto key=vertices;
                    std::sort(key.begin(),key.end());
                    auto it=seen.find(key);
                    if (it==seen.end()) {
                        const std::size_t id=mesh.faces().n_faces();
                        mesh.faces().push_face({vertices.begin(),vertices.end()});
                        mesh.ownership().resize(mesh.n_faces());
                        mesh.ownership().set_owner(id,cell);
                        mesh.ownership().set_neighbour(id,FaceOwnership::BOUNDARY);
                        seen.emplace(key,id);
                        cf[cell].push_back(id);
                    } else {
                        mesh.ownership().set_neighbour(it->second,
                            static_cast<std::int64_t>(cell));
                        cf[cell].push_back(it->second);
                    }
                };
                const auto a=v(i,j,k),b=v(i+1,j,k),c=v(i+1,j+1,k),d=v(i,j+1,k);
                const auto A=v(i,j,k+1),B=v(i+1,j,k+1),C=v(i+1,j+1,k+1),D=v(i,j+1,k+1);
                add({a,d,c,b}); add({A,B,C,D});
                add({a,b,B,A}); add({d,D,C,c});
                add({a,A,D,d}); add({b,c,C,B});
            }
    for (const auto& faces:cf) mesh.cells().push_cell(faces);
    mesh.boundary().add_patch("wall", PatchType::WALL);
    if (two_dimensional) mesh.boundary().add_patch("frontBack",PatchType::EMPTY);
    for (std::size_t f=0; f<mesh.n_faces(); ++f) {
        if (mesh.ownership().neighbour(f)>=0) continue;
        const auto first = mesh.faces().vertices_data()[mesh.faces().offsets_data()[f]];
        const auto second= mesh.faces().vertices_data()[mesh.faces().offsets_data()[f]+1];
        const std::size_t z_first= first/((nx+1)*(ny+1));
        const std::size_t z_second=second/((nx+1)*(ny+1));
        const bool on_empty = two_dimensional && (z_first == z_second);
        mesh.boundary().patch(on_empty?1:0).face_ids.push_back(f);
    }
    return mesh;
}

double analytical(const Vec3& c) {return 4.0 + 1.25*c.x - 0.8*c.y + 0.3*c.z;}
}

int main() {
    run_case("trapezoid_face_integral_centroid", [] {
        const std::vector<double> x={0.,3.,2.,0.};
        const std::vector<double> y={0.,0.,1.,1.};
        const std::vector<double> z={0.,0.,0.,0.};
        const VertexIndex face[]={0,1,2,3};
        const auto g=compute_face_geometry(x,y,z,face,0,4);
        EXPECT_NEAR(g.area,2.5,1.e-13);
        EXPECT_NEAR(g.centre.x,19.0/15.0,1.e-13);
        EXPECT_NEAR(g.centre.y,7.0/15.0,1.e-13);
        EXPECT_NEAR(g.warpage_ratio,0.0,1.e-13);
    });
    run_case("warped_quad_has_independent_area_warning_metric", [] {
        const std::vector<double> x={0.,1.,1.,0.};
        const std::vector<double> y={0.,0.,1.,1.};
        const std::vector<double> z={0.,0.,0.5,0.};
        const VertexIndex face[]={0,1,2,3};
        const auto g=compute_face_geometry(x,y,z,face,0,4);
        EXPECT_TRUE(g.triangulated_area > g.area);
        EXPECT_TRUE(g.warpage_ratio > 0.01);
    });
    run_case("pressure_neumann_is_active_without_skew_option", [] {
        auto mesh=grid(1,1,1);
        mesh.boundary().add_patch("east", PatchType::WALL);
        auto& wall=mesh.boundary().patch(0).face_ids;
        for (auto it=wall.begin();it!=wall.end();) {
            const auto f=*it;
            const auto face=compute_face_geometry(
                mesh.points().x_data(),mesh.points().y_data(),mesh.points().z_data(),
                mesh.faces().vertices_data(),mesh.faces().offsets_data()[f],
                mesh.faces().offsets_data()[f+1]-mesh.faces().offsets_data()[f]);
            if (face.centre.x > 0.99) {
                mesh.boundary().patch(1).face_ids.push_back(f);
                it=wall.erase(it);
            } else ++it;
        }
        const auto g=build_fv_geometry(mesh);
        Field<double,Location::CELL> field(1,"p","Pa"); field(0)=5.0;
        ScalarBoundaryConditions bcs;
        bcs["east"]={ScalarBoundaryType::FIXED_GRADIENT,0,2};
        auto default_grad=gauss_gradient_with_boundary(field,mesh,g,bcs,false);
        EXPECT_NEAR(default_grad(0,0),1.0,1.e-12);
        auto skew_grad=gauss_gradient_with_boundary(field,mesh,g,bcs,true);
        EXPECT_NEAR(default_grad(0,0),skew_grad(0,0),1.e-12);
    });
    run_case("skewed_hexa_affine_boundary_wls_exact", [] {
        auto mesh=grid(2,2,2,0.4);
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::CELL> phi(mesh.n_cells(),"phi","1");
        for (std::size_t c=0;c<mesh.n_cells();++c)
            phi(c)=analytical(geom.cell_centres[c]);
        ScalarBoundaryFaceConditions boundary(mesh.n_faces());
        for (std::size_t f=0;f<mesh.n_faces();++f)
            if (mesh.ownership().neighbour(f)<0) {
                boundary.specified[f]=1;
                boundary.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                                       analytical(geom.face_centres[f])};
            }
        BoundaryWlsDiagnostics metrics;
        auto gradients=compute_boundary_aware_wls_gradient(
            phi,mesh,geom,{},nullptr,&boundary,&metrics,1.e8,true);
        EXPECT_TRUE(metrics.robust_cells==mesh.n_cells());
        EXPECT_TRUE(metrics.fallback_cells==0);
        for (std::size_t c=0;c<mesh.n_cells();++c) {
            EXPECT_NEAR(gradients(c,0),1.25,2.e-12);
            EXPECT_NEAR(gradients(c,1),-0.8,2.e-12);
            EXPECT_NEAR(gradients(c,2),0.3,2.e-12);
        }
        auto deferred=compute_nonorthogonal_diffusion_face_fluxes(
            mesh,geom,phi,0.3,{},nullptr,nullptr,&boundary);
        double max_error=0;
        for (std::size_t f=0;f<mesh.n_faces();++f) {
            const auto o=mesh.ownership().owner(f);
            const auto nr=mesh.ownership().neighbour(f);
            const Vec3 d=nr>=0 ? geom.cell_centres[static_cast<std::size_t>(nr)]-geom.cell_centres[o]
                                 : geom.face_centres[f]-geom.cell_centres[o];
            const auto sf=geom.face_area_vectors[f];
            const double alpha=sf.mag2()/sf.dot(d);
            const Vec3 tangent=sf-d*alpha;
            const double expected=0.3*(1.25*tangent.x-0.8*tangent.y+0.3*tangent.z);
            max_error=std::max(max_error,std::abs(deferred[f]-expected));
        }
        EXPECT_TRUE(max_error < 1.e-10);
    });
    run_case("boundary_conditions_stay_active_with_three_internal_neighbours", [] {
        auto mesh=grid(2,2,2,0.35);
        const auto geom=build_fv_geometry(mesh);
        Field<double,Location::CELL> field(mesh.n_cells(),"T","K"); field.fill(0.0);
        ScalarBoundaryFaceConditions bcs(mesh.n_faces());
        for (std::size_t f=0;f<mesh.n_faces();++f)
            if (mesh.ownership().neighbour(f)<0 &&
                geom.face_centres[f].y > 1.99) {
                bcs.specified[f]=1;
                bcs.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE, 3.0};
            }
        BoundaryWlsDiagnostics diag;
        auto gradient=compute_boundary_aware_wls_gradient(
            field,mesh,geom,{},nullptr,&bcs,&diag);
        EXPECT_TRUE(diag.robust_cells==mesh.n_cells());
        EXPECT_TRUE(gradient(2,1) > 0.1);
        EXPECT_TRUE(gradient(3,1) > 0.1);
    });
    run_case("second_order_and_tvd_assembly_match_independent_face_flux", [] {
        auto mesh=grid(2,2,2,0.35);
        const auto geometry=build_fv_geometry(mesh);
        Field<double,Location::CELL> phi(mesh.n_cells(),"phi","1");
        Field<double,Location::CELL> source(mesh.n_cells(),"source","1");
        Field<double,Location::CELL> implicit(mesh.n_cells(),"implicit","1");
        Field<double,Location::FACE> mass(mesh.n_faces(),"mdot","kg/s");
        source.fill(0.0); implicit.fill(0.0); mass.fill(0.0);
        for (std::size_t c=0;c<mesh.n_cells();++c) {
            const auto& x=geometry.cell_centres[c];
            phi(c)=x.x*x.x+0.3*x.y-0.1*x.z;
        }
        for (std::size_t f=0;f<mesh.n_faces();++f)
            if (mesh.ownership().neighbour(f)>=0)
                mass(f)=0.1+0.02*static_cast<double>(f%7);
        ScalarBoundaryFaceConditions bc(mesh.n_faces());
        for (std::size_t f=0;f<mesh.n_faces();++f)
            if (mesh.ownership().neighbour(f)<0) {
                const auto& x=geometry.face_centres[f];
                bc.specified[f]=1;
                bc.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                    x.x*x.x+0.3*x.y-0.1*x.z};
            }
        for (auto scheme:{ConvectionScheme::SECOND_ORDER_UPWIND,ConvectionScheme::TVD}) {
            const auto eq=assemble_scalar_equation(
                mesh,geometry,mass,0.2,source,implicit,{},false,nullptr,
                nullptr,nullptr,nullptr,scheme,&phi,&bc);
            const auto flux=reconstruct_scalar_transport_flux(
                mesh,geometry,mass,phi,0.2,{},false,scheme,&bc);
            std::vector<double> face_div(mesh.n_cells(),0.0);
            for (std::size_t f=0;f<mesh.n_faces();++f) {
                const auto o=mesh.ownership().owner(f);
                face_div[o]+=flux(f);
                const auto n=mesh.ownership().neighbour(f);
                if (n>=0) face_div[static_cast<std::size_t>(n)]-=flux(f);
            }
            const auto* row=eq.matrix.row_offsets_data();
            const auto* col=eq.matrix.columns_data();
            const auto* a=eq.matrix.values_data();
            for (std::size_t c=0;c<mesh.n_cells();++c) {
                double residual=-eq.rhs(c);
                for (auto k=row[c];k<row[c+1];++k)
                    residual+=a[k]*phi(col[k]);
                EXPECT_NEAR(residual,face_div[c],1.e-11);
            }
        }
    });
    run_case("smooth_skewed_manufactured_diffusion_spatial_refinement", [] {
        constexpr double pi=3.14159265358979323846;
        std::vector<double> errors;
        std::vector<double> interior_errors;
        for (std::size_t n:{8u,16u,32u}) {
            auto mesh=grid(n,n,1,0.35,true,true);
            const auto geometry=build_fv_geometry(mesh);
            Field<double,Location::CELL> phi(mesh.n_cells(),"manufactured","1");
            Field<double,Location::FACE> mass(mesh.n_faces(),"mass","kg/s");
            mass.fill(0.0);
            const auto exact=[&](const Vec3& x) {
                return std::sin(pi*x.x)*std::sin(pi*x.y);
            };
            const auto gradient=[&](const Vec3& x) {
                return Vec3{pi*std::cos(pi*x.x)*std::sin(pi*x.y),
                            pi*std::sin(pi*x.x)*std::cos(pi*x.y),0.0};
            };
            for (std::size_t c=0;c<mesh.n_cells();++c)
                phi(c)=exact(geometry.cell_centres[c]);
            ScalarBoundaryFaceConditions bc(mesh.n_faces());
            for (std::size_t f=0;f<mesh.n_faces();++f)
                if (mesh.ownership().neighbour(f)<0 &&
                    mesh.boundary().patch(geometry.face_patch[f]).type!=PatchType::EMPTY) {
                    bc.specified[f]=1;
                    bc.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                                      exact(geometry.face_centres[f])};
                }
            const auto flux=reconstruct_scalar_transport_flux(
                mesh,geometry,mass,phi,1.0,{},false,ConvectionScheme::UPWIND,&bc);
            double squared=0.0, squared_core=0.0;
            std::size_t count=0, count_core=0;
            for (std::size_t f=0;f<mesh.n_faces();++f)
                if (mesh.ownership().neighbour(f)>=0) {
                    const Vec3 sf=geometry.face_area_vectors[f];
                    const double actual=-gradient(geometry.face_centres[f]).dot(sf);
                    const double error=(flux(f)-actual)/sf.mag();
                    squared+=error*error;
                    ++count;
                    const auto owner=mesh.ownership().owner(f);
                    const auto neighbour=static_cast<std::size_t>(mesh.ownership().neighbour(f));
                    const auto away_from_wall=[=](std::size_t c) {
                        const auto i=c%n, j=c/n;
                        return i>0 && i+1<n && j>0 && j+1<n;
                    };
                    if (away_from_wall(owner) && away_from_wall(neighbour)) {
                        squared_core+=error*error;
                        ++count_core;
                    }
                }
            EXPECT_TRUE(count>0 && count_core>0);
            errors.push_back(std::sqrt(squared/static_cast<double>(count)));
            interior_errors.push_back(std::sqrt(squared_core/static_cast<double>(count_core)));
            std::cout<<"MMS skewed diffusion n="<<n<<" flux_L2="<<errors.back()
                     <<" interior_L2="<<interior_errors.back()<<"\n";
        }
        const double order_1=std::log(errors[0]/errors[1])/std::log(2.0);
        const double order_2=std::log(errors[1]/errors[2])/std::log(2.0);
        const double core_order_1=std::log(interior_errors[0]/interior_errors[1])/std::log(2.0);
        const double core_order_2=std::log(interior_errors[1]/interior_errors[2])/std::log(2.0);
        std::cout<<"MMS full orders="<<order_1<<","<<order_2
                 <<" core orders="<<core_order_1<<","<<core_order_2<<"\n";
        // This is a flux verification (not SST flow validation). Keep the
        // measured order visible; only demand reliable improvement here.
        EXPECT_TRUE(errors[1]<errors[0]);
        EXPECT_TRUE(errors[2]<errors[1]);
    });
    run_case("two_dimensional_extrusion_rank_two_accepted_explicitly", [] {
        auto mesh=grid(2,2,1,0.25,true);
        auto geom=build_fv_geometry(mesh);
        Field<double,Location::CELL> phi(mesh.n_cells(),"T","K");
        for (std::size_t c=0;c<mesh.n_cells();++c)
            phi(c)=1.1*geom.cell_centres[c].x-0.75*geom.cell_centres[c].y;
        ScalarBoundaryFaceConditions bcs(mesh.n_faces());
        for (std::size_t f=0;f<mesh.n_faces();++f)
            if (mesh.ownership().neighbour(f)<0 &&
                mesh.boundary().patch(geom.face_patch[f]).type!=PatchType::EMPTY) {
                bcs.specified[f]=1;
                const auto cf=geom.face_centres[f];
                bcs.conditions[f]={ScalarBoundaryFaceCondition::Type::FIXED_VALUE,
                                    1.1*cf.x-0.75*cf.y};
            }
        BoundaryWlsDiagnostics d;
        const auto gr=compute_boundary_aware_wls_gradient(
            phi,mesh,geom,{},nullptr,&bcs,&d,1.e8,true);
        EXPECT_TRUE(d.robust_cells==mesh.n_cells());
        for (std::size_t c=0;c<mesh.n_cells();++c) {
            EXPECT_NEAR(gr(c,0),1.1,1.e-12);
            EXPECT_NEAR(gr(c,1),-0.75,1.e-12);
            EXPECT_NEAR(gr(c,2),0.0,1.e-12);
        }
    });
    return run_all();
}
