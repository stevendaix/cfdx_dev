#include "cfdx/physics/wall_quadratic_stencil.h"
#include "cfdx/physics/wall_quadratic_forces.h"
#include "common/test_harness.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <stdexcept>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;

static Mesh make_unit_cube()
{
    Mesh mesh;
    mesh.points().resize(8);
    mesh.points().set(0,0,0,0); mesh.points().set(1,1,0,0);
    mesh.points().set(2,1,1,0); mesh.points().set(3,0,1,0);
    mesh.points().set(4,0,0,1); mesh.points().set(5,1,0,1);
    mesh.points().set(6,1,1,1); mesh.points().set(7,0,1,1);
    mesh.faces().push_face({0,3,2,1});
    mesh.faces().push_face({4,5,6,7});
    mesh.faces().push_face({0,1,5,4});
    mesh.faces().push_face({3,7,6,2});
    mesh.faces().push_face({0,4,7,3});
    mesh.faces().push_face({1,2,6,5});
    mesh.ownership().resize(6);
    for (std::size_t face=0; face<6; ++face) {
        mesh.ownership().set_owner(face,0);
        mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
    }
    mesh.cells().push_cell({0,1,2,3,4,5});
    const char* names[]={"bottom","top","front","back","left","right"};
    const PatchType types[]={PatchType::WALL,PatchType::WALL,PatchType::INLET,
        PatchType::OUTLET,PatchType::SYMMETRY,PatchType::WALL};
    for (std::size_t face=0; face<6; ++face) {
        Patch patch;
        patch.name=names[face]; patch.type=types[face];
        patch.face_ids={static_cast<std::uint32_t>(face)};
        mesh.boundary().add_patch(patch);
    }
    return mesh;
}

static Mesh make_cartesian_3x3x3()
{
    constexpr std::size_t n=3, points_per_axis=n+1;
    Mesh mesh;
    mesh.points().resize(points_per_axis*points_per_axis*points_per_axis);
    const auto point_id=[](std::size_t i,std::size_t j,std::size_t k) {
        return static_cast<std::uint32_t>((i*points_per_axis+j)*points_per_axis+k);
    };
    for(std::size_t i=0;i<=n;++i)
        for(std::size_t j=0;j<=n;++j)
            for(std::size_t k=0;k<=n;++k)
                mesh.points().set(point_id(i,j,k),static_cast<double>(i),
                    static_cast<double>(j),static_cast<double>(k));

    std::map<std::array<std::uint32_t,4>,std::size_t> face_lookup;
    for(std::size_t i=0;i<n;++i)
        for(std::size_t j=0;j<n;++j)
            for(std::size_t k=0;k<n;++k) {
                const std::array<std::array<std::uint32_t,4>,6> faces{{
                    {{point_id(i,j,k),point_id(i,j+1,k),point_id(i,j+1,k+1),point_id(i,j,k+1)}},
                    {{point_id(i+1,j,k),point_id(i+1,j,k+1),point_id(i+1,j+1,k+1),point_id(i+1,j+1,k)}},
                    {{point_id(i,j,k),point_id(i,j,k+1),point_id(i+1,j,k+1),point_id(i+1,j,k)}},
                    {{point_id(i,j+1,k),point_id(i+1,j+1,k),point_id(i+1,j+1,k+1),point_id(i,j+1,k+1)}},
                    {{point_id(i,j,k),point_id(i+1,j,k),point_id(i+1,j+1,k),point_id(i,j+1,k)}},
                    {{point_id(i,j,k+1),point_id(i,j+1,k+1),point_id(i+1,j+1,k+1),point_id(i+1,j,k+1)}}
                }};
                std::vector<FaceIndex> cell_faces;
                for(auto vertices:faces) {
                    auto key=vertices;
                    std::sort(key.begin(),key.end());
                    const std::size_t cell=(i*n+j)*n+k;
                    const auto found=face_lookup.find(key);
                    if(found==face_lookup.end()) {
                        const std::size_t face=mesh.faces().n_faces();
                        mesh.faces().push_face({vertices[0],vertices[1],vertices[2],vertices[3]});
                        mesh.ownership().resize(face+1);
                        mesh.ownership().set_owner(face,static_cast<CellIndex>(cell));
                        mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
                        face_lookup.emplace(key,face);
                        cell_faces.push_back(static_cast<FaceIndex>(face));
                    } else {
                        const std::size_t face=found->second;
                        mesh.ownership().set_neighbour(face,static_cast<std::int64_t>(cell));
                        cell_faces.push_back(static_cast<FaceIndex>(face));
                    }
                }
                mesh.cells().push_cell(cell_faces);
            }

    Patch wall; wall.name="x_min"; wall.type=PatchType::WALL;
    Patch other; other.name="other_boundaries"; other.type=PatchType::SYMMETRY;
    for(std::size_t face=0;face<mesh.n_faces();++face) {
        if(mesh.ownership().neighbour(face)!=FaceOwnership::BOUNDARY) continue;
        const auto off=mesh.faces().offsets_data()[face];
        bool on_x_min=true;
        for(std::size_t v=off;v<mesh.faces().offsets_data()[face+1];++v) {
            const auto point=mesh.faces().vertices_data()[v];
            if(std::abs(mesh.points().x_data()[point])>1e-14) on_x_min=false;
        }
        (on_x_min ? wall.face_ids : other.face_ids).push_back(
            static_cast<FaceIndex>(face));
    }
    mesh.boundary().add_patch(wall);
    mesh.boundary().add_patch(other);
    return mesh;
}

static double quadratic_pressure(const Vec3& x)
{
    return 2.0+0.3*x.x-0.2*x.y+0.4*x.z+0.05*x.x*x.x+
        0.02*x.x*x.y-0.03*x.x*x.z+0.04*x.y*x.y+0.01*x.y*x.z+
        0.06*x.z*x.z;
}

static Vec3 quadratic_velocity(const Vec3& x)
{
    return {1.0+0.2*x.x+0.3*x.y-0.1*x.z+0.04*x.x*x.x+0.02*x.y*x.z,
        -0.5+0.1*x.x-0.2*x.y+0.4*x.z+0.03*x.y*x.y-0.01*x.x*x.z,
        0.7-0.3*x.x+0.2*x.y+0.1*x.z+0.02*x.z*x.z+0.015*x.x*x.y};
}

int main()
{
    run_case("mesh_adjacency_includes_exact_wall_sample", [] {
        const Mesh mesh=make_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        const std::vector<double> pressure{2.5};
        const std::vector<Vec3> velocity{{1.0,2.0,3.0}};
        const Vec3 wall_velocity{0.0,0.0,0.0};
        const auto samples=build_wall_quadratic_stencil(
            mesh,geometry,0,pressure,velocity,wall_velocity,2.5,2);
        EXPECT_TRUE(samples.size()==2);
        EXPECT_NEAR(samples.front().position.z,0.0,1e-14);
        EXPECT_NEAR(samples.front().velocity.x,0.0,1e-14);
        EXPECT_NEAR(samples.front().velocity.y,0.0,1e-14);
        EXPECT_NEAR(samples.front().velocity.z,0.0,1e-14);
        EXPECT_NEAR(samples.back().velocity.x,1.0,1e-14);
        bool rank_rejected=false;
        try {
            (void)reconstruct_wall_face_quadratic(
                samples,geometry.face_centres[0],geometry.face_area_vectors[0],1.0);
        } catch (const std::invalid_argument&) {
            rank_rejected=true;
        }
        EXPECT_TRUE(rank_rejected);
    });

    run_case("mesh_adjacency_reconstructs_quadratic_field_on_3d_mesh", [] {
        const Mesh mesh=make_cartesian_3x3x3();
        const auto geometry=build_fv_geometry(mesh);
        std::vector<double> pressure(mesh.n_cells());
        std::vector<Vec3> velocity(mesh.n_cells());
        for(std::size_t cell=0;cell<mesh.n_cells();++cell) {
            pressure[cell]=quadratic_pressure(geometry.cell_centres[cell]);
            velocity[cell]=quadratic_velocity(geometry.cell_centres[cell]);
        }
        const std::size_t wall_patch=mesh.boundary().find("x_min");
        EXPECT_TRUE(wall_patch<mesh.boundary().n_patches());
        std::size_t wall_face=mesh.n_faces();
        for(const auto face:mesh.boundary().patch(wall_patch).face_ids) {
            const auto& centre=geometry.face_centres[face];
            if(std::abs(centre.y-1.5)<1e-12 && std::abs(centre.z-1.5)<1e-12) {
                wall_face=face;
                break;
            }
        }
        EXPECT_TRUE(wall_face<mesh.n_faces());
        const Vec3 wall_velocity=quadratic_velocity(geometry.face_centres[wall_face]);
        const double wall_pressure=quadratic_pressure(geometry.face_centres[wall_face]);
        const auto samples=build_wall_quadratic_stencil(
            mesh,geometry,wall_face,pressure,velocity,wall_velocity,wall_pressure,6);
        EXPECT_TRUE(samples.size()==mesh.n_cells()+1);
        const auto reconstructed=reconstruct_wall_face_quadratic(
            samples,geometry.face_centres[wall_face],
            geometry.face_area_vectors[wall_face],0.7);
        const Vec3 x=geometry.face_centres[wall_face];
        EXPECT_NEAR(reconstructed.pressure,quadratic_pressure(x),2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[0].x,0.2+0.08*x.x,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[0].y,0.3+0.02*x.z,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[0].z,-0.1+0.02*x.y,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[1].x,0.1-0.01*x.z,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[1].y,-0.2+0.06*x.y,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[1].z,0.4-0.01*x.x,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[2].x,-0.3+0.015*x.y,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[2].y,0.2+0.015*x.x,2e-10);
        EXPECT_NEAR(reconstructed.velocity_gradient[2].z,0.1+0.04*x.z,2e-10);

        std::vector<double> wall_pressure_samples(mesh.n_faces(),0.0);
        std::vector<Vec3> wall_velocity_samples(mesh.n_faces());
        for(const auto face:mesh.boundary().patch(wall_patch).face_ids) {
            wall_pressure_samples[face]=quadratic_pressure(geometry.face_centres[face]);
            wall_velocity_samples[face]=quadratic_velocity(geometry.face_centres[face]);
        }
        const auto surface=reconstruct_quadratic_wall_surface(
            mesh,geometry,wall_patch,pressure,velocity,wall_pressure_samples,
            wall_velocity_samples,0.7,6);
        EXPECT_TRUE(surface.size()==mesh.boundary().patch(wall_patch).face_ids.size());
        for(const auto& sample:surface) {
            EXPECT_TRUE(sample.normal.x>0.999999);
            EXPECT_NEAR(sample.pressure,quadratic_pressure(
                Vec3{sample.position.x,sample.position.y,sample.position.z}),2e-10);
            EXPECT_TRUE(std::isfinite(sample.viscous_traction.x));
            EXPECT_TRUE(std::isfinite(sample.viscous_traction.y));
            EXPECT_TRUE(std::isfinite(sample.viscous_traction.z));
        }
        forces::Reference reference;
        reference.area=9.0;
        const auto integrated=forces::integrate(surface,reference);
        EXPECT_TRUE(std::isfinite(integrated.cd));
    });

    run_case("mesh_adjacency_rejects_non_wall_patch", [] {
        const Mesh mesh=make_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        bool rejected=false;
        try {
            (void)build_wall_quadratic_stencil(
                mesh,geometry,2,std::vector<double>{1.0},
                std::vector<Vec3>{{0,0,0}},Vec3{},1.0);
        } catch (const std::invalid_argument&) {
            rejected=true;
        }
        EXPECT_TRUE(rejected);
    });

    run_case("mesh_adjacency_rejects_invalid_ring_depth", [] {
        const Mesh mesh=make_unit_cube();
        const auto geometry=build_fv_geometry(mesh);
        bool rejected=false;
        try {
            (void)build_wall_quadratic_stencil(
                mesh,geometry,0,std::vector<double>{1.0},
                std::vector<Vec3>{{0,0,0}},Vec3{},1.0,0);
        } catch (const std::invalid_argument&) {
            rejected=true;
        }
        EXPECT_TRUE(rejected);
    });
    return run_all();
}
