#include "cfdx/physics/thermal_zones.h"
#include "common/test_harness.h"
#include <iostream>
#include <vector>
#include <cmath>
using namespace cfdx::core;
using namespace cfdx::physics;
using namespace cfdx::testing;
static Mesh chain_mesh(std::size_t n, double dx)
{
    Mesh m;
    m.points().resize(4 * (n + 1));
    for (std::size_t i = 0; i <= n; ++i) {
        const double x = static_cast<double>(i) * dx;
        m.points().set(4 * i + 0, x, 0.0, 0.0);
        m.points().set(4 * i + 1, x, 1.0, 0.0);
        m.points().set(4 * i + 2, x, 1.0, 1.0);
        m.points().set(4 * i + 3, x, 0.0, 1.0);
    }

    // x-normal faces: boundary i uses column i's four points.
    for (std::size_t i = 0; i <= n; ++i)
        m.faces().push_face({4 * i + 0, 4 * i + 1, 4 * i + 2, 4 * i + 3});
    // wall faces for cell i (columns i and i+1).
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t a = 4 * i;
        const std::size_t b = 4 * (i + 1);
        m.faces().push_face({a + 0, b + 0, b + 1, a + 1});  // y = 0
        m.faces().push_face({a + 3, a + 2, b + 2, b + 3});  // y = 1
        m.faces().push_face({a + 0, a + 3, b + 3, b + 0});  // z = 0
        m.faces().push_face({a + 1, b + 1, b + 2, a + 2});  // z = 1
    }

    const std::size_t nf = (n + 1) + 4 * n;
    m.ownership().resize(nf);
    // x-face 0: boundary of cell 0; x-face i (1..n-1): cells i-1 | i;
    // x-face n: boundary of cell n-1.
    m.ownership().set_owner(0, 0);
    m.ownership().set_neighbour(0, FaceOwnership::BOUNDARY);
    for (std::size_t i = 1; i < n; ++i) {
        m.ownership().set_owner(i, i - 1);
        m.ownership().set_neighbour(i, static_cast<int>(i));
    }
    m.ownership().set_owner(n, n - 1);
    m.ownership().set_neighbour(n, FaceOwnership::BOUNDARY);
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
            const std::size_t f = (n + 1) + 4 * i + j;
            m.ownership().set_owner(f, i);
            m.ownership().set_neighbour(f, FaceOwnership::BOUNDARY);
        }
    }

    for (std::size_t i = 0; i < n; ++i)
        m.cells().push_cell({i, i + 1,
                             (n + 1) + 4 * i + 0, (n + 1) + 4 * i + 1,
                             (n + 1) + 4 * i + 2, (n + 1) + 4 * i + 3});

    Patch left;
    left.name = "left";
    left.type = PatchType::WALL;
    left.face_ids = {0};
    Patch right;
    right.name = "right";
    right.type = PatchType::WALL;
    right.face_ids = {n};
    Patch walls;
    walls.name = "walls";
    walls.type = PatchType::WALL;
    walls.face_ids.reserve(4 * n);
    for (std::size_t f = n + 1; f < nf; ++f) walls.face_ids.push_back(f);
    m.boundary().add_patch(left);
    m.boundary().add_patch(right);
    m.boundary().add_patch(walls);
    return m;
}


// Two valid, geometrically closed, disconnected Cartesian components.
// A single cell (anchored) and two coupled cells (pure Neumann component).
static Mesh combine(const Mesh& a, const Mesh& b)
{
    Mesh m;
    const std::size_t np=a.n_points()+b.n_points();
    const std::size_t nf=a.n_faces()+b.n_faces();
    m.points().resize(np);
    for (std::size_t i=0;i<a.n_points();++i)
        m.points().set(i,a.points().x_data()[i],a.points().y_data()[i],a.points().z_data()[i]);
    for (std::size_t i=0;i<b.n_points();++i)
        m.points().set(a.n_points()+i,b.points().x_data()[i]+10.0,
                       b.points().y_data()[i],b.points().z_data()[i]);
    for (int part=0;part<2;++part) {
        const Mesh& x=part?b:a;
        const std::size_t po=part?a.n_points():0;
        for (std::size_t f=0;f<x.n_faces();++f) {
            std::vector<FaceIndex> vertices;
            for (auto k=x.faces().face_offset(f);k<x.faces().offsets_data()[f+1];++k)
                vertices.push_back(static_cast<FaceIndex>(x.faces().vertices_data()[k]+po));
            m.faces().push_face(vertices);
        }
    }
    m.ownership().resize(nf);
    for (int part=0;part<2;++part) {
        const Mesh& x=part?b:a;
        const std::size_t fo=part?a.n_faces():0;
        const std::size_t co=part?a.n_cells():0;
        for (std::size_t f=0;f<x.n_faces();++f) {
            m.ownership().set_owner(fo+f,static_cast<CellIndex>(x.ownership().owner(f)+co));
            const auto neighbour=x.ownership().neighbour(f);
            m.ownership().set_neighbour(fo+f,
                neighbour<0?FaceOwnership::BOUNDARY:neighbour+static_cast<std::int64_t>(co));
        }
        for (std::size_t c=0;c<x.n_cells();++c) {
            std::vector<FaceIndex> cf;
            for (auto k=x.cells().cell_offset(c);k<x.cells().offsets_data()[c+1];++k)
                cf.push_back(static_cast<FaceIndex>(x.cells().faces_data()[k]+fo));
            m.cells().push_cell(cf);
        }
    }
    Patch hot; hot.name="hot"; hot.type=PatchType::WALL;
    Patch insulated; insulated.name="insulated"; insulated.type=PatchType::WALL;
    for(std::size_t f=0;f<nf;++f) {
        if(m.ownership().neighbour(f)!=FaceOwnership::BOUNDARY) continue;
        if(f==0) hot.face_ids.push_back(static_cast<FaceIndex>(f));
        else insulated.face_ids.push_back(static_cast<FaceIndex>(f));
    }
    m.boundary().add_patch(hot);
    m.boundary().add_patch(insulated);
    return m;
}

int main()
{
    run_case("topologically_valid_disconnected_component_requires_own_dirichlet", [] {
        const Mesh mesh=combine(chain_mesh(1,1.0),chain_mesh(2,1.0));
        EXPECT_TRUE(mesh.topo_validate().ok);
        const auto geometry=build_fv_geometry(mesh);
        const auto zones=build_thermal_zones(mesh,std::vector<int>(3,0),
                                             {ThermalZoneSpec{"solid",0,1.0,0.0,0.0,true}});
        ScalarBoundaryConditions bc;
        bc["hot"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        EXPECT_THROW_WITH(assemble_thermal_conduction(mesh,geometry,zones,bc),
                          std::runtime_error, "unanchored connected component");
    });
    run_case("both_components_with_dirichlet_are_assembled", [] {
        const Mesh mesh=combine(chain_mesh(1,1.0),chain_mesh(2,1.0));
        const auto geometry=build_fv_geometry(mesh);
        const auto zones=build_thermal_zones(mesh,std::vector<int>(3,0),
                                             {ThermalZoneSpec{"solid",0,1.0,0.0,0.0,true}});
        ScalarBoundaryConditions bc;
        bc["hot"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        bc["insulated"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        const auto eq=assemble_thermal_conduction(mesh,geometry,zones,bc);
        EXPECT_TRUE(eq.matrix.n_rows()==3);
        for (std::size_t i=0;i<3;++i)
            EXPECT_TRUE(eq.matrix(i,i)>0.0);
    });
    run_case("single_connected_component_with_one_dirichlet_is_valid", [] {
        const Mesh mesh=chain_mesh(2,1.0);
        const auto geometry=build_fv_geometry(mesh);
        const auto zones=build_thermal_zones(mesh,std::vector<int>(2,0),
                                             {ThermalZoneSpec{"solid",0,1.0,0.0,0.0,true}});
        ScalarBoundaryConditions bc;
        bc["left"]={ScalarBoundaryType::FIXED_VALUE,300.0,0.0};
        const auto eq=assemble_thermal_conduction(mesh,geometry,zones,bc);
        EXPECT_TRUE(eq.matrix.n_rows()==2);
    });
    run_case("pure_neumann_rejection_unchanged", [] {
        const Mesh mesh=chain_mesh(2,1.0);
        const auto geometry=build_fv_geometry(mesh);
        const auto zones=build_thermal_zones(mesh,std::vector<int>(2,0),
                                             {ThermalZoneSpec{"solid",0,1.0,0.0,0.0,true}});
        EXPECT_THROW_WITH(assemble_thermal_conduction(mesh,geometry,zones,{}),
                          std::runtime_error, "no fixed-temperature boundary");
    });
    return run_all();
}
