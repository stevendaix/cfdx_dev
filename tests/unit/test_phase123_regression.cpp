#include "cfdx/io/mesh/mesh_importer.h"
#include "cfdx/core/numerics/flux.h"
#include "cfdx/core/numerics/divergence.h"
#include "cfdx/core/mesh/mesh.h"
#include "common/test_harness.h"

using namespace cfdx::core;
using namespace cfdx::io::mesh;
using namespace cfdx::testing;

static Mesh cube() {
    Mesh m;
    m.points().resize(8);
    m.points().set(0,0,0,0); m.points().set(1,1,0,0);
    m.points().set(2,1,1,0); m.points().set(3,0,1,0);
    m.points().set(4,0,0,1); m.points().set(5,1,0,1);
    m.points().set(6,1,1,1); m.points().set(7,0,1,1);
    m.faces().push_face({0,3,2,1}); m.faces().push_face({4,5,6,7});
    m.faces().push_face({0,1,5,4}); m.faces().push_face({3,7,6,2});
    m.faces().push_face({0,4,7,3}); m.faces().push_face({1,2,6,5});
    m.ownership().resize(6);
    for(std::size_t f=0; f<6; ++f) {
        m.ownership().set_owner(f,0);
        m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);
    }
    m.cells().push_cell({0,1,2,3,4,5});
    return m;
}

int main() {
    run_case("universal_import_format_detection", [] {
        EXPECT_TRUE(detect_format("constant/polyMesh") == MeshFormat::MESHIO ||
                    detect_format("constant/polyMesh") == MeshFormat::UNKNOWN);
        EXPECT_TRUE(detect_format("mesh.msh") == MeshFormat::MESHIO);
        EXPECT_TRUE(detect_format("mesh.vtu") == MeshFormat::MESHIO);
    });

    run_case("phase3_closed_volume_flux_conservation", [] {
        Mesh m = cube();
        Field<double,Location::FACE> u(m.n_faces(),"U","m/s",3);
        u.fill(0.0);
        for(std::size_t f=0; f<m.n_faces(); ++f) u.set(f,1.0,0.0,0.0);
        auto phi=compute_flux(u,m);
        auto div=compute_divergence(phi,m);
        EXPECT_NEAR(phi(4)+phi(5),0.0,1e-12);
        EXPECT_NEAR(div(0),0.0,1e-12);
    });

    return run_all();
}
