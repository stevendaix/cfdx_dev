#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "common/test_harness.h"
#include <cmath>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

static double volume_from_faces(
    const std::vector<Vec3>& points,
    const std::vector<std::vector<std::size_t>>& faces)
{
    std::vector<double> x(points.size()), y(points.size()), z(points.size());
    for (std::size_t i=0;i<points.size();++i) {
        x[i]=points[i].x; y[i]=points[i].y; z[i]=points[i].z;
    }
    std::vector<Vec3> centres(faces.size()), sf(faces.size());
    std::vector<FaceIndex> ids(faces.size());
    for (std::size_t f=0;f<faces.size();++f) {
        std::vector<FaceIndex> v;
        for (const auto i:faces[f]) v.push_back(i);
        const auto g=compute_face_geometry(
            x.data(),y.data(),z.data(),v.data(),0,v.size());
        centres[f]=g.centre; sf[f]=g.Sf; ids[f]=static_cast<FaceIndex>(f);
    }
    return compute_cell_geometry(
        centres.data(),sf.data(),ids.data(),ids.size()).volume;
}

int main()
{
    run_case("tetrahedron_volume", [] {
        const std::vector<Vec3> p={
            {0,0,0},{1,0,0},{0,1,0},{0,0,1}};
        const std::vector<std::vector<std::size_t>> f={
            {0,2,1},{0,1,3},{0,3,2},{1,2,3}};
        EXPECT_NEAR(volume_from_faces(p,f),1.0/6.0,1e-12);
    });

    run_case("pyramid_volume", [] {
        const std::vector<Vec3> p={
            {0,0,0},{1,0,0},{1,1,0},{0,1,0},{0.5,0.5,1}};
        const std::vector<std::vector<std::size_t>> f={
            {0,3,2,1},{0,1,4},{1,2,4},{2,3,4},{3,0,4}};
        EXPECT_NEAR(volume_from_faces(p,f),1.0/3.0,1e-12);
    });

    run_case("triangular_prism_volume", [] {
        const std::vector<Vec3> p={
            {0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}};
        const std::vector<std::vector<std::size_t>> f={
            {0,2,1},{3,4,5},{0,1,4,3},{1,2,5,4},{2,0,3,5}};
        EXPECT_NEAR(volume_from_faces(p,f),0.5,1e-12);
    });

    run_case("octahedron_arbitrary_polyhedron_volume", [] {
        const std::vector<Vec3> p={
            {1,0,0},{0,1,0},{-1,0,0},{0,-1,0},
            {0,0,1},{0,0,-1}};
        const std::vector<std::vector<std::size_t>> f={
            {0,1,4},{1,2,4},{2,3,4},{3,0,4},
            {1,0,5},{2,1,5},{3,2,5},{0,3,5}};
        EXPECT_NEAR(volume_from_faces(p,f),4.0/3.0,1e-12);
    });

    return run_all();
}
