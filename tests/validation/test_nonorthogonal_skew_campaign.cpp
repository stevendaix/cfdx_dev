// Phase 3.6 — quantitative non-orthogonal/skewed Laplacian campaign.

#include "cfdx/core/numerics/laplacian.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "common/test_harness.h"
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::testing;

namespace {
Mesh make_two_cell_unit_cubes()
{
    Mesh m;
    m.points().resize(12);
    const double p[12][3] = {{0,0,0},{1,0,0},{2,0,0},{0,1,0},{1,1,0},{2,1,0},
                             {0,0,1},{1,0,1},{2,0,1},{0,1,1},{1,1,1},{2,1,1}};
    for (std::size_t i=0;i<12;++i) m.points().set(i,p[i][0],p[i][1],p[i][2]);
    m.faces().push_face({0,6,9,3}); m.faces().push_face({0,1,7,6});
    m.faces().push_face({3,9,10,4}); m.faces().push_face({0,3,4,1});
    m.faces().push_face({6,7,10,9}); m.faces().push_face({7,10,4,1});
    m.faces().push_face({2,5,11,8}); m.faces().push_face({1,2,8,7});
    m.faces().push_face({4,10,11,5}); m.faces().push_face({1,4,5,2});
    m.faces().push_face({7,8,11,10});
    m.ownership().resize(11);
    for(std::size_t f=0;f<5;++f){m.ownership().set_owner(f,0);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.ownership().set_owner(5,0); m.ownership().set_neighbour(5,1);
    for(std::size_t f=6;f<11;++f){m.ownership().set_owner(f,1);m.ownership().set_neighbour(f,FaceOwnership::BOUNDARY);}
    m.cells().push_cell({0,1,2,3,4,5}); m.cells().push_cell({5,6,7,8,9,10});
    return m;
}
struct Metrics { double orth; double corrected; double conservation; };
Metrics evaluate(double skew)
{
    const Mesh mesh=make_two_cell_unit_cubes();
    ScalarCellField field(2,"phi","1",1); field(0)=0.5; field(1)=1.5;
    GeometryCache geometry=make_geometry_cache(mesh);
    geometry.face_Sf[5].y=skew;
    geometry.face_Sf[2].y=1.0+skew;
    const auto orth=compute_laplacian(field,mesh,geometry,LaplacianScheme::ORTHOGONAL);
    const auto corrected=compute_laplacian(field,mesh,geometry,LaplacianScheme::CORRECTED);
    return {orth(0),corrected(0),std::abs(geometry.cell_volumes[0]*corrected(0)+geometry.cell_volumes[1]*corrected(1))};
}
}
int main()
{
    const std::vector<double> skews={0.0,0.05,0.10,0.20,0.30,0.40};
    std::cout<<std::setprecision(17);
    double previous=-1.0;
    for(const double skew:skews){
        run_case("phase3_6_skew_"+std::to_string(skew),[skew,&previous](){
            const Metrics m=evaluate(skew); const double correction=std::abs(m.corrected-m.orth);
            EXPECT_TRUE(std::isfinite(m.orth)); EXPECT_TRUE(std::isfinite(m.corrected));
            EXPECT_TRUE(std::isfinite(m.conservation)); EXPECT_NEAR(m.conservation,0.0,1e-12);
            if(skew==0.0) EXPECT_NEAR(correction,0.0,1e-12);
            else { EXPECT_TRUE(correction>0.0); if(previous>=0.0) EXPECT_TRUE(correction>=previous); }
            previous=correction;
            std::cout<<"PHASE3_6 skew="<<skew<<" orth="<<m.orth<<" corrected="<<m.corrected
                     <<" correction="<<correction<<" conservation="<<m.conservation<<"\\n";
        });
    }
    return run_all();
}
