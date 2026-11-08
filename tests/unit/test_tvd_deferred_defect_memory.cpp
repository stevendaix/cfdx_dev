#include "cfdx/physics/tvd_deferred_defect_memory.h"
#include "cfdx/core/mesh/mesh.h"
// The regression assertions must remain active in Release (-DNDEBUG).
#define REQUIRE(expr) do { \
    if (!(expr)) throw std::runtime_error("REQUIRE failed: " #expr); \
} while (false)
#include <cmath>
#include <iostream>
#include <stdexcept>
using namespace cfdx::core;
using namespace cfdx::physics;
int main() {
    Mesh m;
    m.points().resize(4);
    m.faces().push_face({0,1,2,3});
    m.ownership().resize(1);
    m.ownership().set_owner(0,0);
    m.ownership().set_neighbour(0,1);
    m.cells().push_cell({0});
    m.cells().push_cell({0});
    ScalarEquation eq{SparseMatrix(2,2),Vector(2,0.0),{1.0,1.0},0.0};
    eq.rhs(0)=-3.0; eq.rhs(1)=3.0; // explicit conservative TVD defect
    const std::vector<double> current{3.0}, previous{1.0};
    auto result=apply_tvd_deferred_defect_memory(m,eq,current,&previous,0.5);
    REQUIRE(std::abs(eq.rhs(0)+2.0)<1e-14);
    REQUIRE(std::abs(eq.rhs(1)-2.0)<1e-14);
    REQUIRE(std::abs(eq.rhs(0)+eq.rhs(1))<1e-14);
    REQUIRE(std::abs(result.used[0]-2.0)<1e-14);
    REQUIRE(std::abs(result.current_raw[0]-3.0)<1e-14);
    REQUIRE(std::abs(result.max_cell_source_difference-1.0)<1e-14);
    // At a fixed point, old and current defect are identical, so the exact
    // full-order (unfiltered) equation is preserved to machine precision.
    ScalarEquation fixed{SparseMatrix(2,2),Vector(2,0.0),{1.,1.},0.0};
    fixed.rhs(0)=-3.; fixed.rhs(1)=3.;
    auto identical=apply_tvd_deferred_defect_memory(m,fixed,current,&current,0.99);
    REQUIRE(identical.max_face_difference==0.0);
    REQUIRE(fixed.rhs(0)==-3. && fixed.rhs(1)==3.);
    ScalarEquation original{SparseMatrix(2,2),Vector(2,0.0),{1.,1.},0.0};
    original.rhs(0)=-3.; original.rhs(1)=3.;
    auto off=apply_tvd_deferred_defect_memory(m,original,current,nullptr,0.5);
    REQUIRE(off.max_cell_source_difference==0 && original.rhs(0)==-3.);
    // Face-wise oscillation filtering must NOT modify a monotone sequence.
    ScalarEquation oscillating{SparseMatrix(2,2),Vector(2,0.0),{1.,1.},0.0};
    oscillating.rhs(0)=-3.; oscillating.rhs(1)=3.;
    const std::vector<double> two_ago{3.0}; // 3 -> 1 -> 3 : reversal
    auto smoothed = apply_tvd_deferred_defect_memory(
        m,oscillating,current,&previous,0.25,&two_ago,true);
    REQUIRE(smoothed.filtered_faces==1);
    REQUIRE(std::abs(oscillating.rhs(0)+2.5)<1e-14);
    REQUIRE(std::abs(oscillating.rhs(1)-2.5)<1e-14);
    ScalarEquation monotone{SparseMatrix(2,2),Vector(2,0.0),{1.,1.},0.0};
    monotone.rhs(0)=-3.; monotone.rhs(1)=3.;
    const std::vector<double> further_old{0.5}; // 0.5 -> 1 -> 3
    auto unchanged=apply_tvd_deferred_defect_memory(
        m,monotone,current,&previous,0.25,&further_old,true);
    REQUIRE(unchanged.filtered_faces==0);
    REQUIRE(monotone.rhs(0)==-3.0 && monotone.rhs(1)==3.0);
    bool threw=false;
    try { (void)apply_tvd_deferred_defect_memory(m,original,{3.},&previous,1.0); }
    catch (const std::invalid_argument&) { threw=true; }
    REQUIRE(threw);
    TvdDeferredDefectMemoryControls invalid_combo;
    invalid_combo.recursive=true; invalid_combo.oscillation_only=true;
    bool combo_rejected=false;
    try { validate_tvd_deferred_defect_memory(invalid_combo); }
    catch (const std::invalid_argument&) { combo_rejected=true; }
    REQUIRE(combo_rejected);
    std::cout<<"PASS conservative TVD memory, fixed-point invariance, fail closed\n";
}

