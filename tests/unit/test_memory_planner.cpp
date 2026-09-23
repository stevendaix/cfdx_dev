#include "cfdx/runtime/execution/memory_planner.h"
#include "common/test_harness.h"
#include <limits>

using namespace cfdx::testing;
using namespace cfdx::runtime;

int main() {
    run_case("planner_accounts_for_all_runtime_components", [] {
        MemoryPlanner planner({100, 50, 0.0});
        planner.add(MemoryKind::Mesh, 1000, "mesh");
        planner.add(MemoryKind::Fields, 2000, "fields");
        planner.add(MemoryKind::LinearAlgebra, 500, "matrix");
        const auto p = planner.plan(4000);
        EXPECT_TRUE(p.estimated_bytes == 3500);
        EXPECT_TRUE(p.runtime_reserve_bytes == 100);
        EXPECT_TRUE(p.safety_margin_bytes == 50);
        EXPECT_TRUE(p.peak_bytes == 3650);
        EXPECT_TRUE(p.fits());
    });

    run_case("planner_supports_fractional_safety_margin", [] {
        MemoryPlanner planner({0, 0, 0.10});
        planner.add(MemoryKind::Mesh, 1000);
        planner.add(MemoryKind::Fields, 1000);
        const auto p = planner.plan(2200);
        EXPECT_TRUE(p.estimated_bytes == 2000);
        EXPECT_TRUE(p.safety_margin_bytes == 200);
        EXPECT_TRUE(p.peak_bytes == 2200);
        EXPECT_TRUE(p.fits());
    });

    run_case("planner_reports_budget_exhaustion_without_backend_fallback", [] {
        MemoryPlanner planner({100, 100, 0.0});
        planner.add(MemoryKind::Mesh, 1000);
        const auto p = planner.plan(1100);
        EXPECT_TRUE(!p.fits());
        EXPECT_TRUE(p.peak_bytes == 1200);
    });

    run_case("planner_preserves_component_breakdown", [] {
        MemoryPlanner planner;
        planner.add(MemoryKind::Mesh, 1, "mesh");
        planner.add(MemoryKind::Geometry, 2, "geometry");
        planner.add(MemoryKind::Temporary, 3, "temporary");
        const auto p = planner.plan(10);
        EXPECT_TRUE(p.components.size() == 3);
        EXPECT_TRUE(p.components[1].kind == MemoryKind::Geometry);
        EXPECT_TRUE(p.components[2].bytes == 3);
    });

    run_case("planner_rejects_invalid_fraction", [] {
        EXPECT_THROW(MemoryPlanner({0, 0, 1.1}), std::invalid_argument);
    });

    run_case("planner_detects_component_overflow", [] {
        MemoryPlanner planner;
        planner.add(MemoryKind::Mesh, std::numeric_limits<std::size_t>::max());
        planner.add(MemoryKind::Fields, 1);
        const auto p = planner.plan(std::numeric_limits<std::size_t>::max());
        EXPECT_TRUE(p.overflow);
        EXPECT_TRUE(!p.fits());
    });

    run_case("planner_detects_peak_overflow", [] {
        MemoryPlanner planner({1, std::numeric_limits<std::size_t>::max(), 0.0});
        planner.add(MemoryKind::Mesh, 1);
        const auto p = planner.plan(std::numeric_limits<std::size_t>::max());
        EXPECT_TRUE(p.overflow);
        EXPECT_TRUE(!p.fits());
    });

    run_case("planner_converts_to_execution_workload", [] {
        MemoryPlanner planner({10, 20, 0.0});
        planner.add(MemoryKind::Mesh, 100);
        const auto p = planner.plan(1000);
        const auto w = make_runtime_workload(p);
        EXPECT_TRUE(w.estimated_bytes == 100);
        EXPECT_TRUE(w.runtime_reserve_bytes == 10);
        EXPECT_TRUE(w.safety_margin_bytes == 20);
    });

    return run_all();
}
