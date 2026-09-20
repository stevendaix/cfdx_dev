// M0.15-T04 — Visualization Module Tests

#include "cfdx/core/visualization/mesh_visualizer.h"
#include "common/test_harness.h"

using namespace cfdx::core::visualization;
using namespace cfdx::testing;

int main() {
    run_case("visualization_scalar_image", []() {
        MeshVisualizer viz;
        bool result = viz.visualize_scalar("mesh.msh", "field.h5", "poisson_result.png");
        EXPECT_TRUE(result);
    });

    run_case("visualization_mesh_only", []() {
        VizConfig cfg{"svg", "/tmp/cfdx_viz/", 256, true, false, true};
        MeshVisualizer viz(cfg);
        bool result = viz.visualize_mesh("cavity_with_patches.msh", "mesh_geometry.svg");
        EXPECT_TRUE(result);
    });

    run_case("visualization_poisson_comparison", []() {
        MeshVisualizer viz;
        bool result = viz.visualize_poisson_comparison("exact_16x16.dat", "computed_16x16.dat", "poisson_comparison.png");
        EXPECT_TRUE(result);
        std::printf("  PASS: Module visualisation Poisson opérationnel (M0.15-T04)\n");
    });

    return run_all();
}
