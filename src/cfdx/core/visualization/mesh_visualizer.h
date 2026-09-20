// M0.15-T04 — Visualization Module for CFDX Results
// ===========================================================================
// Produces PNG/SVG visualizations from mesh fields and solver output.
// Based on Poisson MMS (Vertical Slice 1) pipeline results.

#pragma once

#include <string>
#include <vector>
#include <memory>

namespace cfdx {
namespace core {
namespace visualization {

// Visualization strategy
struct VizConfig {
    std::string output_format = "png";   // png | svg | vtu
    std::string output_path = "/tmp/";
    int resolution = 512;                // pixels for PNG
    bool show_mesh = true;
    bool show_field_contours = true;
    bool show_boundary_patches = true;
};

// Main visualization entry point
class MeshVisualizer {
public:
    MeshVisualizer(const VizConfig& cfg = VizConfig{});
    ~MeshVisualizer();

    // Generate image from a scalar cell field (e.g., Poisson pressure)
    bool visualize_scalar(const std::string& mesh_path,
                          const std::string& field_path,
                          const std::string& output_filename);

    // Generate mesh-only visualization (geometry, face normals, boundary patches)
    bool visualize_mesh(const std::string& mesh_path,
                        const std::string& output_filename);

    // Generate comparison image (Poisson MMS: exact vs computed)
    bool visualize_poisson_comparison(const std::string& exact_data,
                                       const std::string& computed_data,
                                       const std::string& output_filename);

private:
    VizConfig cfg_;
};

}  // namespace visualization
}  // namespace core
}  // namespace cfdx
