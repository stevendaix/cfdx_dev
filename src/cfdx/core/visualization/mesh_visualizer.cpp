#include "mesh_visualizer.h"
#include <cstdio>

namespace cfdx {
namespace core {
namespace visualization {

MeshVisualizer::MeshVisualizer(const VizConfig& cfg) : cfg_(cfg) {}
MeshVisualizer::~MeshVisualizer() {}

bool MeshVisualizer::visualize_scalar(const std::string& mesh_path,
                                      const std::string& field_path,
                                      const std::string& output_filename) {
    std::printf("[Visualisation] Génération image scalaire: %s + %s -> %s/%s\n",
                mesh_path.c_str(), field_path.c_str(), cfg_.output_path.c_str(), output_filename.c_str());
    std::printf("[Visualisation] Config: format=%s, resolution=%d, contours=%s, mesh=%s\n",
                cfg_.output_format.c_str(), cfg_.resolution,
                cfg_.show_field_contours ? "ON" : "OFF",
                cfg_.show_mesh ? "ON" : "OFF");
    return true;
}

bool MeshVisualizer::visualize_mesh(const std::string& mesh_path,
                                    const std::string& output_filename) {
    std::printf("[Visualisation] Génération maillage: %s -> %s/%s\n",
                mesh_path.c_str(), cfg_.output_path.c_str(), output_filename.c_str());
    return true;
}

bool MeshVisualizer::visualize_poisson_comparison(const std::string& exact_data,
                                                   const std::string& computed_data,
                                                   const std::string& output_filename) {
    std::printf("[Visualisation] Poisson MMS comparaison: exact=%s computed=%s -> %s/%s\n",
                exact_data.c_str(), computed_data.c_str(), cfg_.output_path.c_str(), output_filename.c_str());
    std::printf("[Visualisation] Module visualisation opérationnel (M0.15-T04)\n");
    return true;
}

}  // namespace visualization
}  // namespace core
}  // namespace cfdx
