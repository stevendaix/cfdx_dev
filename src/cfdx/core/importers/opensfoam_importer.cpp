// OpenFOAM importer for CFDX geometries
// Converts OpenFOAM mesh files to CFDX internal representation
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include <iostream>
#include <fstream>
#include <stdexcept>

namespace cfdx {
namespace core {

struct OpensFoamImporter {
    void run(const std::string& input_file, std::string& output_dir) {
        std::cout << "Opening OpenFOAM importer for " << input_file << "..." << std::endl;
        
        (void)output_dir;
        throw std::runtime_error("OpenFOAM importer is not implemented; refusing to report a successful import");
    }
};

} // namespace core
} // namespace cfdx
