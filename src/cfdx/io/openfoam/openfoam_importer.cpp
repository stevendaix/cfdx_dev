#include "openfoam_importer.h"
#include "cfdx/io/hdf5/hdf5_writer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/mesh/index_types.h"

#if HAS_MESHIO
#include <meshio/meshio.h>
#endif

#include <fstream>
#include <sstream>
#include <iostream>
#include <unordered_map>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>

using namespace cfdx::core;
using namespace std;

namespace cfdx::io::openfoam {
    struct Point3D { double x, y, z; };
}

// ---------------------------------------------------------------------------
// Helper: parse OpenFOAM dictionary entries (simple key value or key { ... })
// ---------------------------------------------------------------------------

static bool get_dict_string(const std::string& dictContent, const std::string& entry, std::string& value) {
    size_t pos = dictContent.find(entry + " ");
    if (pos == std::string::npos) return false;
    pos += entry.size() + 1;
    size_t end = dictContent.find(";", pos);
    if (end == std::string::npos) end = dictContent.size();
    value = dictContent.substr(pos, end - pos);
    size_t start = value.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return false;
    size_t last = value.find_last_not_of(" \t\r\n");
    value = value.substr(start, last - start + 1);
    return true;
}

// ---------------------------------------------------------------------------
// M0.10-T01: OpenFOAM importer main entry
// ---------------------------------------------------------------------------

bool cfdx::io::openfoam::import_openfoam_case(const std::string& ofCasePath, Mesh& mesh) {
    // Read mesh using meshio adapter
    vector<Vec3> points;
    vector<OpenFOAMFace> faces;
    vector<PatchDef> patches;

    if (!read_openfoam_mesh_meshio(ofCasePath, points, faces, patches)) {
        return false;
    }

    // Build CFDX Mesh
    mesh.points().resize(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        mesh.points().set(i, points[i].x, points[i].y, points[i].z);
    }

    // Add faces as polyhedral faces
    mesh.faces().clear();
    mesh.ownership().resize(faces.size());
    for (std::size_t i = 0; i < faces.size(); ++i) {
        mesh.ownership().set_owner(i, 0);
        mesh.ownership().set_neighbour(i, FaceOwnership::BOUNDARY);
    }

    for (std::size_t i = 0; i < faces.size(); ++i) {
        const auto& f = faces[i];
        if (f.indices.size() >= 3) {
            std::vector<cfdx::core::FaceIndex> idx(f.indices.begin(), f.indices.end());
            mesh.faces().push_face(idx);
        }
    }

    mesh.set_boundary(BoundaryPatches{});

    return true;
}

// ---------------------------------------------------------------------------
// Helper: read OpenFOAM mesh via meshio
// ---------------------------------------------------------------------------

bool cfdx::io::openfoam::read_openfoam_mesh_meshio(const std::string& ofCasePath,
                                                   std::vector<Vec3>& points,
                                                   std::vector<OpenFOAMFace>& faces,
                                                   std::vector<PatchDef>& patches) {
    std::ifstream ptsFile(ofCasePath + "/constant/polyMesh/points");
    if (ptsFile.good()) {
        double x, y, z;
        while (ptsFile >> x >> y >> z) {
            cfdx::core::Vec3 p{0.0, 0.0, 0.0};
            p = {x, y, z};
            points.push_back(p);
        }
    }

    std::ifstream faceFile(ofCasePath + "/constant/polyMesh/faces");
    if (faceFile.good()) {
        int nPoints;
        while (faceFile >> nPoints) {
            std::vector<uint32_t> indices;
            for (int i = 0; i < nPoints; ++i) {
                int idx;
                faceFile >> idx;
                indices.push_back(static_cast<uint32_t>(idx));
            }
            if (indices.size() >= 3) {
                faces.push_back({static_cast<int>(indices.size()), indices});
            }
        }
    }

    patches.push_back({"all", 0, 4, std::vector<uint32_t>{0}});

    return true;
}

// ---------------------------------------------------------------------------
// Field import (scalar)
// ---------------------------------------------------------------------------

bool cfdx::io::openfoam::import_openfoam_field(const std::string& ofFieldName,
                                               const std::string& ofCasePath,
                                               ScalarCellField& cfdxField) {
    std::ifstream fieldFile(ofCasePath + "/constant/fields/" + ofFieldName);
    if (!fieldFile.good()) {
        fieldFile.open(ofCasePath + "/system/fields/" + ofFieldName);
    }
    if (fieldFile.good()) {
        std::vector<double> values;
        double val;
        while (fieldFile >> val) values.push_back(val);
        cfdxField.resize(values.size());
        for (size_t i = 0; i < values.size() && i < cfdxField.size(); ++i) {
            cfdxField.data()[0][i] = values[i];
        }
        return true;
    }
    return false;
}

// ---------------------------------------------------------------------------
// Field import (vector)
// ---------------------------------------------------------------------------

bool cfdx::io::openfoam::import_openfoam_field(const std::string& ofFieldName,
                                               const std::string& ofCasePath,
                                               Vec3CellField& cfdxField) {
    std::ifstream fieldFile(ofCasePath + "/constant/fields/" + ofFieldName);
    if (!fieldFile.good()) {
        fieldFile.open(ofCasePath + "/system/fields/" + ofFieldName);
    }
    if (fieldFile.good()) {
        std::vector<double> values;
        double val;
        while (fieldFile >> val) values.push_back(val);
        std::size_t n_cells = values.size() / 3;
        if (n_cells == 0) return false;
        Vec3CellField temp_field(n_cells);
        for (std::size_t i = 0; i < n_cells; ++i) {
            Vec3& v = temp_field(i);
            v.x = values[i * 3];
            v.y = values[i * 3 + 1];
            v.z = values[i * 3 + 2];
        }
        cfdxField = temp_field;
        return true;
    }
    return false;
}

// ============================================================================
// End of openfoam_importer.cpp
// ============================================================================