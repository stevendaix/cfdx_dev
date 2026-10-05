#include "cfdx/core/geometry/mesh_validator.h"
#include "cfdx/io/mesh/mesh_importer.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using cfdx::core::Mesh;
using cfdx::core::MeshQualityReport;

struct FixtureResult {
    std::string id;
    std::string source_path;
    bool imported = false;
    bool topology_valid = false;
    bool geometry_quality_valid = false;
    std::string error;
    cfdx::core::MeshStats stats{};
    MeshQualityReport quality{};
};

static std::string json_escape(const std::string& value) {
    std::string out;
    for (char c : value) {
        switch (c) {
        case '\\': out += "\\\\"; break;
        case '"': out += "\\\""; break;
        case '\n': out += "\\n"; break;
        case '\r': out += "\\r"; break;
        case '\t': out += "\\t"; break;
        default: out += c; break;
        }
    }
    return out;
}

static void write_report(const fs::path& path,
                         const std::vector<FixtureResult>& results) {
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot write report: " + path.string());

    bool all_pass = true;
    for (const auto& result : results)
        all_pass = all_pass && result.imported && result.topology_valid &&
                   result.geometry_quality_valid;

    out << "{\n";
    out << "  \"campaign\": \"N10 production fixture import/geometry qualification\",\n";
    out << "  \"status\": \"" << (all_pass ? "PASS" : "FAIL") << "\",\n";
    out << "  \"fixtures\": [\n";

    for (std::size_t i = 0; i < results.size(); ++i) {
        const auto& r = results[i];
        out << "    {\n";
        out << "      \"id\": \"" << json_escape(r.id) << "\",\n";
        out << "      \"source_path\": \"" << json_escape(r.source_path) << "\",\n";
        out << "      \"imported\": " << (r.imported ? "true" : "false") << ",\n";
        out << "      \"topology_valid\": " << (r.topology_valid ? "true" : "false") << ",\n";
        out << "      \"geometry_quality_valid\": "
            << (r.geometry_quality_valid ? "true" : "false") << ",\n";
        out << "      \"counts\": {\n";
        out << "        \"points\": " << r.stats.n_points << ",\n";
        out << "        \"faces\": " << r.stats.n_faces << ",\n";
        out << "        \"cells\": " << r.stats.n_cells << ",\n";
        out << "        \"boundary_faces\": " << r.stats.n_boundary_faces << ",\n";
        out << "        \"internal_faces\": " << r.stats.n_internal_faces << ",\n";
        out << "        \"patches\": " << r.stats.n_patches << "\n";
        out << "      },\n";
        out << "      \"quality\": {\n";
        out << "        \"max_skewness\": " << r.quality.max_skewness << ",\n";
        out << "        \"max_nonorthogonality_deg\": "
            << r.quality.max_non_orthogonality_deg << ",\n";
        out << "        \"max_aspect_ratio\": " << r.quality.max_aspect_ratio << ",\n";
        out << "        \"min_volume_ratio\": " << r.quality.min_volume_ratio << ",\n";
        out << "        \"min_cell_volume\": " << r.quality.min_cell_volume << ",\n";
        out << "        \"max_cell_volume\": " << r.quality.max_cell_volume << ",\n";
        out << "        \"surface_closure_error\": "
            << r.quality.surface_closure_error << "\n";
        out << "      },\n";
        out << "      \"errors\": [";
        for (std::size_t e = 0; e < r.quality.errors.size(); ++e) {
            if (e) out << ", ";
            out << "\"" << json_escape(r.quality.errors[e]) << "\"";
        }
        out << "]\n";
        out << "    }" << (i + 1 == results.size() ? "" : ",") << "\n";
    }

    out << "  ],\n";
    out << "  \"policy\": {\n";
    out << "    \"uses_production_importer\": true,\n";
    out << "    \"changes_numerical_tolerances\": false,\n";
    out << "    \"silent_fallbacks\": false,\n";
    out << "    \"solver_qualification\": \"deferred_to_next N10 increment\"\n";
    out << "  }\n";
    out << "}\n";
}

static FixtureResult qualify_fixture(const std::string& id, const fs::path& path) {
    FixtureResult result{id, path.string()};
    try {
        Mesh mesh;
        result.imported = cfdx::io::mesh::import_mesh(path.string(), mesh);
        if (!result.imported) {
            result.error = "production mesh importer rejected fixture";
            return result;
        }

        const auto topology = mesh.topo_validate();
        result.topology_valid = topology.ok;
        if (!result.topology_valid) {
            for (const auto& error : topology.errors) {
                if (!result.error.empty()) result.error += "; ";
                result.error += error;
            }
            return result;
        }

        result.stats = mesh.stats();
        result.quality = cfdx::core::validate_mesh(mesh);
        result.geometry_quality_valid = result.quality.ok;
        if (!result.geometry_quality_valid) {
            for (const auto& error : result.quality.errors) {
                if (!result.error.empty()) result.error += "; ";
                result.error += error;
            }
        }
    } catch (const std::exception& exc) {
        result.error = exc.what();
    }
    return result;
}

int main(int argc, char** argv) {
    if (argc != 3) {
        std::cerr << "usage: test_n10_external_fixture_import <fixture-root> <report>\n";
        return 2;
    }

    const fs::path fixture_root = fs::absolute(argv[1]);
    const fs::path report_path = fs::absolute(argv[2]);

    const std::vector<std::pair<std::string, fs::path>> file_fixtures = {
        {"meshio-su2-square", fixture_root / "meshio-su2-square"},
        {"meshio-gmsh-insulated-2-2", fixture_root / "meshio-gmsh-insulated-2-2"},
        {"meshio-vtk-unstructured", fixture_root / "meshio-vtk-unstructured"},
    };

    std::vector<FixtureResult> results;
    for (const auto& [id, path] : file_fixtures)
        results.push_back(qualify_fixture(id, path));

    // The OpenFOAM fixture is a pinned polyMesh subtree, while the production
    // importer expects a case root containing constant/polyMesh.
    const fs::path staging =
        fs::temp_directory_path() / "cfdx_n10_openfoam_fixture_case";
    fs::remove_all(staging);
    fs::create_directories(staging / "constant");
    fs::create_directory_symlink(
        fixture_root / "openfoam-airfoil2d",
        staging / "constant" / "polyMesh");
    results.push_back(qualify_fixture("openfoam-airfoil2d", staging));
    fs::remove_all(staging);

    fs::create_directories(report_path.parent_path());
    write_report(report_path, results);

    bool pass = true;
    for (const auto& result : results) {
        std::cout << "N10_FIXTURE id=" << result.id
                  << " imported=" << (result.imported ? 1 : 0)
                  << " topology_valid=" << (result.topology_valid ? 1 : 0)
                  << " geometry_quality_valid="
                  << (result.geometry_quality_valid ? 1 : 0)
                  << " cells=" << result.stats.n_cells
                  << " faces=" << result.stats.n_faces
                  << " max_skewness=" << result.quality.max_skewness
                  << " max_nonorth_deg=" << result.quality.max_non_orthogonality_deg
                  << " max_aspect=" << result.quality.max_aspect_ratio
                  << " min_volume=" << result.quality.min_cell_volume
                  << " max_volume=" << result.quality.max_cell_volume
                  << "\n";
        if (!result.imported || !result.topology_valid ||
            !result.geometry_quality_valid) {
            pass = false;
            std::cerr << "N10 fixture failed: " << result.id;
            if (!result.error.empty()) std::cerr << ": " << result.error;
            std::cerr << "\n";
        }
    }

    std::cout << "N10 fixture qualification report: " << report_path << "\n";
    return pass ? 0 : 1;
}
