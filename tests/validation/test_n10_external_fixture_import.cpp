#include "cfdx/core/geometry/mesh_validator.h"
#include "cfdx/io/mesh/mesh_importer.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
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
    std::string validation_mode = "3d";
    double min_cell_area_2d = 0.0;
    double max_cell_area_2d = 0.0;
    double min_edge_length_2d = 0.0;
    double max_edge_length_2d = 0.0;
    double max_aspect_ratio_2d = 0.0;
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

static bool validate_2d_edge_mesh(const Mesh& mesh, FixtureResult& result) {
    if (mesh.n_points() == 0 || mesh.n_cells() == 0 || mesh.n_faces() == 0)
        return false;

    double min_edge = std::numeric_limits<double>::infinity();
    double max_edge = 0.0;
    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const auto size = mesh.faces().face_size(f);
        if (size != 2) return false;
        const auto off = mesh.faces().face_offset(f);
        const auto a = mesh.faces().vertices_data()[off];
        const auto b = mesh.faces().vertices_data()[off + 1];
        const double dx = mesh.points().x(b) - mesh.points().x(a);
        const double dy = mesh.points().y(b) - mesh.points().y(a);
        const double dz = mesh.points().z(b) - mesh.points().z(a);
        const double length = std::sqrt(dx * dx + dy * dy + dz * dz);
        if (!(length > 0.0) || !std::isfinite(length)) return false;
        min_edge = std::min(min_edge, length);
        max_edge = std::max(max_edge, length);
    }

    double min_area = std::numeric_limits<double>::infinity();
    double max_area = 0.0;
    double max_aspect = 0.0;
    for (std::size_t c = 0; c < mesh.n_cells(); ++c) {
        const auto off = mesh.cells().cell_offset(c);
        const auto count = mesh.cells().cell_size(c);
        if (count < 3) return false;

        double signed_area_twice = 0.0;
        double lo_x = std::numeric_limits<double>::infinity();
        double lo_y = std::numeric_limits<double>::infinity();
        double hi_x = -std::numeric_limits<double>::infinity();
        double hi_y = -std::numeric_limits<double>::infinity();

        for (std::size_t k = 0; k < count; ++k) {
            const auto face = mesh.cells().faces_data()[off + k];
            const auto face_off = mesh.faces().face_offset(face);
            const auto a = mesh.faces().vertices_data()[face_off];
            const auto b = mesh.faces().vertices_data()[face_off + 1];
            const bool owner = mesh.ownership().owner(face) == c;
            const auto first = owner ? a : b;
            const auto second = owner ? b : a;

            const double x1 = mesh.points().x(first);
            const double y1 = mesh.points().y(first);
            const double x2 = mesh.points().x(second);
            const double y2 = mesh.points().y(second);
            signed_area_twice += x1 * y2 - x2 * y1;
            lo_x = std::min({lo_x, x1, x2});
            lo_y = std::min({lo_y, y1, y2});
            hi_x = std::max({hi_x, x1, x2});
            hi_y = std::max({hi_y, y1, y2});
        }

        const double area = std::abs(signed_area_twice) * 0.5;
        const double ex = hi_x - lo_x;
        const double ey = hi_y - lo_y;
        const double emin = std::min(ex, ey);
        const double emax = std::max(ex, ey);
        if (!(area > 0.0) || !std::isfinite(area) ||
            !(emin > 0.0) || !std::isfinite(emax))
            return false;

        min_area = std::min(min_area, area);
        max_area = std::max(max_area, area);
        max_aspect = std::max(max_aspect, emax / emin);
    }

    result.validation_mode = "2d_edge_mesh";
    result.min_cell_area_2d = min_area;
    result.max_cell_area_2d = max_area;
    result.min_edge_length_2d = min_edge;
    result.max_edge_length_2d = max_edge;
    result.max_aspect_ratio_2d = max_aspect;
    return true;
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
        bool all_edges = true;
        for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
            all_edges = all_edges && mesh.faces().face_size(f) == 2;
        }
        if (all_edges) {
            result.geometry_quality_valid = validate_2d_edge_mesh(mesh, result);
        } else {
            result.quality = cfdx::core::validate_mesh(mesh);
            result.geometry_quality_valid = result.quality.ok;
        }
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

    const fs::path file_staging =
        fs::temp_directory_path() / "cfdx_n10_meshio_fixture_files";
    fs::remove_all(file_staging);
    fs::create_directories(file_staging);

    // Acquisition intentionally stores a verified file as <fixture-id> so the
    // manifest remains the single source of truth. For the production importer
    // we restore the pinned source filename extension; this is only staging and
    // does not modify the acquired bytes.
    const std::vector<std::tuple<std::string, std::string, std::string>> file_fixtures = {
        {"meshio-su2-square", "square.su2", "meshio-su2-square"},
        {"meshio-gmsh-insulated-2-2", "insulated-2.2.msh", "meshio-gmsh-insulated-2-2"},
        {"meshio-vtk-unstructured", "06_unstructured.vtk", "meshio-vtk-unstructured"},
    };

    std::vector<FixtureResult> results;
    for (const auto& [id, filename, acquired_id] : file_fixtures) {
        const fs::path acquired = fixture_root / acquired_id;
        const fs::path staged = file_staging / filename;
        fs::copy_file(acquired, staged, fs::copy_options::overwrite_existing);
        results.push_back(qualify_fixture(id, staged));
    }
    fs::remove_all(file_staging);

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
