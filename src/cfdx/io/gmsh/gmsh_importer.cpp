#include "gmsh_importer.h"
#include "cfdx/io/mesh/mesh_importer.h"
#include "cfdx/io/hdf5/hdf5_reader.h"
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <sstream>

namespace cfdx::io::gmsh {
namespace {
std::string quote(const std::string& value)
{
    std::string out="'";
    for(char c:value) {
        if(c=='\'') out+="'\\''";
        else out+=c;
    }
    return out+"'";
}
}

bool import_gmsh_mesh(const std::string& path, cfdx::core::Mesh& mesh)
{
    return cfdx::io::mesh::import_mesh(path,mesh);
}

bool import_gmsh_scalar(const std::string& path, const std::string& field_name,
                        cfdx::core::ScalarCellField& field)
{
#ifndef CFDX_SOURCE_DIR
    (void)path; (void)field_name; (void)field;
    return false;
#else
    namespace fs=std::filesystem;
    if(!fs::is_regular_file(path) || field_name.empty()) return false;
    const auto stamp=std::chrono::steady_clock::now().time_since_epoch().count();
    const fs::path output=fs::temp_directory_path()/
        ("cfdx_gmsh_field_"+std::to_string(stamp)+".h5");
    const fs::path script=fs::path(CFDX_SOURCE_DIR)/"scripts"/"gmsh_scalar_import.py";

    std::ostringstream command;
    command<<
#ifdef CFDX_PYTHON_EXECUTABLE
           quote(CFDX_PYTHON_EXECUTABLE)
#else
           "python3"
#endif
           <<" "<<quote(script.string())<<" "<<quote(path)<<" "
           <<quote(output.string())<<" "<<quote(field_name);
    const int rc=std::system(command.str().c_str());
    if(rc!=0 || !fs::exists(output)) {
        std::error_code ec; fs::remove(output,ec);
        return false;
    }

    const bool ok=cfdx::io::read_field_hdf5(output.string(),field);
    std::error_code ec; fs::remove(output,ec);
    return ok && field.dimension()==1;
#endif
}

} // namespace cfdx::io::gmsh
