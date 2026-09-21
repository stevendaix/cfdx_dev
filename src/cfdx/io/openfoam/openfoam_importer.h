#pragma once
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/field/field.h"
#include <string>
#include <vector>
#include <cstdint>
#include <utility>

namespace cfdx::io::openfoam {
using ScalarCellField = cfdx::core::ScalarCellField;
using Vec3CellField = cfdx::core::Vec3CellField;

struct PatchDef {
    std::string name;
    int type;
    int n_elements;
    std::vector<std::uint32_t> face_ids;
    PatchDef(std::string n, int t, int ne, std::vector<std::uint32_t> ids)
        : name(std::move(n)), type(t), n_elements(ne), face_ids(std::move(ids)) {}
};

struct OpenFOAMFace { int n; std::vector<std::uint32_t> indices; };

bool import_openfoam_case(const std::string&, cfdx::core::Mesh&);
bool import_openfoam_field(const std::string&, const std::string&, ScalarCellField&);
bool import_openfoam_field(const std::string&, const std::string&, Vec3CellField&);
bool read_openfoam_mesh_meshio(const std::string&, std::vector<cfdx::core::Vec3>&,
                               std::vector<OpenFOAMFace>&, std::vector<PatchDef>&);
}
