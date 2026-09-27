#pragma once
#include "cfdx/core/boundary/boundary_role.h"
#include <cstddef>
#include <string>
#include <utility>
#include <vector>
namespace cfdx::core {
class Boundary {
public:
    Boundary() = default;
    Boundary(std::string name, BoundaryRole role, std::vector<std::size_t> face_ids)
        : name_(std::move(name)), role_(role), face_ids_(std::move(face_ids)) {}
    const std::string& name() const noexcept { return name_; }
    BoundaryRole role() const noexcept { return role_; }
    const std::vector<std::size_t>& face_ids() const noexcept { return face_ids_; }
private:
    std::string name_;
    BoundaryRole role_ = BoundaryRole::NONE;
    std::vector<std::size_t> face_ids_;
};
}