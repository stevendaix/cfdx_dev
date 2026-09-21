#pragma once

#include <cstddef>
#include <optional>
#include <string>
#include <vector>

namespace cfdx::application {

enum class SceneObjectType { Body, Face, Edge, Region, Patch, Cell };

struct SceneObject {
    std::size_t id{0};
    SceneObjectType type{SceneObjectType::Body};
    std::string name;
    bool visible{true};
};

class SceneModel {
public:
    std::size_t add(SceneObject object) {
        objects_.push_back(std::move(object));
        return objects_.back().id;
    }
    bool select(std::size_t id) {
        for (const auto& object : objects_) {
            if (object.id == id) {
                selected_ = id;
                return true;
            }
        }
        return false;
    }
    void clear_selection() noexcept { selected_.reset(); }
    std::optional<std::size_t> selected() const noexcept { return selected_; }
    const std::vector<SceneObject>& objects() const noexcept { return objects_; }

private:
    std::vector<SceneObject> objects_;
    std::optional<std::size_t> selected_;
};

} // namespace cfdx::application
