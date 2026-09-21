#pragma once

#include "cfdx/application/case_model.h"
#include "cfdx/application/monitor.h"
#include "cfdx/application/post.h"
#include "cfdx/application/scene.h"
#include "cfdx/application/simulation_controller.h"

#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::application {

enum class TreeNodeType { Case, Geometry, Mesh, Physics, Materials, Boundaries, Numerics, Solver, Monitors, Reports };

struct TreeNode {
    std::uint64_t id{0};
    TreeNodeType type{TreeNodeType::Case};
    std::uint64_t parent{0};
    std::string name;
    bool expanded{true};
};

class CaseTree {
public:
    explicit CaseTree(const CaseModel& model) {
        reset(model);
    }

    void reset(const CaseModel& model) {
        nodes_.clear();
        next_id_ = 1;
        add(TreeNodeType::Case, 0, model.name);
        const auto root = nodes_.begin()->first;
        add(TreeNodeType::Geometry, root, "Geometry");
        add(TreeNodeType::Mesh, root, "Mesh");
        add(TreeNodeType::Physics, root, "Physics");
        add(TreeNodeType::Materials, root, "Materials");
        add(TreeNodeType::Boundaries, root, "Boundaries");
        add(TreeNodeType::Numerics, root, "Numerics");
        add(TreeNodeType::Solver, root, "Solver");
        add(TreeNodeType::Monitors, root, "Monitors");
        add(TreeNodeType::Reports, root, "Reports");
        selected_ = root;
    }

    std::uint64_t add(TreeNodeType type, std::uint64_t parent, std::string name) {
        const auto id = next_id_++;
        nodes_.emplace(id, TreeNode{id, type, parent, std::move(name), true});
        return id;
    }

    bool select(std::uint64_t id) {
        if (!nodes_.count(id)) return false;
        selected_ = id;
        return true;
    }

    std::optional<TreeNode> selected() const {
        if (!selected_) return std::nullopt;
        return nodes_.at(*selected_);
    }

    const std::map<std::uint64_t, TreeNode>& nodes() const noexcept { return nodes_; }

private:
    std::map<std::uint64_t, TreeNode> nodes_;
    std::uint64_t next_id_{1};
    std::optional<std::uint64_t> selected_;
};

class PropertyEditor {
public:
    explicit PropertyEditor(SimulationController& controller) : controller_(controller) {}

    std::optional<Parameter> get(const std::string& key) const {
        return controller_.model().get_parameter(key);
    }

    ChangeImpact set(const std::string& key, Parameter value) {
        return controller_.edit(key, std::move(value));
    }

private:
    SimulationController& controller_;
};

struct RunToolbarState {
    bool can_run{false};
    bool can_pause{false};
    bool can_stop{false};
    bool can_edit{false};
    bool restart_required{false};
};

class Workbench {
public:
    explicit Workbench(CaseModel model = {})
        : controller_(std::move(model)), tree_(controller_.model()), properties_(controller_) {}

    SimulationController& controller() noexcept { return controller_; }
    const SimulationController& controller() const noexcept { return controller_; }
    CaseTree& tree() noexcept { return tree_; }
    const CaseTree& tree() const noexcept { return tree_; }
    PropertyEditor& properties() noexcept { return properties_; }
    SceneModel& scene() noexcept { return scene_; }
    MonitorManager& monitors() noexcept { return monitors_; }
    PostProcessor& post() noexcept { return post_; }

    bool select_scene_object(std::size_t id) {
        if (!scene_.select(id)) return false;
        selected_scene_object_ = id;
        return true;
    }

    std::optional<std::size_t> selected_scene_object() const noexcept {
        return selected_scene_object_;
    }

    RunToolbarState toolbar() const {
        const auto state = controller_.state();
        const bool idle = state == SimulationState::Ready ||
                          state == SimulationState::Paused ||
                          state == SimulationState::Stopped;
        return {
            idle,
            state == SimulationState::Running,
            state == SimulationState::Running || state == SimulationState::Paused,
            idle,
            controller_.requires_restart()
        };
    }

private:
    SimulationController controller_;
    CaseTree tree_;
    PropertyEditor properties_;
    SceneModel scene_;
    MonitorManager monitors_;
    PostProcessor post_;
    std::optional<std::size_t> selected_scene_object_;
};

} // namespace cfdx::application
