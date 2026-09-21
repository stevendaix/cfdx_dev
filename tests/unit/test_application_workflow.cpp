#include "cfdx/application/checkpoint_store.h"
#include "cfdx/application/command_bus.h"
#include "cfdx/application/solver_adapter.h"
#include "cfdx/application/workbench.h"

#include <cassert>
#include <filesystem>
#include <string>

using namespace cfdx::application;

int main() {
    CaseModel model;
    model.name = "complete-workflow";
    model.mode = SimulationMode::Transient;
    model.time.start = 0.0;
    model.time.end = 2.0;
    model.time.step = 0.5;

    Workbench workbench(model);

    // Case tree + common selection model.
    assert(workbench.tree().selected().has_value());
    assert(workbench.tree().select(2));
    workbench.scene().add({10, SceneObjectType::Patch, "inlet", true});
    workbench.scene().add({11, SceneObjectType::Patch, "outlet", true});
    assert(workbench.select_scene_object(11));
    assert(workbench.selected_scene_object().value() == 11);

    // Configure before running.
    assert(workbench.properties().set(
        "physics.density",
        {1.0, "kg/m3", ChangeImpact::Hot}) == ChangeImpact::Hot);

    workbench.controller().validate();
    assert(workbench.controller().state() == SimulationState::Ready);
    assert(workbench.toolbar().can_run);

    // Solver adapter is invoked by the controller, not by the UI.
    std::size_t solver_iterations = 0;
    CallbackSolverAdapter solver;
    solver.validate_fn = [](const CaseModel&) {};
    solver.begin_fn = [](const CaseModel&, const Checkpoint&) {};
    solver.iterate_fn = [&](std::size_t, double) { ++solver_iterations; return true; };
    solver.end_fn = [] {};
    workbench.controller().set_solver_adapter(&solver);

    // Transient run-to-time case.
    workbench.controller().run({0, 1.0, false});
    assert(workbench.controller().state() == SimulationState::Converged);
    assert(workbench.controller().checkpoint().time >= 1.0);
    assert(solver_iterations == 2);

    // Restart/rebuild edit is explicitly visible.
    assert(workbench.properties().set(
        "mesh.size", {std::int64_t{32}, "", ChangeImpact::RequiresRebuild})
        == ChangeImpact::RequiresRebuild);
    assert(workbench.toolbar().restart_required);
    workbench.controller().acknowledge_restart();
    assert(!workbench.toolbar().restart_required);

    // Command/TUI workflow.
    CommandBus bus(workbench);
    assert(bus.execute("set output.frequency 20 hot").ok);
    assert(bus.execute("checkpoint").ok);
    assert(bus.execute("restore").ok);
    assert(bus.restored_snapshot().has_value());

    // Persistent checkpoint round-trip.
    const auto path = std::filesystem::temp_directory_path() / "cfdx_application_checkpoint.txt";
    TextCheckpointStore disk(path.string());
    disk.save({workbench.controller().model(), workbench.controller().latest_checkpoint()});
    assert(disk.contains());
    const auto restored = disk.load();
    assert(restored.model.name == model.name);
    assert(restored.checkpoint.iteration == workbench.controller().checkpoint().iteration);
    std::filesystem::remove(path);

    // Solver adapter contract is UI-independent.
    std::size_t iterations = 0;
    CallbackSolverAdapter adapter;
    adapter.validate_fn = [](const CaseModel&) {};
    adapter.begin_fn = [](const CaseModel&, const Checkpoint&) {};
    adapter.iterate_fn = [&](std::size_t, double) { ++iterations; return true; };
    adapter.end_fn = [] {};
    adapter.validate(workbench.controller().model());
    adapter.begin(workbench.controller().model(), workbench.controller().checkpoint());
    assert(adapter.iterate(1, 1.0));
    adapter.end();
    assert(iterations == 1);

    // Steady workflow remains supported.
    CaseModel steady;
    steady.name = "steady";
    steady.mode = SimulationMode::Steady;
    SimulationController steady_controller(steady);
    steady_controller.validate();
    steady_controller.run({5, -1.0, false});
    assert(steady_controller.state() == SimulationState::Paused);
    assert(steady_controller.checkpoint().iteration == 5);

    return 0;
}
