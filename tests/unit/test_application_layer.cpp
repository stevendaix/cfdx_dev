#include "cfdx/application/case_model.h"
#include "cfdx/application/monitor.h"
#include "cfdx/application/post.h"
#include "cfdx/application/scene.h"
#include "cfdx/application/simulation_controller.h"
#include "cfdx/application/tui.h"

#include <cassert>
#include <cmath>
#include <string>

using namespace cfdx::application;

int main() {
    CaseModel model;
    model.name = "fluent-like-demo";
    model.mode = SimulationMode::Transient;
    model.time.start = 0.0;
    model.time.end = 1.0;
    model.time.step = 0.1;
    model.time.max_iterations_per_step = 10;

    SimulationController controller(model);
    controller.validate();
    assert(controller.state() == SimulationState::Ready);

    std::size_t callbacks = 0;
    controller.set_iteration_callback([&](std::size_t iteration, double time) {
        ++callbacks;
        assert(iteration >= 1);
        assert(time >= 0.0);
        if (iteration == 3) controller.pause();
        return true;
    });
    controller.run({10, -1.0, false});
    assert(controller.state() == SimulationState::Paused);
    assert(controller.checkpoint().iteration == 3);

    const auto hot = controller.edit("output.frequency", {std::int64_t{10}, "", ChangeImpact::Hot});
    assert(hot == ChangeImpact::Hot);
    const auto restart = controller.edit("numerics.pressure_solver", {std::string{"gmres"}, "", ChangeImpact::RequiresRebuild});
    assert(restart == ChangeImpact::RequiresRebuild);
    assert(controller.requires_restart());
    controller.acknowledge_restart();
    assert(!controller.requires_restart());

    controller.set_iteration_callback([&](std::size_t, double) {
        controller.stop();
        return true;
    });
    controller.run({10, -1.0, false});
    assert(controller.state() == SimulationState::Stopped);

    const std::string tui = TuiRenderer::render(controller);
    assert(tui.find("CFDX") != std::string::npos);
    assert(tui.find("STOPPED") != std::string::npos);
    assert(tui.find("[RUN]") != std::string::npos);

    MonitorManager monitors;
    auto& residual = monitors.create("continuity");
    residual.add(1, 0.1, 1e-2);
    residual.add(2, 0.2, 1e-4);
    assert(monitors.find("continuity") != nullptr);
    assert(monitors.find("continuity")->samples().size() == 2);

    SceneModel scene;
    scene.add({1, SceneObjectType::Patch, "inlet", true});
    scene.add({2, SceneObjectType::Patch, "outlet", true});
    assert(scene.select(2));
    assert(scene.selected().value() == 2);
    assert(!scene.select(42));

    PostProcessor post;
    post.add_field("Mach", "mag(U)/a");
    post.add_report("mass_balance", [] { return 1.0e-12; });
    assert(post.fields().size() == 1);
    const auto reports = post.evaluate_reports();
    assert(reports.size() == 1);
    assert(std::abs(reports[0].second) < 1e-10);

    assert(callbacks == 3);
    return 0;
}
