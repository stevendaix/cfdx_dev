#include "cfdx/application/simulation_controller.h"
#include "common/test_harness.h"

#include <atomic>
#include <condition_variable>
#include <mutex>
#include <thread>

using namespace cfdx::application;
using namespace cfdx::testing;

int main()
{
    run_case("snapshot_matches_iteration_callback_state", [] {
        CaseModel model;
        model.mode = SimulationMode::Transient;
        model.time.start = 0.0;
        model.time.step = 0.25;
        model.time.max_iterations_per_step = 10;
        model.time.end = 10.0;

        SimulationController controller(model);
        std::mutex mutex;
        std::condition_variable cv;
        bool entered = false;
        bool release = false;
        bool coherent = true;

        controller.set_iteration_callback([&](std::size_t iteration, double time) {
            const auto snapshot = controller.snapshot();
            coherent = coherent &&
                       snapshot.checkpoint.iteration == iteration &&
                       snapshot.checkpoint.time == time &&
                       snapshot.model.revision == snapshot.checkpoint.case_revision;
            {
                std::lock_guard<std::mutex> lock(mutex);
                entered = true;
            }
            cv.notify_one();
            std::unique_lock<std::mutex> lock(mutex);
            cv.wait(lock, [&] { return release; });
            return true;
        });

        std::thread runner([&] { controller.run({3, -1.0, false}); });

        {
            std::unique_lock<std::mutex> lock(mutex);
            cv.wait(lock, [&] { return entered; });
        }

        EXPECT_THROW(
            controller.edit("pressure", Parameter{1.0, "Pa", ChangeImpact::RequiresRestart}),
            std::logic_error);

        {
            std::lock_guard<std::mutex> lock(mutex);
            release = true;
        }
        cv.notify_one();
        controller.stop();
        runner.join();

        EXPECT_TRUE(coherent);
        EXPECT_TRUE(controller.state() == SimulationState::Stopped ||
                    controller.state() == SimulationState::Paused);
    });

    run_case("restart_edit_and_checkpoint_revisions_are_coherent", [] {
        SimulationController controller;
        EXPECT_TRUE(controller.snapshot().checkpoint.case_revision == 0);

        const ChangeImpact impact = controller.edit(
            "mesh", Parameter{std::int64_t{2}, "", ChangeImpact::RequiresRebuild});
        EXPECT_TRUE(impact == ChangeImpact::RequiresRebuild);
        EXPECT_TRUE(controller.requires_restart());

        const auto edited = controller.snapshot();
        EXPECT_TRUE(edited.checkpoint.case_revision == edited.model.revision);
        EXPECT_TRUE(edited.model.revision == 1);

        controller.acknowledge_restart();
        const auto acknowledged = controller.snapshot();
        EXPECT_TRUE(!controller.requires_restart());
        EXPECT_TRUE(acknowledged.checkpoint.numerics_revision ==
                    acknowledged.model.numerics_revision);
    });

    return run_all();
}
