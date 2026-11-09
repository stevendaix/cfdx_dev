#include "cfdx/physics/cell_convergence_tracker.h"
#include <array>
#ifdef NDEBUG
#undef NDEBUG
#endif
#include <cassert>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

int main() {
    using cfdx::physics::CellConvergenceTracker;
    const std::string output = "cell_convergence_tracker_regression.jsonl";
    const std::vector<std::array<double, 3>> positions{
        {0, 0, 0}, {1, -2, 3}, {4, 5, 6}, {7, 8, 9}};
    std::vector<std::array<double, 3>> momentum{
        {0.01, 0, 0}, {2, 0, 0}, {3, 0, 0}, {0.3, 0, 0}};
    std::vector<double> continuity{0.00001, 0.001, 0.02, 0.03};
    {
        CellConvergenceTracker tracer(output, "/tmp/mesh.h5", 4, 9, 2, 2, 0.0001);
        tracer.observe(1, momentum, continuity, {1, 1, 1}, positions, false);
        // One rejected trial is deliberately not observed: the unchanged
        // committed residual persists, and streaks count accepted states only.
        momentum[2][0] = 3.01;
        tracer.observe(2, momentum, continuity, {1, 1, 1}, positions);
        tracer.observe(3, momentum, continuity, {1, 1, 1}, positions);
        // A changing RHS must not change the per-cell trend denominator.
        tracer.observe(4, momentum, continuity, {10, 1, 1}, positions);
        tracer.observe(5, momentum, continuity, {1, 1, 1}, positions, false, true); // force final, operator changed
        tracer.finish(5);
    }
    std::ifstream in(output);
    assert(in.good());
    std::string line;
    std::vector<std::string> lines;
    while (std::getline(in, line)) lines.push_back(line);
    assert(lines.size() == 6); // header + frames 1/2/4/5 + end
    assert(lines.at(0).find("CFDX-CELL-CONVERGENCE") != std::string::npos);
    assert(lines.at(1).find("\"iteration\":1") != std::string::npos);
    assert(lines.at(2).find("\"iteration\":2") != std::string::npos);
    assert(lines.at(3).find("\"iteration\":4") != std::string::npos);
    assert(lines.at(3).find("\"nondecrease_streak\":3") != std::string::npos);
    assert(lines.at(4).find("\"iteration\":5") != std::string::npos);
    assert(lines.at(4).find("\"nondecrease_streak\":0") != std::string::npos);
    assert(lines.at(5).find("\"iterations\":5") != std::string::npos);
    assert(lines.at(2).find("\"cell_id\":2") != std::string::npos);
    bool failed = false;
    try {
        CellConvergenceTracker bad("", "mesh", 4, 8);
    } catch (const std::invalid_argument&) { failed = true; }
    assert(failed);
    failed = false;
    try {
        CellConvergenceTracker bad(output, "mesh", 5'000'001, 8);
    } catch (const std::invalid_argument&) { failed = true; }
    assert(failed);
    // Regression: a failed run can leave no committed iteration at all.
    {
        CellConvergenceTracker empty(output, "/tmp/mesh.h5", 4, 9, 3, 2);
        empty.finish(0);
    }
    {
        std::ifstream zero(output);
        assert(std::getline(zero, line));
        assert(std::getline(zero, line));
        assert(line.find("\"iterations\":0") != std::string::npos);
        assert(!std::getline(zero, line));
    }
    // Regression: stopping between scheduled snapshots must not fabricate a frame.
    {
        CellConvergenceTracker tracer(output, "/tmp/mesh.h5", 4, 9, 5, 2);
        tracer.observe(1, momentum, continuity, {1, 1, 1}, positions, false);
        tracer.observe(2, momentum, continuity, {1, 1, 1}, positions);
        tracer.observe(3, momentum, continuity, {1, 1, 1}, positions);
        tracer.finish(3);
    }
    {
        std::ifstream early(output);
        std::vector<std::string> rows;
        while (std::getline(early, line)) rows.push_back(line);
        assert(rows.size() == 3);
        assert(rows[2].find("\"iterations\":3") != std::string::npos);
        assert(rows[2].find("\"last_frame_iteration\":1") != std::string::npos);
    }
    std::remove(output.c_str());
}
