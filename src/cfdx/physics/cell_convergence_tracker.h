#pragma once

// Opt-in, observational per-cell residual trace. The values are reassembled
// nonlinear momentum A(U,p)*U-b(U,p), not the pressure Krylov residual. They
// must NEVER be used as a solver acceptance criterion. Samples are emitted
// only after the outer nonlinear transaction has committed.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#include <memory>
#include "cell_field_map_writer.h"

namespace cfdx::physics {

inline std::string cell_trace_json_string(const std::string& input) {
    constexpr char hex[] = "0123456789abcdef";
    std::string result = "\"";
    for (char byte : input) {
        const auto c = static_cast<unsigned char>(byte);
        switch (c) {
            case '"': result += "\\\""; break;
            case '\\': result += "\\\\"; break;
            case '\n': result += "\\n"; break;
            case '\r': result += "\\r"; break;
            case '\t': result += "\\t"; break;
            default:
                if (c < 0x20) {
                    result += "\\u00";
                    result += hex[c >> 4];
                    result += hex[c & 0xf];
                } else result += static_cast<char>(c);
        }
    }
    return result + '"';
}

struct CellTracePoint {
    std::size_t cell = 0;
    std::array<double, 3> center{};
    std::array<double, 3> momentum_abs{};
    double momentum_scaled = 0.0;
    double mass_imbalance_abs = 0.0;
    std::uint16_t nondecrease_streak = 0;
};

// Directly usable by small deterministic C++ unit tests, independently of a
// mesh or CFD equation. All cell IDs refer to their original input order.
class CellConvergenceTracker {
public:
    CellConvergenceTracker(
        const std::string& path, const std::string& mesh_source,
        std::size_t n_cells, std::size_t n_faces,
        std::size_t every = 5, std::size_t top_k = 24,
        double reference_relative_tolerance = 1e-6,
        const std::string& full_map_path = "")
        : previous_(n_cells, std::numeric_limits<double>::quiet_NaN()),
          streak_(n_cells, 0), every_(every), top_k_(top_k),
          reference_tolerance_(reference_relative_tolerance)
    {
        if (path.empty() || mesh_source.empty() || n_cells == 0 ||
            n_cells > 5'000'000 || n_faces == 0 || every == 0 || every > 100'000 ||
            top_k == 0 || top_k > 64 ||
            !std::isfinite(reference_tolerance_) || reference_tolerance_ <= 0.0)
            throw std::invalid_argument("invalid cell convergence diagnostic configuration or memory budget");
        if (std::filesystem::is_symlink(path))
            throw std::invalid_argument("cell convergence trace must not follow a symlink");
        if (!full_map_path.empty())
            full_map_ = std::make_unique<CellFieldMapWriter>(
                full_map_path, mesh_source, n_cells, n_faces, reference_tolerance_);
        stream_.open(path, std::ios::out | std::ios::trunc);
        if (!stream_) throw std::runtime_error("cannot open cell convergence diagnostic stream");
        stream_ << std::setprecision(17)
                << "{\"format\":\"CFDX-CELL-CONVERGENCE\",\"schema_version\":1,"
                << "\"kind\":\"header\",\"mesh_source\":" << cell_trace_json_string(mesh_source)
                << ",\"n_cells\":" << n_cells << ",\"n_faces\":" << n_faces
                << ",\"every\":" << every_ << ",\"top_k\":" << top_k_
                << ",\"reference_relative_tolerance\":" << reference_tolerance_
                << ",\"momentum_scale\":\"FIXED per-component global max abs RHS from first accepted equation assembly\""
                << ",\"stagnation_rule\":\"at least three consecutive accepted transitions with score >= 0.98 times previous and score > reported reference tolerance\""
                << ",\"coverage\":\"bounded top-K union: worst momentum, worst mass, persistent momentum; NOT all cells in each frame\""
                << ",\"convergence_certificate\":false}\n";
        verify_write();
    }

    // Call exactly once per *accepted* outer iteration. A change of nonlinear
    // operator (e.g. TVD continuation) invalidates trend comparability.
    void observe(
        std::size_t iteration,
        const std::vector<std::array<double, 3>>& momentum_abs,
        const std::vector<double>& mass_imbalance_abs,
        const std::array<double, 3>& rhs_scales,
        const std::vector<std::array<double, 3>>& centers,
        bool same_operator_as_previous = true,
        bool force_frame = false)
    {
        const auto count = previous_.size();
        if (iteration == 0 || iteration <= last_iteration_ ||
            momentum_abs.size() != count || mass_imbalance_abs.size() != count ||
            centers.size() != count)
            throw std::invalid_argument("cell convergence trace has inconsistent iteration or cell topology");
        for (double scale : rhs_scales)
            if (!std::isfinite(scale) || scale <= 0.0)
                throw std::invalid_argument("cell convergence requires finite positive RHS scales");
        // Freeze the denominator at the first accepted iteration: comparing
        // residual/RHS ratios with a different RHS at every nonlinear step can
        // manufacture an apparent plateau or improvement in one cell.
        if (!fixed_rhs_scales_valid_) {
            fixed_rhs_scales_ = rhs_scales;
            fixed_rhs_scales_valid_ = true;
        }
        std::vector<double> scores(count, std::numeric_limits<double>::quiet_NaN());
        std::size_t above_reference = 0;
        std::size_t persistent = 0;
        std::size_t nonfinite = 0;
        using Pair = std::pair<double, std::size_t>;
        using Heap = std::priority_queue<Pair, std::vector<Pair>, std::greater<Pair>>;
        Heap worst_momentum, worst_mass, worst_persistent;
        auto offer = [&](Heap& heap, const double priority, std::size_t c) {
            if (heap.size() < top_k_) heap.emplace(priority, c);
            else if (Pair{priority, c} > heap.top()) {
                heap.pop();
                heap.emplace(priority, c);
            }
        };
        for (std::size_t c = 0; c < count; ++c) {
            const auto& row = momentum_abs[c];
            double score = 0.0;
            bool finite = std::isfinite(mass_imbalance_abs[c]) &&
                          mass_imbalance_abs[c] >= 0.0;
            for (std::size_t d = 0; d < 3; ++d) {
                finite = finite && std::isfinite(row[d]) && row[d] >= 0.0 &&
                         std::isfinite(centers[c][d]);
                if (finite) score = std::max(score, row[d] / fixed_rhs_scales_[d]);
            }
            if (!finite || !std::isfinite(score)) {
                ++nonfinite;
                streak_[c] = 0;
                previous_[c] = std::numeric_limits<double>::quiet_NaN();
                continue;
            }
            scores[c] = score;
            if (score > reference_tolerance_) {
                ++above_reference;
                const double prior = previous_[c];
                if (same_operator_as_previous && std::isfinite(prior) &&
                    prior > reference_tolerance_ && score >= 0.98 * prior)
                    streak_[c] = static_cast<std::uint16_t>(std::min<unsigned>(
                        65535, static_cast<unsigned>(streak_[c]) + 1));
                else streak_[c] = 0;
            } else streak_[c] = 0;
            previous_[c] = score;
            offer(worst_momentum, score, c);
            offer(worst_mass, mass_imbalance_abs[c], c);
            if (streak_[c] >= 3) {
                ++persistent;
                offer(worst_persistent, score, c);
            }
        }
        last_iteration_ = iteration;
        if (iteration != 1 && iteration % every_ != 0 && !force_frame)
            return;
        std::vector<std::size_t> chosen;
        auto take = [&](Heap& heap) {
            while (!heap.empty()) {
                chosen.push_back(heap.top().second);
                heap.pop();
            }
        };
        take(worst_momentum);
        take(worst_mass);
        take(worst_persistent);
        std::sort(chosen.begin(), chosen.end());
        chosen.erase(std::unique(chosen.begin(), chosen.end()), chosen.end());
        stream_ << "{\"kind\":\"frame\",\"iteration\":" << iteration
                << ",\"assessed\":" << (count - nonfinite)
                << ",\"nonfinite\":" << nonfinite
                << ",\"momentum_above_reference\":" << above_reference
                << ",\"persistent_momentum_cells\":" << persistent
                << ",\"operator_comparable\":" << (same_operator_as_previous ? "true" : "false")
                << ",\"rhs_scale\":[" << fixed_rhs_scales_[0] << ',' << fixed_rhs_scales_[1] << ',' << fixed_rhs_scales_[2]
                << "],\"cells\":[";
        for (std::size_t i = 0; i < chosen.size(); ++i) {
            const auto c = chosen[i];
            if (i) stream_ << ',';
            stream_ << "{\"cell_id\":" << c << ",\"center\":["
                    << centers[c][0] << ',' << centers[c][1] << ',' << centers[c][2]
                    << "],\"momentum_abs\":["
                    << momentum_abs[c][0] << ',' << momentum_abs[c][1] << ',' << momentum_abs[c][2]
                    << "],\"momentum_scaled\":" << scores[c]
                    << ",\"mass_imbalance_abs\":" << mass_imbalance_abs[c]
                    << ",\"nondecrease_streak\":" << streak_[c] << '}';
        }
        stream_ << "]}\n";
        verify_write();
        last_frame_iteration_ = iteration;
        if (full_map_)
            full_map_->write_frame(iteration, same_operator_as_previous, scores, mass_imbalance_abs, streak_);
    }

    void finish(std::size_t accepted_iterations) {
        if (accepted_iterations != last_iteration_ ||
            (accepted_iterations > 0 && last_frame_iteration_ == 0))
            throw std::logic_error("cell convergence trace end disagrees with accepted iterations");
        // An aborted run can have zero accepted iterations. Never invent a
        // sample or claim convergence; persist an explicit empty end marker.
        // A truncated run can end after the last periodic sample; that
        // sample must keep its actual iteration number.
        stream_ << "{\"kind\":\"end\",\"iterations\":" << accepted_iterations
                << ",\"last_frame_iteration\":" << last_frame_iteration_ << "}\n";
        verify_write();
        if (full_map_) full_map_->finish(accepted_iterations);
    }

private:
    void verify_write() {
        stream_.flush();
        if (!stream_) throw std::runtime_error("failed to persist cell convergence diagnostic stream");
    }
    std::unique_ptr<CellFieldMapWriter> full_map_;
    std::vector<double> previous_;
    std::vector<std::uint16_t> streak_;
    std::array<double, 3> fixed_rhs_scales_{1.0, 1.0, 1.0};
    bool fixed_rhs_scales_valid_ = false;
    std::size_t every_;
    std::size_t top_k_;
    double reference_tolerance_;
    std::ofstream stream_;
    std::size_t last_iteration_ = 0;
    std::size_t last_frame_iteration_ = 0;
};

}  // namespace cfdx::physics
