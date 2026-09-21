#pragma once

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <mutex>
#include <string>

namespace cfdx::utils {

class SimulationLogger {
public:
    explicit SimulationLogger(const std::string& path = "simulation.log")
        : log_file_(path, std::ios::app) {}

    void info(const std::string& msg) { write("INFO", msg); }
    void debug(const std::string& msg) { write("DEBUG", msg); }
    void error(const std::string& msg) { write("ERROR", msg); }

private:
    void write(const char* level, const std::string& msg) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (log_file_) log_file_ << "[" << level << "] " << msg << '\n';
    }
    std::ofstream log_file_;
    std::mutex mutex_;
};

class ProfilingHook {
public:
    using Clock = std::chrono::steady_clock;

    template<class Function>
    static auto profile(Function&& function) {
        const auto start = Clock::now();
        if constexpr (std::is_void_v<std::invoke_result_t<Function>>) {
            std::forward<Function>(function)();
            return std::chrono::duration_cast<std::chrono::microseconds>(
                Clock::now() - start);
        } else {
            auto result = std::forward<Function>(function)();
            return std::pair<decltype(result), std::chrono::microseconds>{
                std::move(result),
                std::chrono::duration_cast<std::chrono::microseconds>(
                    Clock::now() - start)};
        }
    }
};

inline SimulationLogger& simulation_logger() {
    static SimulationLogger logger;
    return logger;
}

inline void initialize_logging_and_profiling() {
    (void)simulation_logger();
}

inline void log_residual(const std::string& residual, double tolerance = 1e-6) {
    const double value = std::stod(residual);
    if (!std::isfinite(value))
        throw std::invalid_argument("log_residual: residual must be finite");
    if (std::abs(value) > tolerance)
        simulation_logger().error(
            "Residual too large: " + residual + " > " + std::to_string(tolerance));
    else
        simulation_logger().info(
            "Residual within tolerance: " + residual + " <= " + std::to_string(tolerance));
}

}  // namespace cfdx::utils
