#pragma once

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>

namespace cfdx::runtime {

enum class ExecutionPolicy : std::uint8_t {
    CPU = 0,
    GPU,
    GPU_OUT_OF_CORE,
    AUTO
};

inline const char* to_string(ExecutionPolicy p) {
    switch (p) {
        case ExecutionPolicy::CPU: return "CPU";
        case ExecutionPolicy::GPU: return "GPU";
        case ExecutionPolicy::GPU_OUT_OF_CORE: return "GPU_OUT_OF_CORE";
        case ExecutionPolicy::AUTO: return "AUTO";
        default: return "UNKNOWN";
    }
}

inline ExecutionPolicy execution_policy_from_string(const std::string& s) {
    if (s == "CPU" || s == "cpu") return ExecutionPolicy::CPU;
    if (s == "GPU" || s == "gpu") return ExecutionPolicy::GPU;
    if (s == "GPU_OUT_OF_CORE" || s == "gpu_out_of_core" || s == "gpu-ooc")
        return ExecutionPolicy::GPU_OUT_OF_CORE;
    if (s == "AUTO" || s == "auto") return ExecutionPolicy::AUTO;
    throw std::invalid_argument("unknown execution policy: " + s);
}

enum class CpuExecutionMode : std::uint8_t {
    SERIAL = 0,
    PARALLEL
};

struct CpuExecutionConfig {
    CpuExecutionMode mode = CpuExecutionMode::PARALLEL;
    std::size_t threads = 0;       // 0 = use the detected CPU capacity
    bool deterministic = false;    // reduction/order constraints are backend-owned
};

struct RuntimeDecision {
    ExecutionPolicy requested = ExecutionPolicy::AUTO;
    ExecutionPolicy selected = ExecutionPolicy::CPU;
    bool requires_gpu = false;
    bool uses_out_of_core = false;
    std::size_t cpu_threads = 1;
    bool deterministic = false;
    std::string reason;
};

struct RuntimeCapabilities {
    bool cpu = true;
    bool cuda = false;
    bool mpi = false;
    std::size_t gpu_memory_bytes = 0;
    std::size_t cpu_threads = 0;  // 0 = query std::thread::hardware_concurrency()
};

struct RuntimeWorkload {
    std::size_t estimated_bytes = 0;
    std::size_t safety_margin_bytes = 0;
    std::size_t runtime_reserve_bytes = 0;
};

inline std::size_t resolve_cpu_threads(
    const CpuExecutionConfig& config,
    const RuntimeCapabilities& caps)
{
    if (config.mode == CpuExecutionMode::SERIAL) {
        if (config.threads > 1) {
            throw std::invalid_argument("serial CPU execution requires threads <= 1");
        }
        return 1;
    }

    const std::size_t available =
        caps.cpu_threads != 0 ? caps.cpu_threads
                              : static_cast<std::size_t>(std::thread::hardware_concurrency());

    if (config.threads != 0) {
        if (available != 0 && config.threads > available) {
            throw std::invalid_argument("requested CPU threads exceed available CPU capacity");
        }
        return config.threads;
    }

    return available == 0 ? 1 : available;
}

inline RuntimeDecision choose_execution_policy(
    ExecutionPolicy requested,
    const RuntimeCapabilities& caps,
    const RuntimeWorkload& workload,
    const CpuExecutionConfig& cpu_config)
{
    if (!caps.cpu && requested == ExecutionPolicy::CPU) {
        throw std::runtime_error("CPU policy requested but CPU backend is unavailable");
    }

    const bool reserve_overflow =
        workload.safety_margin_bytes >
        std::numeric_limits<std::size_t>::max() - workload.runtime_reserve_bytes;
    const std::size_t reserve = reserve_overflow
        ? std::numeric_limits<std::size_t>::max()
        : workload.safety_margin_bytes + workload.runtime_reserve_bytes;

    const bool fits = !reserve_overflow &&
        reserve <= caps.gpu_memory_bytes &&
        workload.estimated_bytes <= caps.gpu_memory_bytes - reserve;

    RuntimeDecision d;
    d.requested = requested;

    auto select_cpu = [&] {
        if (!caps.cpu) {
            throw std::runtime_error(
                "CPU backend is unavailable and no GPU execution policy can be selected");
        }
        d.selected = ExecutionPolicy::CPU;
        d.cpu_threads = resolve_cpu_threads(cpu_config, caps);
        d.deterministic = cpu_config.deterministic;
    };

    if (requested == ExecutionPolicy::CPU) {
        select_cpu();
        d.reason = "explicit CPU policy";
        return d;
    }

    if (requested == ExecutionPolicy::GPU) {
        if (!caps.cuda) {
            throw std::runtime_error("GPU policy requested but CUDA is unavailable");
        }
        if (!fits) {
            throw std::runtime_error(
                "GPU policy requested but workload exceeds GPU memory budget");
        }
        d.selected = ExecutionPolicy::GPU;
        d.requires_gpu = true;
        d.deterministic = cpu_config.deterministic;
        d.reason = "explicit GPU policy and workload fits";
        return d;
    }

    if (requested == ExecutionPolicy::GPU_OUT_OF_CORE) {
        if (!caps.cuda) {
            throw std::runtime_error("GPU OOC policy requested but CUDA is unavailable");
        }
        d.selected = ExecutionPolicy::GPU_OUT_OF_CORE;
        d.requires_gpu = true;
        d.uses_out_of_core = true;
        d.deterministic = cpu_config.deterministic;
        d.reason = fits ? "explicit GPU OOC policy" : "workload exceeds full-GPU budget";
        return d;
    }

    // AUTO is deliberately a pre-run decision: there is no silent backend
    // migration during an iteration.
    if (caps.cuda && fits) {
        d.selected = ExecutionPolicy::GPU;
        d.requires_gpu = true;
        d.deterministic = cpu_config.deterministic;
        d.reason = "AUTO: workload fits the available GPU budget";
    } else if (caps.cuda) {
        d.selected = ExecutionPolicy::GPU_OUT_OF_CORE;
        d.requires_gpu = true;
        d.uses_out_of_core = true;
        d.deterministic = cpu_config.deterministic;
        d.reason = "AUTO: workload exceeds full-GPU budget";
    } else {
        select_cpu();
        d.reason = "AUTO: no CUDA backend available";
    }
    return d;
}

inline RuntimeDecision choose_execution_policy(
    ExecutionPolicy requested,
    const RuntimeCapabilities& caps,
    const RuntimeWorkload& workload)
{
    return choose_execution_policy(
        requested, caps, workload, CpuExecutionConfig{});
}

} // namespace cfdx::runtime
