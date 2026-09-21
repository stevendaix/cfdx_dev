#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>

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

struct RuntimeDecision {
    ExecutionPolicy requested = ExecutionPolicy::AUTO;
    ExecutionPolicy selected = ExecutionPolicy::CPU;
    bool requires_gpu = false;
    bool uses_out_of_core = false;
    std::string reason;
};

struct RuntimeCapabilities {
    bool cpu = true;
    bool cuda = false;
    bool mpi = false;
    std::size_t gpu_memory_bytes = 0;
};

struct RuntimeWorkload {
    std::size_t estimated_bytes = 0;
    std::size_t safety_margin_bytes = 0;
    std::size_t runtime_reserve_bytes = 0;
};

inline RuntimeDecision choose_execution_policy(
    ExecutionPolicy requested,
    const RuntimeCapabilities& caps,
    const RuntimeWorkload& workload)
{
    const std::size_t reserve = workload.safety_margin_bytes + workload.runtime_reserve_bytes;
    const bool fits = workload.estimated_bytes <=
                      (caps.gpu_memory_bytes > reserve ? caps.gpu_memory_bytes - reserve : 0);

    RuntimeDecision d;
    d.requested = requested;

    if (requested == ExecutionPolicy::CPU) {
        d.selected = ExecutionPolicy::CPU;
        d.reason = "explicit CPU policy";
        return d;
    }
    if (requested == ExecutionPolicy::GPU) {
        if (!caps.cuda) throw std::runtime_error("GPU policy requested but CUDA is unavailable");
        if (!fits) throw std::runtime_error("GPU policy requested but workload exceeds GPU memory budget");
        d.selected = ExecutionPolicy::GPU;
        d.requires_gpu = true;
        d.reason = "explicit GPU policy and workload fits";
        return d;
    }
    if (requested == ExecutionPolicy::GPU_OUT_OF_CORE) {
        if (!caps.cuda) throw std::runtime_error("GPU OOC policy requested but CUDA is unavailable");
        d.selected = ExecutionPolicy::GPU_OUT_OF_CORE;
        d.requires_gpu = true;
        d.uses_out_of_core = true;
        d.reason = fits ? "explicit GPU OOC policy" : "workload exceeds full-GPU budget";
        return d;
    }

    // AUTO is deliberately a pre-run decision: there is no silent backend
    // migration during an iteration.
    if (caps.cuda && fits) {
        d.selected = ExecutionPolicy::GPU;
        d.requires_gpu = true;
        d.reason = "AUTO: workload fits the available GPU budget";
    } else if (caps.cuda) {
        d.selected = ExecutionPolicy::GPU_OUT_OF_CORE;
        d.requires_gpu = true;
        d.uses_out_of_core = true;
        d.reason = "AUTO: workload exceeds full-GPU budget";
    } else {
        d.selected = ExecutionPolicy::CPU;
        d.reason = "AUTO: no CUDA backend available";
    }
    return d;
}

} // namespace cfdx::runtime
