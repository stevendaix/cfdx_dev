#pragma once

#include <cstddef>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace cfdx::application {

struct MonitorSample {
    std::size_t iteration{0};
    double time{0.0};
    double value{0.0};
};

class Monitor {
public:
    explicit Monitor(std::string name = {}) : name_(std::move(name)) {}

    const std::string& name() const noexcept { return name_; }
    void add(std::size_t iteration, double time, double value) {
        samples_.push_back({iteration, time, value});
    }
    const std::vector<MonitorSample>& samples() const noexcept { return samples_; }

private:
    std::string name_;
    std::vector<MonitorSample> samples_;
};

class MonitorManager {
public:
    Monitor& create(std::string name) {
        return monitors_.try_emplace(name, Monitor{name}).first->second;
    }
    Monitor* find(const std::string& name) {
        const auto it = monitors_.find(name);
        return it == monitors_.end() ? nullptr : &it->second;
    }
    const std::map<std::string, Monitor>& all() const noexcept { return monitors_; }

private:
    std::map<std::string, Monitor> monitors_;
};

} // namespace cfdx::application
