#pragma once

#include "cfdx/application/monitor.h"

#include <functional>
#include <string>
#include <vector>

namespace cfdx::application {

struct DerivedField {
    std::string name;
    std::string expression;
};

class PostProcessor {
public:
    void add_field(std::string name, std::string expression) {
        fields_.push_back({std::move(name), std::move(expression)});
    }
    void add_report(std::string name, std::function<double()> evaluator) {
        reports_.push_back({std::move(name), std::move(evaluator)});
    }
    std::vector<std::pair<std::string, double>> evaluate_reports() const {
        std::vector<std::pair<std::string, double>> result;
        for (const auto& report : reports_) result.emplace_back(report.first, report.second());
        return result;
    }
    const std::vector<DerivedField>& fields() const noexcept { return fields_; }

private:
    std::vector<DerivedField> fields_;
    std::vector<std::pair<std::string, std::function<double()>>> reports_;
};

} // namespace cfdx::application
