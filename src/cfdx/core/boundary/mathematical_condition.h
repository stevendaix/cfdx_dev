#pragma once
#include "cfdx/core/boundary/value_provider.h"
#include <memory>
#include <stdexcept>
#include <variant>
namespace cfdx::core {
struct Dirichlet { std::shared_ptr<const ValueProvider> value; };
struct Neumann { std::shared_ptr<const ValueProvider> gradient; };
struct Robin { double alpha = 1.0; double beta = 0.0; std::shared_ptr<const ValueProvider> gamma; };
struct Flux { std::shared_ptr<const ValueProvider> value; };
struct Mixed { std::shared_ptr<const ValueProvider> value; std::shared_ptr<const ValueProvider> weight; };
struct Coupled { std::string target_boundary; };
struct Periodic { std::string partner_boundary; };
using MathematicalCondition = std::variant<Dirichlet, Neumann, Robin, Flux, Mixed, Coupled, Periodic>;
inline void validate_condition(const MathematicalCondition& condition) {
    std::visit([](const auto& c) {
        using T = std::decay_t<decltype(c)>;
        if constexpr (std::is_same_v<T, Dirichlet> || std::is_same_v<T, Neumann> ||
                      std::is_same_v<T, Flux>) {
            if (!c.value && !c.gradient) throw std::invalid_argument("boundary condition requires a value provider");
        } else if constexpr (std::is_same_v<T, Robin>) {
            if (!c.gamma) throw std::invalid_argument("Robin condition requires gamma");
            if (c.alpha == 0.0 && c.beta == 0.0) throw std::invalid_argument("Robin alpha and beta cannot both be zero");
        } else if constexpr (std::is_same_v<T, Mixed>) {
            if (!c.value || !c.weight) throw std::invalid_argument("Mixed condition requires value and weight providers");
        } else if constexpr (std::is_same_v<T, Coupled>) {
            if (c.target_boundary.empty()) throw std::invalid_argument("Coupled condition requires a target boundary");
        } else if constexpr (std::is_same_v<T, Periodic>) {
            if (c.partner_boundary.empty()) throw std::invalid_argument("Periodic condition requires a partner boundary");
        }
    }, condition);
}
}