#pragma once
#include "cfdx/thermodynamics/thermo_state.h"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>

namespace cfdx::thermodynamics {

struct ThermoCache {
    std::vector<ThermoState> state;
    std::vector<double> pressure;
    std::vector<double> temperature;
    bool valid = false;

    void resize(std::size_t n) {
        state.resize(n);
        pressure.resize(n);
        temperature.resize(n);
        valid = false;
    }

    void invalidate() noexcept { valid = false; }

    template<class Model>
    void update(const Model& model,
                const std::vector<double>& p,
                const std::vector<double>& T) {
        if (p.size() != T.size()) throw std::invalid_argument("thermo cache size mismatch");
        resize(p.size());
        for (std::size_t i = 0; i < p.size(); ++i) {
            state[i] = model.state(p[i], T[i]);
            pressure[i] = p[i];
            temperature[i] = T[i];
        }
        valid = true;
    }
};

} // namespace cfdx::thermodynamics
