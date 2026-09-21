#pragma once

#include "cfdx/thermodynamics/thermo_state.h"
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cfdx::thermodynamics {

struct ThermoCache {
    std::vector<ThermoState> state;
    std::vector<double> pressure;
    std::vector<double> temperature;
    bool valid = false;

    void resize(std::size_t n)
    {
        state.resize(n);
        pressure.resize(n);
        temperature.resize(n);
        valid = false;
    }

    void invalidate() noexcept { valid = false; }

    template<class Model>
    void update(const Model& model,
                const std::vector<double>& p,
                const std::vector<double>& T)
    {
        if (p.size() != T.size())
            throw std::invalid_argument("thermo cache size mismatch");

        // Compute into temporaries so a failed state evaluation cannot leave a
        // partially updated cache marked as valid.
        std::vector<ThermoState> new_state;
        std::vector<double> new_pressure = p;
        std::vector<double> new_temperature = T;
        new_state.reserve(p.size());

        for (std::size_t i = 0; i < p.size(); ++i)
            new_state.push_back(model.state(p[i], T[i]));

        state.swap(new_state);
        pressure.swap(new_pressure);
        temperature.swap(new_temperature);
        valid = true;
    }
};

} // namespace cfdx::thermodynamics
