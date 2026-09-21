// M0.7-T06 — Source Terms (Chemistry Module)
//
// Architecture réorganisée P2:
//   chemistry/source_terms.h - Source term linearization (Su + Sp*phi)
//   Transport models moved to transport/*/ modules

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <cstddef>
#include <functional>
#include <vector>

namespace cfdx {
namespace chemistry {

using cfdx::core::Location;
using cfdx::core::Field;

struct SourceTerm {
    Field<double, Location::CELL> Su;
    Field<double, Location::CELL> Sp;

    SourceTerm() = default;

    SourceTerm(std::size_t n_cells, const std::string& name = "source")
        : Su(n_cells, name + "_Su", "kg/m3/s", 1),
          Sp(n_cells, name + "_Sp", "1/s", 1)
    {
        Su.fill(0.0);
        Sp.fill(0.0);
    }

    inline void apply_to_system(
        std::vector<double>& diag,
        std::vector<double>& rhs) const
    {
        const std::size_t n = Su.size();
        const double* su_data = Su.component_data(0);
        const double* sp_data = Sp.component_data(0);

        for (std::size_t c = 0; c < n; ++c) {
            diag[c] -= sp_data[c];
            rhs[c] += su_data[c];
        }
    }

    inline void apply_to_cell(std::size_t c, double& diag, double& rhs) const {
        const double* su_data = Su.component_data(0);
        const double* sp_data = Sp.component_data(0);
        diag -= sp_data[c];
        rhs += su_data[c];
    }
};

inline SourceTerm make_constant_source(std::size_t n_cells, double su_value, const std::string& name = "const_source") {
    SourceTerm st(n_cells, name);
    st.Su.fill(su_value);
    return st;
}

inline SourceTerm make_linear_source(const Field<double, Location::CELL>& phi, double k, const std::string& name = "linear_source") {
    SourceTerm st(phi.size(), name);
    double* sp_data = st.Sp.component_data(0);
    for (std::size_t c = 0; c < phi.size(); ++c) {
        sp_data[c] = k;
    }
    return st;
}

inline SourceTerm make_decay_source(const Field<double, Location::CELL>& phi, double lambda, const std::string& name = "decay_source") {
    return make_linear_source(phi, -lambda, name);
}

using LinearizationFunc = std::function<SourceTerm(const Field<double, Location::CELL>&, const cfdx::core::Mesh&)>;

}  // namespace chemistry
}  // namespace cfdx