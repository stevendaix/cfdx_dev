#pragma once
#include "cfdx/core/linalg/vector.h"
#include <cstddef>
#include <vector>

namespace cfdx::core {

struct KrylovWorkspace {
    Vector r, r_hat, p, v, z, s, z_s, t, ap;
    void resize(std::size_t n) {
        r.resize(n); r_hat.resize(n); p.resize(n); v.resize(n);
        z.resize(n); s.resize(n); z_s.resize(n); t.resize(n); ap.resize(n);
    }
};

struct GmresWorkspace {
    std::size_t n = 0;
    std::size_t restart = 0;
    Vector r, w, z, vin, vout, ax;
    std::vector<double> basis, zbasis, h, cs, sn, g, y;

    void resize(std::size_t n_, std::size_t m_) {
        n = n_;
        restart = m_;
        r.resize(n); w.resize(n); z.resize(n); vin.resize(n); vout.resize(n); ax.resize(n);
        basis.assign((m_ + 1) * n_, 0.0);
        zbasis.assign(m_ * n_, 0.0);
        h.assign((m_ + 1) * m_, 0.0);
        cs.assign(m_, 0.0);
        sn.assign(m_, 0.0);
        g.assign(m_ + 1, 0.0);
        y.assign(m_, 0.0);
    }

    double* v(std::size_t j) { return basis.data() + j * n; }
    const double* v(std::size_t j) const { return basis.data() + j * n; }
    double* zv(std::size_t j) { return zbasis.data() + j * n; }
    const double* zv(std::size_t j) const { return zbasis.data() + j * n; }
    double& H(std::size_t i, std::size_t j) { return h[i * restart + j]; }
    const double& H(std::size_t i, std::size_t j) const { return h[i * restart + j]; }
};

} // namespace cfdx::core
