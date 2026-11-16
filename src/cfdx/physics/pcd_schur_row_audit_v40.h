#pragma once

// V40: read-only, scalable row audit of the *frozen* coupled saddle system.
// Compare PCD auxiliary Fp/Kp to a PETSc-selfp-type surrogate
//   Sp = C - D diag(Auu)^(-1) G.
// This is NOT the exact Schur C-D Auu^(-1) G, a physical residual, a
// correction of the pressure BCs, or a nonlinear Jacobian. Nothing here
// modifies either the coupled matrix or the preconditioner.

#include "cfdx/core/linalg/sparse_matrix.h"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <stdexcept>
#include <vector>

namespace cfdx::physics {
struct PcdSchurRowFitV40 {
    std::size_t cell = 0;
    double sp_l2 = 0.0;
    double fp_l2 = 0.0;
    double kp_l2 = 0.0;
    double fp_best_multiplier = 0.0;
    double kp_best_multiplier = 0.0;
    double fp_scaled_shape_defect = 0.0; // [0,1], invariant to scalar units
    double kp_scaled_shape_defect = 0.0;
    std::size_t sp_nonzeros = 0;
    std::size_t fp_nonzeros = 0;
    bool gauge_like = false;
};

// Returns one row per pressure cell. O(nnz(D)*stencil(G)) work,
// O(n_pressure+one pressure row) auxiliary memory; no dense N*N matrix.
inline std::vector<PcdSchurRowFitV40> audit_pcd_schur_rows_v40(
    const cfdx::core::SparseMatrix& A,
    const cfdx::core::SparseMatrix& Fp,
    const cfdx::core::SparseMatrix& Kp,
    std::size_t n_pressure)
{
    if (n_pressure == 0 || n_pressure > std::numeric_limits<std::size_t>::max()/4 ||
        !A.is_consistent() || !Fp.is_consistent() || !Kp.is_consistent() ||
        A.n_rows() != 4*n_pressure || A.n_cols() != 4*n_pressure ||
        Fp.n_rows() != n_pressure || Fp.n_cols() != n_pressure ||
        Kp.n_rows() != n_pressure || Kp.n_cols() != n_pressure)
        throw std::invalid_argument("V40 PCD Schur audit: incomplete coupled/pressure CSR contract");
    const std::size_t nu = 3*n_pressure;
    std::vector<long double> invdiag(nu, 0.0L);
    for (std::size_t u=0; u<nu; ++u) {
        long double diagonal=0.0L;
        for (std::size_t k=A.row_offsets_data()[u]; k<A.row_offsets_data()[u+1]; ++k)
            if (A.columns_data()[k]==u) diagonal+=A.values_data()[k];
        if (!(diagonal>0.0L) || !std::isfinite(diagonal))
            throw std::invalid_argument("V40 PCD Schur audit: nonpositive/nonfinite velocity diagonal");
        invdiag[u]=1.0L/diagonal;
    }
    const auto score = [](const std::map<std::size_t,long double>& s,
                          const cfdx::core::SparseMatrix& p, std::size_t row,
                          double& norm, double& multiplier, std::size_t& nnz) {
        long double sn=0.0L,pn=0.0L,dot=0.0L;
        for (const auto& e:s) sn+=e.second*e.second;
        for (std::size_t k=p.row_offsets_data()[row];k<p.row_offsets_data()[row+1];++k) {
            const long double v=p.values_data()[k];
            if (!std::isfinite(v)) throw std::invalid_argument("V40 PCD Schur audit: nonfinite auxiliary coefficient");
            pn+=v*v; dot+=v*(s.count(p.columns_data()[k])?s.at(p.columns_data()[k]):0.0L);
            if (v != 0.0L) ++nnz;
        }
        if (!(sn>0.0L) || !(pn>0.0L) || !std::isfinite(sn) || !std::isfinite(pn))
            throw std::invalid_argument("V40 PCD Schur audit: zero/nonfinite pressure row");
        norm=static_cast<double>(std::sqrt(pn));
        multiplier=static_cast<double>(dot/pn);
        const long double ratio=(dot*dot)/(sn*pn);
        const double defect=static_cast<double>(std::sqrt(std::max(0.0L,1.0L-std::min(1.0L,ratio))));
        if (!std::isfinite(norm) || !std::isfinite(multiplier) || !std::isfinite(defect))
            throw std::runtime_error("V40 PCD Schur audit: nonfinite fit");
        return defect;
    };
    std::vector<PcdSchurRowFitV40> rows;
    rows.reserve(n_pressure);
    for (std::size_t i=0; i<n_pressure; ++i) {
        std::map<std::size_t,long double> s;
        std::size_t pressure_velocity_entries=0;
        // Row i of the original coupled matrix: [ D_i | C_i ].
        for (std::size_t d=A.row_offsets_data()[nu+i];d<A.row_offsets_data()[nu+i+1];++d) {
            const std::size_t col=A.columns_data()[d];
            const long double value=A.values_data()[d];
            if (!std::isfinite(value))
                throw std::invalid_argument("V40 PCD Schur audit: nonfinite coupled row");
            if (col>=nu) {s[col-nu]+=value; continue;}
            if (value == 0.0L) continue;
            ++pressure_velocity_entries;
            for (std::size_t k=A.row_offsets_data()[col];k<A.row_offsets_data()[col+1];++k)
                if (A.columns_data()[k]>=nu)
                    s[A.columns_data()[k]-nu]-=value*invdiag[col]*A.values_data()[k];
        }
        PcdSchurRowFitV40 r;
        r.cell=i;
        long double sn=0.0L;
        for (auto it=s.begin();it!=s.end();) {
            if (!std::isfinite(it->second))
                throw std::runtime_error("V40 PCD Schur audit: invalid surrogate Schur row");
            if (it->second==0.0L) {it=s.erase(it);continue;}
            sn+=it->second*it->second;
            ++it;
        }
        r.sp_nonzeros=s.size();
        r.sp_l2=static_cast<double>(std::sqrt(sn));
        r.gauge_like=pressure_velocity_entries==0 && s.size()==1 && s.count(i)==1;
        if (!(r.sp_l2>0.0) || !std::isfinite(r.sp_l2))
            throw std::invalid_argument("V40 PCD Schur audit: null/nonfinite surrogate row");
        r.fp_scaled_shape_defect=score(s,Fp,i,r.fp_l2,r.fp_best_multiplier,r.fp_nonzeros);
        std::size_t kp_nnz=0;
        r.kp_scaled_shape_defect=score(s,Kp,i,r.kp_l2,r.kp_best_multiplier,kp_nnz);
        rows.push_back(r);
    }
    return rows;
}
} // namespace cfdx::physics
