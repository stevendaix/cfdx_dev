#pragma once

#include "cfdx/core/linalg/hypre_amg.h"
#include "cfdx/core/linalg/preconditioner.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace cfdx::core {

// Experimental line-smoothed AMG for strongly anisotropic scalar matrices.
// Strong mutual off-diagonal connections form disjoint, non-cyclic lines.
// Each line is solved by a guarded tridiagonal factorization; the existing
// native AMG supplies global coarse-grid corrections. The compound operator
// is generally nonsymmetric: ONLY use with GMRES/FGMRES, never with CG.
// No claim of physical convergence is made by successful preconditioner setup.
class AnisotropicLineAMGPreconditioner final : public Preconditioner {
public:
    bool setup(const SparseMatrix& a) override {
        // Build a complete replacement first. A failed hierarchy rebuild
        // must not destroy the preconditioner used by an accepted step.
        AnisotropicLineAMGPreconditioner next;
        if (!next.setup_candidate(a)) return false;
        std::swap(n_, next.n_);
        std::swap(matrix_, next.matrix_);
        lines_.swap(next.lines_);
        amg_.swap(next.amg_);
        return true;
    }

private:
    bool setup_candidate(const SparseMatrix& a) {
        if (a.n_rows() == 0 || a.n_rows() != a.n_cols() || !a.is_consistent())
            return false;
        const std::size_t n = a.n_rows();
        const auto* ro=a.row_offsets_data();
        const auto* ci=a.columns_data();
        const auto* av=a.values_data();
        std::vector<double> diagonal(n,0.0);
        for (std::size_t i=0;i<n;++i) {
            for (std::size_t k=ro[i];k<ro[i+1];++k) {
                if (!std::isfinite(av[k]) || ci[k]>=n) return false;
                if (ci[k]==i) diagonal[i]+=av[k];
            }
            if (!(diagonal[i]>0.0) || !std::isfinite(diagonal[i])) return false;
        }
        struct Edge { std::size_t i,j; double score; };
        std::vector<Edge> edges;
        for (std::size_t i=0;i<n;++i) {
            for (std::size_t k=ro[i];k<ro[i+1];++k) {
                const std::size_t j=ci[k];
                if (j<=i || av[k]>=0.0) continue;
                const double rev=a(j,i);
                if (!(rev<0.0)) continue;
                const double s=std::min(std::abs(av[k])/diagonal[i],
                                        std::abs(rev)/diagonal[j]);
                if (s>=0.12) edges.push_back({i,j,s});
            }
        }
        std::sort(edges.begin(),edges.end(),[](const Edge& x,const Edge& y){
            if (x.score!=y.score) return x.score>y.score;
            return x.i<y.i || (x.i==y.i && x.j<y.j);
        });
        std::vector<std::size_t> parent(n),size(n,1);
        std::iota(parent.begin(),parent.end(),0);
        const auto find=[&](std::size_t x) {
            while (parent[x]!=x) x=parent[x]=parent[parent[x]];
            return x;
        };
        std::vector<std::vector<std::size_t>> nbr(n);
        for (const auto& e:edges) {
            const auto ri=find(e.i),rj=find(e.j);
            if (ri==rj || nbr[e.i].size()>=2 || nbr[e.j].size()>=2 ||
                size[ri]+size[rj]>1024) continue;
            nbr[e.i].push_back(e.j);
            nbr[e.j].push_back(e.i);
            parent[rj]=ri; size[ri]+=size[rj];
        }
        std::vector<bool> visited(n,false);
        for (std::size_t i=0;i<n;++i) {
            if (visited[i] || nbr[i].size()>1) continue;
            Line line;
            std::size_t current=i,previous=n;
            while (current<n && !visited[current]) {
                line.nodes.push_back(current);
                visited[current]=true;
                std::size_t next=n;
                for (auto j:nbr[current]) if (j!=previous) next=j;
                previous=current; current=next;
            }
            if (line.nodes.empty()) return false;
            const auto m=line.nodes.size();
            line.lower.resize(m,0.0);
            line.upper.resize(m,0.0);
            line.pivot.resize(m);
            for (std::size_t k=0;k<m;++k) {
                const auto cell=line.nodes[k];
                if (k) line.lower[k]=a(cell,line.nodes[k-1]);
                if (k+1<m) line.upper[k]=a(cell,line.nodes[k+1]);
                const double pivot=diagonal[cell]-(k?line.lower[k]*line.upper[k-1]/line.pivot[k-1]:0.0);
                const double reference=diagonal[cell]+std::abs(line.lower[k])+std::abs(line.upper[k]);
                if (!(pivot>128.0*std::numeric_limits<double>::epsilon()*reference) ||
                    !std::isfinite(pivot)) { clear(); return false; }
                line.pivot[k]=pivot;
            }
            lines_.push_back(std::move(line));
        }
        for (const auto hit:visited) if (!hit) { clear(); return false; }
        matrix_=a;
        amg_=std::make_unique<NativeBoomerAMGPreconditioner>(false,false);
        if (!amg_->setup(a)) { clear(); return false; }
        n_=n;
        return true;
    }

public:
    bool apply(const Vector& rhs,Vector& out) const override {
        if (n_==0 || rhs.size()!=n_ || out.size()!=n_ || !amg_ || !rhs.is_valid())
            return false;
        // The full preconditioner uses the right hand side three times.
        // Preserve it before the first line sweep when input/output alias.
        // The ordinary (non-aliased) Krylov path never makes this copy.
        Vector saved_rhs;
        const Vector* source = &rhs;
        if (&rhs == &out) {
            saved_rhs = rhs;
            source = &saved_rhs;
        }
        if (!line_solve(*source,out)) return false;
        // Reuse the AMG correction as the post-line scratch buffer.
        Vector residual(n_,0.0), correction(n_,0.0);
        // Reuse one pre-sized CSR matvec output in both residual evaluations.
        // Avoid the Vector->std::vector copy and two heap-allocated matvecs.
        matrix_.matvec_into(out,residual);
        for (std::size_t i=0;i<n_;++i) residual(i)=(*source)(i)-residual(i);
        if (!amg_->apply(residual,correction)) return false;
        for (std::size_t i=0;i<n_;++i) out(i)+=correction(i);
        matrix_.matvec_into(out,residual);
        for (std::size_t i=0;i<n_;++i) residual(i)=(*source)(i)-residual(i);
        if (!line_solve(residual,correction)) return false;
        for (std::size_t i=0;i<n_;++i) {
            out(i)+=correction(i);
            if (!std::isfinite(out(i))) return false;
        }
        return true;
    }
    const char* name() const override { return "AnisotropicLine-AMG"; }
    std::size_t number_of_lines() const noexcept { return lines_.size(); }
    std::size_t longest_line() const noexcept {
        std::size_t m=0;
        for (const auto& l:lines_) m=std::max(m,l.nodes.size());
        return m;
    }
private:
    struct Line {
        std::vector<std::size_t> nodes;
        std::vector<double> lower,upper,pivot;
    };
    bool line_solve(const Vector& r,Vector& z) const {
        for (const auto& line:lines_) {
            const auto m=line.nodes.size();
            // Store Thomas forward substitution directly in z. Each line is
            // disjoint; no O(line_length) allocation in every Krylov apply.
            for (std::size_t k=0;k<m;++k) {
                const std::size_t i=line.nodes[k];
                const double numerator=r(i)-
                    (k?line.lower[k]*z(line.nodes[k-1]):0.0);
                z(i)=numerator/line.pivot[k];
                if (!std::isfinite(z(i))) return false;
            }
            for (std::size_t k=m;k-->0;) {
                const std::size_t i=line.nodes[k];
                if (k+1<m)
                    z(i)-=line.upper[k]*z(line.nodes[k+1])/line.pivot[k];
                if (!std::isfinite(z(i))) return false;
            }
        }
        return true;
    }
    void clear() { n_=0; lines_.clear(); amg_.reset(); matrix_=SparseMatrix{}; }
    std::size_t n_=0;
    SparseMatrix matrix_;
    std::vector<Line> lines_;
    std::unique_ptr<NativeBoomerAMGPreconditioner> amg_;
};
} // namespace cfdx::core
