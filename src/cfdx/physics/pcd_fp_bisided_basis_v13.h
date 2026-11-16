#pragma once

// V13: exact-operator evidence and reusable algebraic projector for a
// certified LEFT/RIGHT numerical nullspace pair of the ORIGINAL Fp.
// No matrix/gauge/boundary condition is changed. Supply an independently
// certified basis for the CURRENT operator, never reuse one across updates.
// Projection changes a PRECONDITIONER RHS only; the outer coupled physical
// true residual must always be computed against the unprojected equations.

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace cfdx::physics {

inline std::uint64_t pcd_fp_v13_hash_matrix(const cfdx::core::SparseMatrix& A) {
    std::uint64_t hash=14695981039346656037ULL;
    const auto mix=[&](const unsigned char* bytes,std::size_t count) {
        for (std::size_t i=0;i<count;++i) {
            hash^=static_cast<std::uint64_t>(bytes[i]);
            hash*=1099511628211ULL;
        }
    };
    const std::uint64_t n=A.n_rows(), nnz=A.nnz();
    mix(reinterpret_cast<const unsigned char*>(&n),sizeof(n));
    mix(reinterpret_cast<const unsigned char*>(&nnz),sizeof(nnz));
    mix(reinterpret_cast<const unsigned char*>(A.row_offsets_data()),(A.n_rows()+1)*sizeof(cfdx::core::SparseMatrix::Index));
    mix(reinterpret_cast<const unsigned char*>(A.columns_data()),A.nnz()*sizeof(cfdx::core::SparseMatrix::Index));
    mix(reinterpret_cast<const unsigned char*>(A.values_data()),A.nnz()*sizeof(double));
    return hash;
}

struct PcdFpBiNullspaceV13 {
    std::vector<cfdx::core::Vector> left;
    std::vector<cfdx::core::Vector> right;
    std::uint64_t operator_hash = 0;
    std::size_t size = 0;
    double max_left_defect = 0.0;
    double max_right_defect = 0.0;

    static double dot(const cfdx::core::Vector& a, const cfdx::core::Vector& b) {
        if (a.size()!=b.size() || !a.is_valid() || !b.is_valid())
            throw std::invalid_argument("Fp V13 invalid dot arguments");
        long double sum=0.;
        for (std::size_t i=0;i<a.size();++i) sum += static_cast<long double>(a(i))*b(i);
        return static_cast<double>(sum);
    }

    static PcdFpBiNullspaceV13 certify(const cfdx::core::SparseMatrix& A,
        std::vector<cfdx::core::Vector> left_basis,
        std::vector<cfdx::core::Vector> right_basis,
        double normalized_tolerance = 1e-12)
    {
        const auto n=A.n_rows();
        if (n==0 || n!=A.n_cols() || !A.is_consistent() ||
            left_basis.size()!=right_basis.size() ||
            left_basis.empty() || left_basis.size()>8 ||
            !(normalized_tolerance>0. && normalized_tolerance<=1e-10))
            throw std::invalid_argument("Fp V13 invalid two-sided nullspace dimensions");
        const auto* rows=A.row_offsets_data();
        const auto* cols=A.columns_data();
        const auto* val=A.values_data();
        long double frobenius2=0.;
        for (std::size_t j=0;j<A.nnz();++j) {
            if (!std::isfinite(val[j])) throw std::invalid_argument("Fp V13 nonfinite matrix");
            frobenius2+=static_cast<long double>(val[j])*val[j];
        }
        const auto normA=std::sqrt(frobenius2);
        if (!(normA>0.)) throw std::invalid_argument("Fp V13 zero matrix");
        auto result=PcdFpBiNullspaceV13{};
        result.size=n;
        result.operator_hash=pcd_fp_v13_hash_matrix(A);
        result.left=std::move(left_basis);
        result.right=std::move(right_basis);
        const auto check=[&](const std::vector<cfdx::core::Vector>& vectors, bool transpose, double& maximum) {
            for (std::size_t j=0;j<vectors.size();++j) {
                const auto& x=vectors[j];
                if (x.size()!=n || !x.is_valid())
                    throw std::invalid_argument("Fp V13 nonfinite/wrong-sized nullspace vector");
                if (std::abs(dot(x,x)-1.)>1e-8)
                    throw std::invalid_argument("Fp V13 nullspace vector not normalized");
                for (std::size_t k=0;k<j;++k)
                    if (std::abs(dot(x,vectors[k]))>1e-8)
                        throw std::invalid_argument("Fp V13 basis is not orthogonal");
                std::vector<long double> residual(n,0.L);
                for (std::size_t row=0;row<n;++row)
                    for (std::size_t index=rows[row];index<rows[row+1];++index)
                        if (transpose)
                            residual[cols[index]]+=static_cast<long double>(val[index])*x(row);
                        else residual[row]+=static_cast<long double>(val[index])*x(cols[index]);
                long double squared=0.;
                for (long double term:residual) squared+=term*term;
                const auto defect=static_cast<double>(std::sqrt(squared)/normA);
                if (!std::isfinite(defect) || defect>normalized_tolerance)
                    throw std::invalid_argument("Fp V13 candidate is NOT a certified left/right null vector");
                maximum=std::max(maximum,defect);
            }
        };
        check(result.left,true,result.max_left_defect);
        check(result.right,false,result.max_right_defect);
        return result;
    }

    // Low-rank, matrix-free deflation: A_d = A + W Z^T. This supplies a
    // unit singular value for each certified missing singular direction.
    // W and Z are respectively LEFT and RIGHT null bases of this exact A;
    // they need not span the same space for nonsymmetric transport.
    // No sparse fill-in and no modification of the physical coupled matrix.
    cfdx::core::Vector deflated_matvec(
        const cfdx::core::SparseMatrix& A,
        const cfdx::core::Vector& x) const {
        if (pcd_fp_v13_hash_matrix(A) != operator_hash ||
            x.size() != size || !x.is_valid())
            throw std::invalid_argument("V14 PCD deflation: stale matrix or vector");
        cfdx::core::Vector result(size,0.0);
        A.matvec_into(x,result);
        for (std::size_t j=0;j<left.size();++j) {
            const double factor=dot(right[j],x);
            for (std::size_t i=0;i<size;++i)
                result(i)+=left[j](i)*factor;
        }
        if (!result.is_valid())
            throw std::runtime_error("V14 PCD deflation: nonfinite operator image");
        return result;
    }

    // Independent residual of the ORIGINAL operator against the compatible
    // part of b. Not the residual reported by the deflated Krylov solver.
    double projected_residual_relative(
        const cfdx::core::SparseMatrix& A,
        const cfdx::core::Vector& original_rhs,
        const cfdx::core::Vector& solution) const {
        if (pcd_fp_v13_hash_matrix(A) != operator_hash ||
            original_rhs.size() != size || solution.size() != size ||
            !original_rhs.is_valid() || !solution.is_valid())
            throw std::invalid_argument("V14 PCD residual: stale/invalid operands");
        const auto compatible=compatible_rhs(original_rhs);
        cfdx::core::Vector actual(size,0.0);
        A.matvec_into(solution,actual);
        long double defect=0.0L, norm=0.0L;
        for (std::size_t i=0;i<size;++i) {
            defect=std::hypot(defect,static_cast<long double>(compatible(i))-actual(i));
            norm=std::hypot(norm,static_cast<long double>(compatible(i)));
        }
        return static_cast<double>(defect/std::max(norm,1e-300L));
    }

    // Orthogonal projection of a pressure-space RHS onto range(A).
    // The discarded component is the irreducible ORIGINAL Fp residual floor.
    cfdx::core::Vector compatible_rhs(const cfdx::core::Vector& original) const {
        if (original.size()!=size || !original.is_valid())
            throw std::invalid_argument("Fp V13 invalid original RHS");
        auto output=original;
        for (const auto& v:left) {
            const double factor=dot(v,original);
            for (std::size_t i=0;i<size;++i) output(i)-=factor*v(i);
        }
        if (!output.is_valid()) throw std::runtime_error("Fp V13 nonfinite projected RHS");
        return output;
    }

    // Canonical minimal-norm gauge: removes RIGHT null components, not LEFT.
    cfdx::core::Vector canonical_solution(const cfdx::core::Vector& candidate) const {
        if (candidate.size()!=size || !candidate.is_valid())
            throw std::invalid_argument("Fp V13 invalid solution vector");
        auto output=candidate;
        for (const auto& v:right) {
            const double factor=dot(v,candidate);
            for (std::size_t i=0;i<size;++i) output(i)-=factor*v(i);
        }
        if (!output.is_valid()) throw std::runtime_error("Fp V13 nonfinite projected gauge");
        return output;
    }
    double incompatibility_relative(const cfdx::core::Vector& b) const {
        const auto projected=compatible_rhs(b);
        long double removed2=0., rhs2=0.;
        for (std::size_t i=0;i<size;++i) {
            removed2+=static_cast<long double>(b(i)-projected(i))*(b(i)-projected(i));
            rhs2+=static_cast<long double>(b(i))*b(i);
        }
        return rhs2>0. ? static_cast<double>(std::sqrt(removed2/rhs2)) : 0.;
    }
};

// CFPB13 basis produced by scripts/certify_pcd_fp_bi_nullspace_v13.py. The
// independent SVD proof must also be retained; this reader verifies each
// basis numerically against the EXACT current matrix after fingerprint check.
inline PcdFpBiNullspaceV13 load_pcd_fp_bisided_basis_v13(
    const std::string& path, const cfdx::core::SparseMatrix& A) {
    std::ifstream in(path,std::ios::binary);
    if (!in) throw std::runtime_error("Fp V13 basis file unavailable");
    char magic[8];in.read(magic,8);
    if (!in || std::memcmp(magic,"CFDXNB13",8))
        throw std::invalid_argument("Fp V13 incorrect basis file header");
    std::uint64_t header[3];
    in.read(reinterpret_cast<char*>(header),sizeof(header));
    const auto n=static_cast<std::size_t>(header[0]);
    const auto k=static_cast<std::size_t>(header[1]);
    if (!in || n==0 || n!=A.n_rows() || k==0 || k>8 ||
        header[2]!=pcd_fp_v13_hash_matrix(A))
        throw std::invalid_argument("Fp V13 basis does not match current original matrix");
    std::vector<cfdx::core::Vector> left,right;
    const auto read_basis=[&](std::vector<cfdx::core::Vector>& out) {
        for (std::size_t j=0;j<k;++j) {
            cfdx::core::Vector x(n);
            in.read(reinterpret_cast<char*>(x.data()),static_cast<std::streamsize>(n*sizeof(double)));
            if (!in) throw std::invalid_argument("Fp V13 truncated basis");
            out.push_back(std::move(x));
        }
    };
    read_basis(left);read_basis(right);
    if (in.peek()!=std::char_traits<char>::eof())
        throw std::invalid_argument("Fp V13 trailing basis payload");
    return PcdFpBiNullspaceV13::certify(A,std::move(left),std::move(right));
}

} // namespace cfdx::physics
