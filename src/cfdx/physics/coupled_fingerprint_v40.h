#pragma once
// Deterministic read-only FNV-1a/64 fingerprint of the ORIGINAL assembled
// coupled CSR and its b,x0. For cross-run *identity checks*, NOT a
// cryptographic SHA-256, nonlinear convergence certificate, or collision-proof
// proof. Fixed 64-bit little-endian word encoding is independent of host byte
// order; reject NaNs/Inf and unfinalized CSR before hashing.
#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include <cstdint>
#include <cstring>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>

namespace cfdx::physics {
struct CoupledFingerprintV40 {
    std::string A;
    std::string b;
    std::string x0;
};
inline std::uint64_t fnv64_word_v40(std::uint64_t hash, std::uint64_t value) {
    for (unsigned i=0;i<8;++i) {
        hash^=static_cast<unsigned char>((value>>(8*i))&0xffu);
        hash*=1099511628211ULL;
    }
    return hash;
}
inline std::uint64_t fnv64_double_v40(std::uint64_t hash,double value) {
    if (!std::isfinite(value)) throw std::invalid_argument("V40 fingerprint: nonfinite coefficient");
    std::uint64_t bits=0;
    static_assert(sizeof(bits)==sizeof(value),"V40 fingerprint requires binary64");
    std::memcpy(&bits,&value,sizeof(bits));
    return fnv64_word_v40(hash,bits);
}
inline std::string fnv64_hex_v40(std::uint64_t hash) {
    std::ostringstream out;
    out<<std::hex<<std::setfill('0')<<std::setw(16)<<hash;
    return out.str();
}
inline CoupledFingerprintV40 fingerprint_coupled_v40(
    const cfdx::core::SparseMatrix& A,
    const cfdx::core::Vector& b,
    const cfdx::core::Vector& x0) {
    if (!A.is_consistent() || !A.n_rows() || A.n_rows()!=A.n_cols() ||
        b.size()!=A.n_rows() || x0.size()!=A.n_cols())
        throw std::invalid_argument("V40 fingerprint: incompatible coupled matrix/vector contract");
    constexpr std::uint64_t seed=14695981039346656037ULL;
    std::uint64_t h=fnv64_word_v40(seed,A.n_rows());
    h=fnv64_word_v40(h,A.n_cols());
    h=fnv64_word_v40(h,A.nnz());
    for (std::size_t i=0;i<=A.n_rows();++i) h=fnv64_word_v40(h,A.row_offsets_data()[i]);
    for (std::size_t k=0;k<A.nnz();++k) {
        h=fnv64_word_v40(h,A.columns_data()[k]);
        h=fnv64_double_v40(h,A.values_data()[k]);
    }
    auto hv=[&](const cfdx::core::Vector& v) {
        auto out=fnv64_word_v40(seed,v.size());
        for (std::size_t i=0;i<v.size();++i) out=fnv64_double_v40(out,v(i));
        return fnv64_hex_v40(out);
    };
    return {fnv64_hex_v40(h),hv(b),hv(x0)};
}
} // namespace cfdx::physics
