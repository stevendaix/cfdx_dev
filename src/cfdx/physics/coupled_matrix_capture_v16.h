#pragma once
// Opt-in forensic export: exact monolithic matrix, right-hand side and guess
// at the point of the Krylov call. Never used to change the production solve.
// Maximum 6000 unknowns avoids accidental huge exports on physical meshes.

#include "cfdx/core/linalg/sparse_matrix.h"
#include "cfdx/core/linalg/vector.h"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>

namespace cfdx::physics {

inline void capture_coupled_matrix_v16(
    const cfdx::core::SparseMatrix& A,
    const cfdx::core::Vector& b,
    const cfdx::core::Vector& guess,
    const std::filesystem::path& prefix)
{
    if (prefix.empty() || prefix.filename().empty() ||
        !A.is_consistent() || A.n_rows()!=A.n_cols() ||
        A.n_rows()!=b.size() || A.n_cols()!=guess.size() ||
        A.n_rows()==0 || A.n_rows()>6000)
        throw std::invalid_argument("V16 coupled forensic capture: invalid dimensions/path or >6000 unknowns");
    const auto parent = prefix.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
    const auto write = [](const std::filesystem::path& target, const auto& output) {
        const auto tmp=target.string()+".partial";
        {
            std::ofstream f(tmp, std::ios::binary|std::ios::trunc);
            if (!f) throw std::runtime_error("V16 coupled capture: unable to create file");
            f << std::setprecision(17);
            output(f);
            f.flush();
            if (!f) throw std::runtime_error("V16 coupled capture: file I/O failed");
        }
        std::filesystem::rename(tmp,target);
    };
    write(prefix.string()+"_A.mtx", [&](std::ofstream& f) {
        f << "%%MatrixMarket matrix coordinate real general\n";
        f << A.n_rows() << ' ' << A.n_cols() << ' ' << A.nnz() << '\n';
        for (std::size_t row=0;row<A.n_rows();++row)
            for (std::size_t k=A.row_offsets_data()[row];k<A.row_offsets_data()[row+1];++k)
                f << row+1 << ' ' << A.columns_data()[k]+1 << ' '
                  << A.values_data()[k] << '\n';
    });
    const auto vector_file=[&](const char* suffix,const cfdx::core::Vector& value) {
        write(prefix.string()+suffix,[&](std::ofstream& f) {
            f << "%%MatrixMarket matrix array real general\n";
            f << value.size() << " 1\n";
            for(std::size_t i=0;i<value.size();++i) f << value(i) << '\n';
        });
    };
    vector_file("_b.mtx",b);
    vector_file("_x0.mtx",guess);
}

}  // namespace cfdx::physics
