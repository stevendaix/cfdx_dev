#include "cfdx/physics/pcd_pressure_operator_capture_v38.h"

#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

using namespace cfdx::core;
using namespace cfdx::physics;

static void require(bool condition, const char* reason) {
    if (!condition) throw std::runtime_error(reason);
}

static SparseMatrix diagonal(double multiplier) {
    SparseMatrix result(3,3);
    for (std::size_t i=0;i<3;++i)
        result.push_back(i,i,multiplier*(static_cast<double>(i)+1.0));
    result.finalize();
    return result;
}

static void check(const std::filesystem::path& path, double expected) {
    const auto capture=read_pcd_fp_capture(path);
    require(capture.matrix.n_rows()==3 && capture.matrix.nnz()==3,
            "pressure operator capture dimension changed");
    for (std::size_t i=0;i<3;++i) {
        require(std::abs(capture.matrix(i,i)-expected*static_cast<double>(i+1))<1e-15,
                "pressure operator value changed");
        require(capture.rhs(i)==0.0 && capture.initial(i)==0.0,
                "V38 matrix-only placeholder vectors not zero");
    }
}

int main() {
    const auto temp=std::filesystem::temp_directory_path() /
        ("cfdx_pcd_capture_v38_" + std::to_string(
            std::chrono::steady_clock::now().time_since_epoch().count()));
    try {
        std::filesystem::create_directories(temp);
        const auto prefix=(temp/"triplet").string();
        const auto mp=diagonal(2.0), kp=diagonal(3.0), fp=diagonal(5.0);
        try {
            (void)capture_pcd_pressure_operators_v38("",mp,kp,fp);
            throw std::runtime_error("empty prefix was accepted");
        } catch (const std::invalid_argument&) {}
        for (std::size_t slot=1;slot<=3;++slot) {
            require(capture_pcd_pressure_operators_v38(prefix,mp,kp,fp)==slot,
                    "capture slot numbering is incorrect");
            const auto base=prefix+"_state"+std::to_string(slot);
            check(base+"_mp.cfp",2.0);
            check(base+"_kp.cfp",3.0);
            check(base+"_fp.cfp",5.0);
        }
        require(capture_pcd_pressure_operators_v38(prefix,mp,kp,fp)==0,
                "unbounded PCD operator capture");
        require(!std::filesystem::exists(prefix+"_state4_mp.cfp"),
                "a fourth pressure operator sample was written");
        std::filesystem::remove_all(temp);
        std::cout << "PCD_PRESSURE_OPERATORS_V38_PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::filesystem::remove_all(temp);
        std::cerr << "PCD_PRESSURE_OPERATORS_V38_FAIL " << e.what() << '\n';
        return 2;
    }
}
