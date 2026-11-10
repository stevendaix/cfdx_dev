#include "cfdx/physics/coupled_matrix_capture_v16.h"
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

#define REQUIRE(v) do { if (!(v)) throw std::runtime_error("V16 coupled capture failed: " #v); } while (false)

int main() {
    using namespace cfdx::core;
    const auto root=std::filesystem::temp_directory_path()/"cfdx_coupled_v16_capture_test";
    std::filesystem::create_directories(root);
    SparseMatrix A(2,2);
    A.push_back(0,0,2); A.push_back(0,1,-1); A.push_back(1,0,-3); A.push_back(1,1,4);
    A.finalize();
    Vector b(2),x0(2);b(0)=1;b(1)=2;
    const auto prefix=root/"case";
    cfdx::physics::capture_coupled_matrix_v16(A,b,x0,prefix);
    std::ifstream f(prefix.string()+"_A.mtx");
    std::string line;
    std::getline(f,line);REQUIRE(line=="%%MatrixMarket matrix coordinate real general");
    std::getline(f,line);REQUIRE(line=="2 2 4");
    std::getline(f,line);REQUIRE(line=="1 1 2");
    std::getline(f,line);REQUIRE(line=="1 2 -1");
    std::ifstream r(prefix.string()+"_b.mtx");
    std::getline(r,line); REQUIRE(line=="%%MatrixMarket matrix array real general");
    std::filesystem::remove_all(root);
}
