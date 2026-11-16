#include "cfdx/physics/pressure_split_oracle_v79.h"
#include <cstdio>
#include <stdexcept>
#include <cmath>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string>
static void need(bool condition, const char* msg) {
    if(!condition) throw std::runtime_error(msg);
}
int main() {
    using namespace cfdx::core;
    using namespace cfdx::physics;
    SparseMatrix a(2,2);
    a.push_back(0,0,2.0);a.push_back(0,1,-1.0);
    a.push_back(1,0,-1.0);a.push_back(1,1,2.0);a.finalize();
    Vector b(2,0.0);b(0)=1.0;b(1)=0.0;
    PressureSplitOracleV79 trial(2);
    trial.high(0)=2.0/3.;trial.high(1)=1.0/3.;
    const auto bad=audit_pressure_split_oracle_v79(a,b,trial,1e-17);
    need(!bad.csr_closed,"false pressure proof accepted");
    need(!bad.production_commit_authorized,"production promotion forbidden");
    trial.low(0)=static_cast<double>((static_cast<long double>(2)/3)-trial.high(0));
    trial.low(1)=static_cast<double>((static_cast<long double>(1)/3)-trial.high(1));
    const auto good=audit_pressure_split_oracle_v79(a,b,trial,1e-15);
    need(good.csr_closed && !good.production_commit_authorized,"good algebraic test failed");
    const std::string path="v79_test_oracle.split";
    {
        std::ofstream out(path);
        out << "CFDX_PRESSURE_SPLIT_V79\n2 " << a.nnz() << "\n";
        out.precision(17);
        for(std::size_t i=0;i<2;++i)out << trial.high(i) << ' ' << trial.low(i) << '\n';
    }
    auto loaded=read_pressure_split_oracle_v79(path,2,a.nnz());
    need(std::abs(loaded.low(0)-trial.low(0))<1e-30,"split low lost");
    bool rejected=false;
    try { (void)read_pressure_split_oracle_v79(path,3,a.nnz()); }
    catch(const std::runtime_error&) { rejected=true; }
    need(rejected,"bad or malformed split accepted");
    rejected=false;
    try { (void)audit_pressure_split_oracle_v79(a,b,loaded,0.0); }
    catch(const std::invalid_argument&) { rejected=true; }
    need(rejected,"bad or malformed split accepted");
    std::remove(path.c_str());
    return 0;
}