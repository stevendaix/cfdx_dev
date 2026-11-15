#include "cfdx/physics/pressure_native_defect_v80.h"
#include <cmath>
#include <iostream>
#include <limits>
#include <stdexcept>

using namespace cfdx::core;
using namespace cfdx::physics;
static void need(bool ok, const char* message) {
    if (!ok) throw std::runtime_error(message);
}
int main() {
    // A nearly singular 2x2 face-Laplacian: the single-double pressure
    // correction loses digits, but high/low retains the Newton increment.
    SparseMatrix a(2,2);
    const double tiny=1.e-8;
    a.push_back(0,0,1.0); a.push_back(0,1,-1.0);
    a.push_back(1,0,-1.0);a.push_back(1,1,1.0+tiny);
    a.finalize();
    Vector b(2,0.0);b(0)=1.0;
    Vector high(2,0.0),low(2,0.0);
    const auto inner=[&](const Vector& rhs,Vector& d) {
        const long double eps=static_cast<long double>(a(1,1))-1.0L;
        d(1)=static_cast<double>((static_cast<long double>(rhs(0))+rhs(1))/eps);
        d(0)=static_cast<double>(static_cast<long double>(rhs(0))+d(1));
        return true;
    };
    const auto report=refine_native_pressure_algebraic_v80(
        a,b,high,low,1.e-11,8,inner);
    need(report.converged && report.final_relative < 1.e-11,
         "V80 compensated native defect must close the true CSR residual");
    need(report.accepted>0 && !report.production_commit_authorized,
         "V80 CSR-only evidence must never grant a physical state commit");
    need(std::isfinite(high(0)) && std::isfinite(low(0)),
         "V80 correction must remain finite");
    Vector failed_hi(2,0.0), failed_lo(2,0.0);
    const auto fail=refine_native_pressure_algebraic_v80(
        a,b,failed_hi,failed_lo,1.e-11,8,
        [](const Vector&, Vector& delta){
            delta(0)=std::numeric_limits<double>::quiet_NaN();
            delta(1)=1.; return true;
        });
    need(!fail.converged && fail.accepted==0 &&
         failed_hi(0)==0.0 && failed_lo(1)==0.0,
         "V80 nonfinite directions must fail closed without mutation");
    bool rejected=false;
    try { refine_native_pressure_algebraic_v80(a,b,high,low,0.,8,inner); }
    catch(const std::invalid_argument&) {rejected=true;}
    need(rejected,"V80 must refuse invalid target tolerance");
    std::cout << "NATIVE_PRESSURE_DEFECT_V80_PASS relative="
              << report.final_relative << " accepted=" << report.accepted
              << " committed=0\n";
}
