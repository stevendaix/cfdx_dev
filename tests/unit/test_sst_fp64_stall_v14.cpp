#include "cfdx/physics/sst_fp64_stall_v14.h"
#include <iostream>
#include <stdexcept>
#include <limits>
using namespace cfdx::physics;
int main() {
 try {
  SSTFp64StallV14 s;
  if (s.observe(1e-12,1e-13,false) || s.observe(0,0,false) ||
      s.observe(0,0,false) || !s.observe(0,0,false))
    throw std::runtime_error("precision plateau not detected correctly");
  if (s.observe(0,0,true) || s.identical_unconverged!=0)
    throw std::runtime_error("converged step incorrectly marked stalled");
  if (s.observe(0,0,false) || s.observe(1e-300,0,false) || s.identical_unconverged!=0)
    throw std::runtime_error("progress failed to reset counter");
  bool reject=false;
  try {s.observe(std::numeric_limits<double>::quiet_NaN(),0,false);}
  catch(const std::invalid_argument&){reject=true;}
  if(!reject) throw std::runtime_error("NaN escaped stall contract");
  std::cout<<"SST_FP64_STALL_V14_PASS same_tolerance=1\n";
  return 0;
 } catch(const std::exception& e) {
  std::cerr<<"SST_FP64_STALL_V14_FAIL "<<e.what()<<"\n";return 1;
 }
}
