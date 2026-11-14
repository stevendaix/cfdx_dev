#include "cfdx/physics/sst_solver.h"
int main(){ cfdx::physics::TurbulenceTransportControls ctl; ctl.sst_anderson_m2_v20=true; return ctl.sst_anderson_m2_v20 ? 0:1; }
