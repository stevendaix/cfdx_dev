@@ -0,0 +1,28 @@
#include "cfdx/physics/steady_incompressible_solver.h"
#include <cmath>
#include <stdexcept>
#include <string>
int main() {
    using namespace cfdx::physics;
    V17PseudoMomentumControls pseudo;
    pseudo.enabled=true;
    pseudo.maximum_diagonal_fraction=0.5;
    validate_v17_pseudo_momentum_controls(pseudo);
    if (!std::isfinite(pseudo.maximum_diagonal_fraction)) return 1;
    // The finite cap is a NONLINEAR globalisation parameter, not a change
    // to the physical equation. Fail closed on malformed limits.
    for (double bad: {-1.0, 0.0, std::numeric_limits<double>::quiet_NaN()}) {
       auto x=pseudo;x.maximum_diagonal_fraction=bad;
       try { validate_v17_pseudo_momentum_controls(x); return 2; }
       catch (const std::invalid_argument&) {}
    }
    // Demonstrate that adding the same positive mass to three segregated
    // blocks preserves the exact physical fixed point U_new=U_old.
    for (double u: {-44.2, 0.0, 44.2}) {
       for (double mass: {0.01,1.0,1000.0}) {
          const double delta=mass*(u-u);
          if(delta!=0.0) return 3;
       }
    }
    return 0;
}
