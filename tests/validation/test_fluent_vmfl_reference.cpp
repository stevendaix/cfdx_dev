#include "common/test_harness.h"

#include <cmath>
#include <iostream>

int main()
{
    try {
        // VMFL004 — Plain Couette flow with pressure gradient.
        // Fluent 2026 R1 setup:
        // L=1.5 m, W=1 m, rho=1 kg/m3, mu=1 Pa.s,
        // U_wall=3 m/s, dp/dx=-12 Pa/m.
        // Exact solution for H=1 m:
        //   u(y) = U_wall*y/H + (dp/dx)/(2*mu)*y*(y-H)
        //        = 9*y - 6*y^2.
        {
            constexpr double H = 1.0;
            constexpr double U = 3.0;
            constexpr double mu = 1.0;
            constexpr double dpdx = -12.0;

            for (const double y : {0.0, 0.25, 0.5, 0.75, 1.0}) {
                const double u =
                    U*y/H + dpdx/(2.0*mu)*y*(y-H);
                const double expected = 9.0*y - 6.0*y*y;
                if (std::abs(u - expected) > 1e-14)
                    throw std::runtime_error("VMFL004 reference mismatch");
            }

            // Integrated wall-normal profile average:
            // integral_0^H u dy / H = U/2 - dpdx*H^2/(12 mu) = 1.5+1 = 2.5 m/s.
            const double ubar = U/2.0 - dpdx*H*H/(12.0*mu);
            if (std::abs(ubar - 2.5) > 1e-14)
                throw std::runtime_error("VMFL004 bulk velocity mismatch");
        }

        // VMFL005 — Poiseuille flow in a pipe.
        // Fluent 2026 R1 setup:
        // Re_D=500, Ubar=2 m/s, L=0.1 m, R=0.00125 m,
        // rho=1 kg/m3, mu=1e-5 Pa.s.
        // Hagen-Poiseuille pressure drop:
        //   Delta-p = 8 mu L Ubar / R^2 = 10.24 Pa.
        {
            constexpr double Ubar = 2.0;
            constexpr double L = 0.1;
            constexpr double R = 0.00125;
            constexpr double mu = 1e-5;
            const double dp = 8.0*mu*L*Ubar/(R*R);
            if (std::abs(dp - 10.24) > 1e-12)
                throw std::runtime_error("VMFL005 pressure-drop oracle mismatch");

            const double Re = 2.0*R*Ubar/(mu/1.0);
            if (std::abs(Re - 500.0) > 1e-12)
                throw std::runtime_error("VMFL005 Reynolds-number mismatch");
        }

        // VMFL003 — turbulent pipe pressure drop.
        // Fluent reports the Moody-chart target as 21744 Pa for:
        // L=2 m, R=0.002 m, rho=1.225 kg/m3,
        // mu=1.7894e-5 Pa.s, U=50 m/s, Re=1.37e4.
        // Keep the published target as a reference datum; CFDX must
        // eventually obtain it from the coupled turbulent pipe solver.
        {
            constexpr double target_dp = 21744.0;
            constexpr double L = 2.0;
            constexpr double D = 0.004;
            constexpr double rho = 1.225;
            constexpr double mu = 1.7894e-5;
            constexpr double U = 50.0;
            const double Re = rho*U*D/mu;
            if (std::abs(Re - 13700.0) > 80.0)
                throw std::runtime_error("VMFL003 Reynolds setup mismatch");
            if (!(target_dp > 0.0))
                throw std::runtime_error("VMFL003 invalid reference");
            (void)L;
        }

        std::cout << "FLUENT_VMFL_REFERENCE: PASS\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "FLUENT_VMFL_REFERENCE: FAIL: " << e.what() << "\n";
        return 1;
    }
}
