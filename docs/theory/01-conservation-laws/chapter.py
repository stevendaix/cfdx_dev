# %% [markdown]
# Conservation Laws

## 1.1 Generic balance

\[
\frac{\partial q}{\partial t}+\nabla\cdot\mathbf F=s.
\]
Integration gives
\[
\frac{d}{dt}\int_{V_P}q\,dV+
\oint_{\partial V_P}\mathbf F\cdot\mathbf n\,dA
=\int_{V_P}s\,dV.
\]

## 1.2 Mass

\[
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf u)=0.
\]
The cell equation is
\[
\frac{d}{dt}(\rho_PV_P)+\sum_f\dot m_f=0,
\qquad
\dot m_f=\rho_f\mathbf u_f\cdot\mathbf S_f.
\]

## 1.3 Momentum

\[
\frac{\partial(\rho\mathbf u)}{\partial t}
+\nabla\cdot(\rho\mathbf u\otimes\mathbf u)
=-\nabla p+\nabla\cdot\boldsymbol\tau+\rho\mathbf f.
\]
The integral equation contains convective momentum, pressure traction, viscous traction and body force.

## 1.4 Energy

With total specific energy \(E=e+|\mathbf u|^2/2\),
\[
\frac{\partial(\rho E)}{\partial t}
+\nabla\cdot[(\rho E+p)\mathbf u]
=\nabla\cdot(\boldsymbol\tau\mathbf u-\mathbf q)
+\rho\mathbf f\cdot\mathbf u+\dot q_v.
\]

## 1.5 Species

\[
\frac{\partial(\rho Y_i)}{\partial t}
+\nabla\cdot(\rho\mathbf uY_i)
=-\nabla\cdot\mathbf J_i+\dot\omega_i.
\]
For a conservative mixture,
\[
\sum_iY_i=1,\qquad \sum_i\dot\omega_i=0.
\]

## 1.6 Internal-face cancellation

For two adjacent cells sharing one physical face,
\[
\mathbf S_{f,N}=-\mathbf S_{f,P}.
\]
A single conservative interface flux therefore satisfies
\[
F_{N,f}=-F_{P,f}.
\]
Summing all cell equations cancels internal faces exactly, leaving only physical-domain boundary fluxes and sources.

## 1.7 Global conservation diagnostic

\[
R_C=
\frac{d}{dt}\sum_Pq_PV_P+
\sum_{f\in\partial\Omega}\Phi_f-
\sum_Ps_PV_P.
\]
A normalized defect can be
\[
\epsilon_C=\frac{|R_C|}{Q_{\mathrm{ref}}},
\]
where the reference scale must be explicitly documented.

## 1.8 CFDX traceability

src/cfdx/core/numerics/conservation.h  
src/cfdx/core/numerics/flux.h  
src/cfdx/core/numerics/integrate.h  
src/cfdx/physics/finite_volume_transport.h

Verification must test face antisymmetry and global balances independently of solver residuals.

# %%
from __future__ import annotations
import numpy as np
S=np.array([1.2,-0.4,0.7])
assert np.allclose(S+(-S),0.0)
