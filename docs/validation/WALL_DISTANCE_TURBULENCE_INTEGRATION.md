# Wall-distance integration and turbulence coupling

This PR closes the integration contract between the wall-distance service from #470 and the wall-sensitive turbulence kernels qualified by #474.

## Mathematical contract

For a positive cell-centred distance d:
- y+ = d u_tau / nu
- DES: ell_DES = min(d, C_DES Delta), with transition ratio d/(C_DES Delta)
- DDES: ell_DDES = d - f_d max(0, d-C_DES Delta), f_d = 1-tanh[(8 r_d)^3]
- IDDES uses the same distance-controlled hybrid scale with shielding/stress blending.
- SST uses d in F1/F2 and therefore in the blended closure.
- SA uses d in modified vorticity and wall-destruction terms.

## Verification matrix

| Path | Evidence |
|---|---|
| Exact/BVH distance | Existing #470 point-to-triangle and complex-geometry oracle |
| Distance field | finite/positive fluid-field contract |
| y+ | independent recomputation |
| SST | F1/F2 sensitivity to d |
| SA | modified-vorticity sensitivity to d |
| DES | d/(C_DES Delta) and length-scale contract |
| DDES | distance-driven shielding and length-scale contract |
| IDDES | distance-driven hybrid length-scale contract |

## Scope boundary

This is integration and closure verification, not physical validation. It does not promote any model from KERNEL_ONLY to SOLVER_READY. Boundary-layer, separation, wake, mesh-refinement and transient validation remain under #461.

No tolerance is relaxed and no validation is disabled.
