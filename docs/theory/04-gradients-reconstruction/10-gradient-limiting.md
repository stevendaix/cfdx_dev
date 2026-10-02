# 10. Gradient Limiting

High-order reconstruction can create nonphysical extrema near discontinuities or strong gradients. A limiter modifies the reconstructed gradient so that extrapolated values satisfy local admissibility constraints.

A limiter therefore introduces a deliberate trade-off between formal smooth-region accuracy and boundedness near non-smooth regions.

The documentation must distinguish the unlimited gradient, limiter activation criterion, limiter coefficient, admissibility constraints and smooth-region preservation.

A limiter must not be silently enabled in a verification experiment intended to measure the formal order of the underlying gradient reconstruction.
