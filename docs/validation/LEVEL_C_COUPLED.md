# Level C — coupled verification

Level C verifies interaction between the M1-M4 components.

Gates currently included:
- turbulence transport preserves positive k and epsilon;
- radiation and energy converge through the radiative source term;
- two-region CHT converges and drives interface heat-flux imbalance below tolerance.

These are integration consistency tests, not analytical solutions. Production-level extension is multi-cell coupled mesh/time-step convergence.
