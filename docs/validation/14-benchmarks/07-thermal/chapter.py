# %% [markdown]
# # Thermal benchmark family
#
# Thermal evidence is split into conduction, convection, CHT and coupled thermal-radiation populations.
#
# ## 1-D conduction
#
# \[
# \frac{d^2T}{dx^2}=0,
# \qquad
# T(x)=T_0+\frac{T_L-T_0}{L}x,
# \]
# \[
# q_x=-k\frac{T_L-T_0}{L}.
# \]
#
# Verify both temperature error and integrated heat-flux balance. Mesh refinement should recover the expected spatial order of the selected discretisation.
#
# ## Convection/energy balance
#
# For a fully developed heated channel,
# \[
# \dot m c_p(T_b-T_{b,0})=Q_{wall}.
# \]
# The bulk-temperature definition and wall heat-input convention must be frozen.
#
# ## CHT
#
# Perfect contact requires
# \[
# T_s=T_f,\qquad q_s''+q_f''=0.
# \]
# Finite contact conductance gives
# \[
# q''=G(T_1-T_2).
# \]
# Report interface temperature jump, heat flux on both sides, global energy imbalance and coupling-iteration error independently.
#
# ## Qualification boundary
#
# A converged temperature field is not sufficient. Property-model correctness, energy conservation, interface conservation and mesh/time sensitivity are separate evidence layers.
