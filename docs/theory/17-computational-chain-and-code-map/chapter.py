# %% [markdown]
"""# CFDX Computational Chain and Code Map

The traceability chain is case -> schema/import -> topology -> geometry -> fields/BC -> numerical selection -> FVM operators -> physics -> linearisation -> sparse system -> solver/preconditioner -> updated fields -> diagnostics -> V&V. Each arrow is a contract for data, units, orientation or mathematics.
"""

# %%
CHAIN = ('case', 'mesh', 'geometry', 'fields', 'numerics', 'physics', 'linear_algebra', 'diagnostics', 'V&V')
assert CHAIN[0] == 'case' and CHAIN[-1] == 'V&V'
