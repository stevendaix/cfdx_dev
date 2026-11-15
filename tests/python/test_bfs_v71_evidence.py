import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from audit_bfs_v71_joint import audit

hist = "\n".join([
    "GLOBAL_SAME_OPERATOR_TRIAL_V65 iteration=23 accepted=0",
    "GLOBAL_SAME_OPERATOR_COMMITTED_V65 iteration=22 committed=1",
    "BFS_SST_RC_PICARD_V69 iteration=23 private_transport=1",
    "INCOMPRESSIBLE_ITER iteration=22 converged=0",
]) + "\n"
valid = ("BFS_JOINT_PICARD_V71 iteration=23 retry=8 gamma=0.0 pressure_fraction=0.5 "
         "mass_before=1 mass_after_predictor=1 mass_after_frozen_pressure=0 "
         "mass_after_sst=0 merit_ratio=0.95 worst_growth=1.1 omega_after=0 "
         "physical_mass_closed=1 original_guards_ok=1 physically_eligible=1 "
         "private_trial=1 committed=0 qualified=0")
face = ("BFS_JOINT_FACE3_V71 iteration=23 gamma=0.0 pressure_fraction=0.5 "
        "phi3_before=0 phi3_after=1 rAUx_before=1 rAUx_after=1 "
        "nut_before=0.1 nut_after=0.1 omega_raw_before=1 omega_raw_after=0 "
        "convection_delta=-1 diffusion_delta=0 production_delta=0 "
        "destruction_delta=0 cross_delta=0 committed=0 qualified=0 cell=")
log = hist + valid + "\n" + face + "0\n" + face + "12\n"
r = audit(hist, log)
assert r["joint_trial_attempts"] == 1
assert r["face_cell_budget_rows"] == 2
assert not r["physical_convergence_qualified"]

def reject(candidate):
    try:
        audit(hist, candidate)
    except ValueError:
        return
    raise AssertionError("falsified V71 certificate was accepted")

reject(log.replace("committed=0", "committed=1"))
reject(log.replace("qualified=0", "qualified=1"))
reject(log.replace("phi3_after=1", "phi3_after=NaN"))
reject(log.replace("rAUx_before=1", "rAUx_before=inf"))
reject(log.replace("convection_delta=-1", "convection_delta=2"))
reject(log.replace("cell=12", "cell=0"))
reject(log.replace("cell=12", "cell=13"))
reject(log.replace("GLOBAL_SAME_OPERATOR_TRIAL_V65", "CHANGED_GLOBAL"))
reject(log.replace("physically_eligible=1", "physically_eligible=1 physical_mass_closed=0"))
reject(hist)
print("V71_FORENSIC_FALSIFICATION_PASS cases=10")