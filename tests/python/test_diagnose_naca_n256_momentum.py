"""Counterexamples for the NACA same-state residual / farfield audit."""
import importlib.util
from pathlib import Path

import pytest

MODULE = Path(__file__).resolve().parents[2] / "scripts" / "diagnose_naca_n256_momentum.py"
spec = importlib.util.spec_from_file_location("naca_n256_diagnosis", MODULE)
mod = importlib.util.module_from_spec(spec)
assert spec and spec.loader
spec.loader.exec_module(mod)


def test_mass_flux_kinematically_valid_but_forced_wake():
    cells = [{"cell":"0","ux":".94","uy":"0"}]
    faces = [{"face":"9","owner":"0","neighbour":"-1","cx":"20","cy":"0",
              "sx":".25","sy":"0","sz":"0","mass_flux":".25"}]
    a = mod.mass_flux_audit(cells, faces)
    assert a["boundary_flux_contract_error"] == 0
    assert a["wake_to_prescribed_flux_gap"] == pytest.approx(.015)
    assert a["wake_to_prescribed_relative_gap"] == pytest.approx(.06)


def test_refuse_mixed_iteration_and_malformed_history():
    history = [{"iteration":"1"},{"iteration":"2"}]
    budget = [{"iteration":"1"}]
    with pytest.raises(ValueError, match="different accepted iterations"):
        mod.iteration_rows(history,budget)
    with pytest.raises(ValueError, match="not strictly increasing"):
        mod.iteration_rows([history[1],history[0]],budget)


def test_mass_correction_is_not_the_physical_residual():
    b={"iteration":"10","component":"0","cell":"0","x":"20","y":"0",
       "physical_residual":"-4e-8","bounded_operator_residual":"-8e-9",
       "bounded_divergence":"3.2e-8","mass_divergence":str(-3.2e-8/.8),
       "convection":"-.003","pressure":".003","viscous":"0","transpose_viscous":"0"}
    out = mod.worst_budget([b],0,[{"cell":"0","ux":".8"}])
    assert out["fraction_physical_from_mass"] == pytest.approx(.8)
    bad = {**b,"bounded_operator_residual":"-4e-8"}
    with pytest.raises(ValueError,match="algebra"):
        mod.worst_budget([bad],0,[{"cell":"0","ux":".8"}])


def test_decay_fit_requires_improving_residual_and_many_steps():
    hist=[{"iteration":str(i),"res":str(1e-4 * .99**i)} for i in range(1,161)]
    out=mod.log_decay(hist,"res")
    assert out["fit_r2"] > .999999
    assert out["hypothetical_additional_iterations"] > 0
    bad=[{"iteration":str(i),"res":"1e-5"} for i in range(1,161)]
    assert mod.log_decay(bad,"res")["hypothetical_additional_iterations"] is None


def test_invalid_duplicate_cell_and_nonfinite_flux():
    cells=[{"cell":"0","ux":"1","uy":"0"}]
    f=[{"face":"0","owner":"0","neighbour":"-1","cx":"20","cy":"0",
        "sx":".5","sy":"0","sz":"0","mass_flux":"NaN"}]
    with pytest.raises(ValueError,match="nonfinite"):
        mod.mass_flux_audit(cells,f)
    f[0]["mass_flux"]=".5"
    with pytest.raises(ValueError,match="duplicated"):
        mod.mass_flux_audit(cells+cells,f)


def test_manifest_requires_matching_cell_flux_pair(tmp_path):
    manifest=tmp_path/"checkpoint.manifest"
    manifest.write_text("CFDX_NACA_PAIRED_V3\n675,1,1,0.1,1,state.csv,flux.csv,a,b,c\n")
    out=mod.checkpoint_metadata(manifest,tmp_path/"state.csv",tmp_path/"flux.csv",1,1)
    assert out["accepted_iteration_absolute"] == 675
    with pytest.raises(ValueError,match="does not match"):
        mod.checkpoint_metadata(manifest,tmp_path/"wrong.csv",tmp_path/"flux.csv",1,1)
    with pytest.raises(ValueError,match="does not match"):
        mod.checkpoint_metadata(manifest,tmp_path/"state.csv",tmp_path/"flux.csv",1,2)


def test_spatial_signature_detects_mirror_pair():
    rows=[]
    for i,(y,r) in enumerate(((-.5,1.1e-8),(.5,-1.1e-8),(.75,-2e-9))):
        rows.append({"component":"1","cell":str(i),"x":"19.8","y":str(y),
                     "physical_residual":str(r)})
    sig=mod.spatial_residual_signature(rows,1)
    assert sig["count_at_x_ge_19_5"] == 3
    assert sig["worst_two_mirror_symmetric"]

