"""BFS six-Schur software gate, selection and physically unqualified status."""
from __future__ import annotations

import copy
import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import h5py
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))
sys.path.insert(0, str(ROOT / "src/cfdx/python"))
from generate_th03_bfs_case import write_case  # noqa: E402
from qualify_bfs_schur_evidence import validate_bfs_evidence, strict_load  # noqa: E402
spec = importlib.util.spec_from_file_location(
    "bfs_v9_campaign", ROOT / "scripts/run_cfdx_coupled_campaign.py"
)
assert spec and spec.loader
campaign = importlib.util.module_from_spec(spec)
spec.loader.exec_module(campaign)


def _args(tmp_path):
    return SimpleNamespace(
        build=tmp_path / "build", meshes=tmp_path / "meshes", out=tmp_path / "out",
        only="bfs", models=[], n9_models=[], levels=[], iterations=10,
        bfs_variants=list(campaign.BFS_VARIANTS), wall_timeout=None,
    )


def test_six_schur_plus_block_preconditioners_have_real_bfs_commands(tmp_path):
    variants = list(campaign.BFS_VARIANTS)
    args = _args(tmp_path)
    commands = list(campaign.build_commands(args))
    assert len(commands) == len(variants) + 1  # quick smoke remains separate
    assert {t for t, _, _ in commands if t in {'bfs_' + v for v in variants}} == {
        'bfs_' + v for v in variants
    }
    fingerprints = set()
    selections = set()
    for tag, cmd, _ in commands:
        if not tag.removeprefix("bfs_") in campaign.BFS_VARIANTS:
            continue
        assert cmd[1] == "--mesh"
        assert cmd[-2:] == ["--experimental-sst", "--benchmark-export"]
        assert cmd[cmd.index("--iterations") + 1] == "10"
        with h5py.File(cmd[2]) as h5:
            fingerprints.add((h5.attrs["topology_hash"], h5.attrs["geometry_hash"]))
            setup = json.loads(h5.attrs["case_setup_json"])
            assert setup["numerics"]["max_iterations"] == 10  # not clipped to 250
            assert setup["numerics"]["raw_settings"]["v8_global_nonlinear_residual"] == "true"
            entries = setup["numerics"]["selection"]["entries"]
            selections.add(tuple((e["family"], e["configuration_key"]) for e in entries))
    assert len(fingerprints) == 1  # identical geometry, not just same cell count
    assert len(selections) == len(variants)
    assert {x[1] for item in selections for x in item if x[0] == 'schur'} == {
        f'schur.{name}' for name in campaign.SCHUR_MODELS
    }


@pytest.mark.parametrize("algorithm,schur,pc", [
    ("coupled", "pcd", "braess_sarazin"),
    ("simple", "pcd", None),
    ("coupled", "not_schur", "mgr"),
    ("coupled", "pcd", "native_amg"),
    ("coupled", None, "auto"),
])
def test_incompatible_selections_rejected_before_h5_write(tmp_path, algorithm, schur, pc):
    path = tmp_path / "must_not_exist.h5"
    with pytest.raises(ValueError):
        write_case(path, algorithm=algorithm, schur=schur, preconditioner=pc)
    assert not path.exists()


def _evidence(schur="pcd", pc="pcd"):
    v8 = dict(evaluated=True, algorithm="coupled", schur=schur,
              momentum_global_relative=1e-12,
              momentum_max_active_component_relative=1e-12,
              momentum_component_relative=[1e-12, 1e-12, 0.],
              continuity=1e-12, k=1e-12, omega=1e-12,
              eddy_viscosity_closure=1e-12, merit=1e-12,
              sst=True, sst_transport_committed=True, bounded=True,
              all_finite=True, converged=True)
    summary = {"format": "CFDX-EXECUTION", "process_exit_code": 0, "converged": True}
    history = {"format": "CFDX-CONVERGENCE", "converged": True, "iterations": 1,
               "linear_plan": {"schur_resolved": True, "requested_schur": schur,
                               "resolved_schur": schur,
                               "coupled": {"resolved": True, "krylov": "fgmres",
                                           "preconditioner": pc}},
               "history": [{"iteration": 1, "high_order_fraction": 1.0,
                            "full_order_eligible": True,
                            "mass_normalized_imbalance": 1e-12,
                            "global_nonlinear_v8": v8}]}
    return summary, history


def _write(tmp_path, pair):
    summary, history = pair
    (tmp_path / "execution.json").write_text(json.dumps(summary))
    (tmp_path / "convergence.json").write_text(json.dumps(history))


def _gate(tmp_path):
    return validate_bfs_evidence(tmp_path, 0, algorithm="coupled",
                                 schur="pcd", preconditioner="pcd")


def test_strict_realized_plan_and_v8_gate_never_promote_physics(tmp_path):
    _write(tmp_path, _evidence())
    got = _gate(tmp_path)
    assert got["execution_proven"] and got["nonlinear_v8_proven"]
    assert got["physical_validation"] != "PASS"
    assert got["measured_reattachment_x_over_h"] is None
    assert got["actual_schur"] == "pcd"
    assert got["actual_preconditioner"] == "pcd"


@pytest.mark.parametrize("mutate", [
    lambda s, h: h.pop("linear_plan"),
    lambda s, h: h["linear_plan"]["coupled"].update(preconditioner="mgr"),
    lambda s, h: h["linear_plan"].update(resolved_schur="lsc"),
    lambda s, h: h["history"][-1].pop("global_nonlinear_v8"),
    lambda s, h: h["history"][-1]["global_nonlinear_v8"].update(k=None),
    lambda s, h: h["history"][-1]["global_nonlinear_v8"].update(sst_transport_committed=False),
    lambda s, h: h["history"][-1]["global_nonlinear_v8"].update(merit=float("nan")),
    lambda s, h: h["history"][-1].update(high_order_fraction=.6),
    lambda s, h: s.update(converged=False),
])
def test_no_missing_mismatched_or_nonfinite_evidence_can_pass(tmp_path, mutate):
    pair = _evidence()
    mutate(*pair)
    _write(tmp_path, pair)
    assert not _gate(tmp_path)["nonlinear_v8_proven"]


def test_duplicate_json_keys_are_not_tolerated(tmp_path):
    path = tmp_path / "duplicates.json"
    path.write_text('{"converged": true, "converged": true}')
    with pytest.raises(ValueError, match="duplicate"):
        strict_load(path)


def test_execution_status_not_proven_from_just_exit_zero(tmp_path):
    assert not _gate(tmp_path)["nonlinear_v8_proven"]


def test_imported_reference_candidate_keeps_boundary_and_field_datasets(tmp_path):
    from generate_th03_bfs_case import reconfigure_bfs_reference_case
    source = tmp_path / 'reference.h5'
    write_case(source, n_up=5, n_down=20, ny=12, step_cells=3,
               algorithm='simple', iterations=12)
    output = tmp_path / 'variant.h5'
    result = reconfigure_bfs_reference_case(
        source, output, algorithm='coupled', schur='lsc',
        preconditioner='coupled_block_schur', iterations=23)
    assert result['physical_qualification'] == 'NOT_QUALIFIED'
    assert result['input_kind'] == 'EXTERNAL_BFS_REFERENCE_CANDIDATE'
    with h5py.File(source) as before, h5py.File(output) as after:
        for key in before.keys():
            assert (before[key][()] == after[key][()]).all()
        assert before.attrs['geometry_hash'] == after.attrs['geometry_hash']
        assert before.attrs['topology_hash'] == after.attrs['topology_hash']
        a = json.loads(before.attrs['case_setup_json'])
        b = json.loads(after.attrs['case_setup_json'])
        assert a['boundary_conditions'] == b['boundary_conditions']
        assert a['initial_condition'] == b['initial_condition']
        assert a['materials'] == b['materials']
        entries = b['numerics']['selection']['entries']
        assert any(e['configuration_key'] == 'schur.lsc' for e in entries)
        assert any(e['configuration_key'] == 'preconditioner.coupled_block_schur'
                   for e in entries)
        assert b['numerics']['max_iterations'] == 23


def test_case_sha_input_and_numerical_selection_are_checked(tmp_path):
    from hashlib import sha256
    case_path = tmp_path / 'bfs.h5'
    write_case(case_path, algorithm='coupled', schur='pcd',
               preconditioner='pcd')
    _write(tmp_path, _evidence())
    digest = sha256(case_path.read_bytes()).hexdigest()
    kw = dict(algorithm='coupled', schur='pcd', preconditioner='pcd',
              case_path=case_path, case_sha256=digest)
    result = validate_bfs_evidence(tmp_path, 0, **kw)
    assert result['nonlinear_v8_proven']
    assert 'case_sha256_mesh_and_numerics_verified' in result['checks']
    assert not validate_bfs_evidence(tmp_path, 0, **{**kw, 'case_sha256': 'fake'})[
        'nonlinear_v8_proven']
    assert not validate_bfs_evidence(tmp_path, 0, **{**kw, 'preconditioner': 'mgr'})[
        'nonlinear_v8_proven']


def test_upwind_cannot_pass_second_order_gate_even_with_v8_certificate(tmp_path):
    from hashlib import sha256
    case_path = tmp_path / 'upwind.h5'
    write_case(case_path, algorithm='coupled', schur='pcd', preconditioner='pcd')
    with h5py.File(case_path, 'r+') as h5:
        setup = json.loads(h5.attrs['case_setup_json'])
        for entry in setup['numerics']['selection']['entries']:
            if entry['family'] == 'convection':
                entry['configuration_key'] = 'numerics.convection.upwind'
        h5.attrs['case_setup_json'] = json.dumps(setup)
    _write(tmp_path, _evidence())
    gate = validate_bfs_evidence(tmp_path, 0, algorithm='coupled', schur='pcd',
                                 preconditioner='pcd', case_path=case_path,
                                 case_sha256=sha256(case_path.read_bytes()).hexdigest())
    assert not gate['nonlinear_v8_proven']
    assert 'case_second_order_convection_not_selected' in gate['defects']


def test_invalid_segregated_resolved_plan_is_fail_closed(tmp_path):
    summary, history = _evidence()
    v8 = history['history'][0]['global_nonlinear_v8']
    v8['algorithm'] = 'simple'
    v8['schur'] = None
    history['linear_plan'] = {'coupled': None}
    _write(tmp_path, (summary, history))
    result = validate_bfs_evidence(tmp_path, 0, algorithm='simple',
                                   schur=None, preconditioner=None)
    assert result['nonlinear_v8_proven']
    history['linear_plan']['coupled'] = {'resolved': True}
    _write(tmp_path, (summary, history))
    result = validate_bfs_evidence(tmp_path, 0, algorithm='simple',
                                   schur=None, preconditioner=None)
    assert not result['nonlinear_v8_proven']
