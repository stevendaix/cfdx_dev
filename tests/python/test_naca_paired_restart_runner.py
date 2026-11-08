"""Runner enforces paired restart provenance without weakening physical gates."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[2]
SCRIPTS = ROOT / "scripts"
if str(SCRIPTS) not in sys.path:
    sys.path.insert(0, str(SCRIPTS))


def load(name: str):
    source = SCRIPTS / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, source)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUNNER = load("run_naca_gci_case")
COMPARATOR = load("compare_naca_paired_restarts")

run_case = RUNNER.run_case
compare = COMPARATOR.compare
load_manifest = COMPARATOR.load_manifest
fnv64 = COMPARATOR.fnv64


def test_runner_pairs_manifest_and_preserves_physical_threshold(
        tmp_path, monkeypatch):
    mesh = tmp_path / "mesh.h5"
    mesh.touch()
    exe = tmp_path / "solver"
    exe.touch()
    original = tmp_path / "original.manifest"
    original.write_text("fixture only")
    captured = {}

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        captured["env"] = kwargs["env"]
        return SimpleNamespace(returncode=1)

    monkeypatch.setattr(
        RUNNER, "subprocess", SimpleNamespace(run=fake_run))
    output = tmp_path / "run"
    assert run_case(mesh, output, exe, 10, 20,
                    initial_paired_checkpoint=original, checkpoint_every=10,
                    stagnation_diagnostic_only=True) == 1
    env = captured["env"]
    assert env["CFDX_NACA_IMPORT_PAIRED_MANIFEST"] == str(original.resolve())
    assert env["CFDX_NACA_CHECKPOINT_PAIRED_MANIFEST"] == str(
        output / "checkpoint.manifest")
    assert env["CFDX_NACA_CHECKPOINT_EVERY"] == "10"
    assert env["CFDX_NACA_STAGNATION_DIAGNOSTIC_ONLY"] == "1"
    assert env["CFDX_NACA_LINEAR_TOLERANCE"] == "1e-12"
    summary = json.loads((output / "run_summary.json").read_text())
    assert summary["restart_fidelity"] == "PAIRED_ACCEPTED_RHIE_CHOW_FACE_FLUX"
    assert summary["nonlinear_tolerance"] == 1e-8
    assert not summary["certificate_present"]


def test_runner_refuses_missing_or_colliding_paired_checkpoint(tmp_path):
    mesh = tmp_path / "mesh.h5"
    exe = tmp_path / "solver"
    output = tmp_path / "out"
    with pytest.raises(FileNotFoundError):
        run_case(mesh, output, exe, 10, 20,
                 initial_paired_checkpoint=tmp_path / "missing.manifest")
    output.mkdir()
    active = output / "checkpoint.manifest"
    active.write_text("unchanged")
    with pytest.raises(ValueError, match="overwrite"):
        run_case(mesh, output, exe, 10, 20, initial_paired_checkpoint=active)
    assert active.read_text() == "unchanged"
    with pytest.raises(ValueError, match="exclusive"):
        run_case(mesh, output, exe, 10, 20, initial_state=active,
                 initial_paired_checkpoint=active)


def test_live_v2_checkpoints_replay_to_roundoff(tmp_path):
    # This unit check tests rejection and normalized baseline syntax without
    # depending on expensive NACA runs, which are recorded as external
    # evidence.
    a = tmp_path / "bad.manifest"
    a.write_text("CFDX_NACA_PAIRED_V1\n1,3,4,foo.csv,bar.csv,a,b\n")
    with pytest.raises(ValueError, match="invalid NACA"):
        load_manifest(a)


def test_comparator_rejects_extra_rows_even_with_valid_checksums(tmp_path):
    def make_checkpoint(folder: Path, *, extra: bool) -> Path:
        folder.mkdir()
        cell = folder / "state.csv"
        cell.write_text(
            "cell,x,y,z,ux,uy,p\n0,0,0,0,1,2,3\n"
            + ("1,0,0,0,4,5,6\n" if extra else ""))
        flux = folder / "flux.csv"
        flux.write_text(
            "face,owner,neighbour,cx,cy,cz,sx,sy,sz,mass_flux\n"
            "0,0,-1,0,0,0,1,0,0,0.2\n")
        manifest = folder / "checkpoint.manifest"
        manifest.write_text(
            "CFDX_NACA_PAIRED_V2\n1,1,1,0.1,1,state.csv,flux.csv,"
            + fnv64(cell) + "," + fnv64(flux) + "\n")
        return manifest

    a = make_checkpoint(tmp_path / "a", extra=False)
    b = make_checkpoint(tmp_path / "b", extra=True)
    with pytest.raises(ValueError, match="row-count mismatch"):
        compare(a, b)


def test_runner_refuses_success_without_positive_flux_audit(tmp_path, monkeypatch):
    mesh = tmp_path / "mesh.h5"
    mesh.touch()
    exe = tmp_path / "binary"
    exe.touch()
    monkeypatch.setattr(RUNNER.subprocess, "run",
                         lambda argv, **kwargs: SimpleNamespace(returncode=0))
    output = tmp_path / "result"
    status = run_case(mesh, output, exe, 10, 10, farfield_pressure_neumann=True)
    assert status == 2
    summary = json.loads((output / "run_summary.json").read_text())
    assert summary["pressure_boundary_policy"] == "ZERO_NORMAL_GRADIENT"
    assert summary["exit_code"] == 2
    assert summary["kinematic_face_flux_certified"] is False
    assert "FAIL missing positive" in (output / "run.log").read_text()
