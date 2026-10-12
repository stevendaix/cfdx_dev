from pathlib import Path
import json

import pytest

from cfdx.checkpoint import CHECKPOINT_SCHEMA_VERSION, load_checkpoint, save_checkpoint
from cfdx import Checkpoint


def test_checkpoint_roundtrip(tmp_path: Path) -> None:
    original = Checkpoint(2, 3, 4, 5, 120, 1.25)
    path = tmp_path / "checkpoint.json"
    save_checkpoint(original, path)
    assert load_checkpoint(path) == original
    assert json.loads(path.read_text())["schema_version"] == CHECKPOINT_SCHEMA_VERSION


def test_checkpoint_schema_is_strict(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="fields"):
        load_checkpoint(path)


def test_checkpoint_rejects_unsupported_schema_version(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({
        "schema_version": 99, "case_revision": 0, "mesh_revision": 0,
        "physics_revision": 0, "numerics_revision": 0, "iteration": 0, "time": 0.0
    }))
    with pytest.raises(ValueError, match="schema version"):
        load_checkpoint(path)


@pytest.mark.parametrize("field", ["case_revision", "mesh_revision", "physics_revision",
                                   "numerics_revision", "iteration"])
def test_checkpoint_rejects_negative_integer_metadata(tmp_path: Path, field: str) -> None:
    payload = {
        "schema_version": CHECKPOINT_SCHEMA_VERSION, "case_revision": 0,
        "mesh_revision": 0, "physics_revision": 0, "numerics_revision": 0,
        "iteration": 0, "time": 0.0,
    }
    payload[field] = -1
    path = tmp_path / "bad.json"
    path.write_text(json.dumps(payload))
    with pytest.raises(ValueError, match="non-negative"):
        load_checkpoint(path)


@pytest.mark.parametrize("time_value", ["nan", "inf", "-inf"])
def test_checkpoint_rejects_non_finite_time(tmp_path: Path, time_value: str) -> None:
    path = tmp_path / "bad.json"
    path.write_text(json.dumps({
        "schema_version": CHECKPOINT_SCHEMA_VERSION, "case_revision": 0,
        "mesh_revision": 0, "physics_revision": 0, "numerics_revision": 0,
        "iteration": 0, "time": time_value,
    }))
    with pytest.raises(ValueError, match="time"):
        load_checkpoint(path)


def test_checkpoint_publish_failure_preserves_previous_checkpoint(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    path = tmp_path / "checkpoint.json"
    original = Checkpoint(1, 2, 3, 4, 5, 0.5)
    save_checkpoint(original, path)

    def fail_replace(_source: str, _destination: Path) -> None:
        raise OSError("simulated publish failure")

    monkeypatch.setattr("cfdx.checkpoint.os.replace", fail_replace)
    with pytest.raises(OSError, match="publish"):
        save_checkpoint(Checkpoint(6, 7, 8, 9, 10, 1.0), path)

    assert load_checkpoint(path) == original
    assert not list(tmp_path.glob(".checkpoint.json.*.tmp"))
