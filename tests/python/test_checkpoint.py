from pathlib import Path

import pytest

from cfdx.checkpoint import load_checkpoint, save_checkpoint
from cfdx import Checkpoint


def test_checkpoint_roundtrip(tmp_path: Path) -> None:
    original = Checkpoint(2, 3, 4, 5, 120, 1.25)
    path = tmp_path / "checkpoint.json"
    save_checkpoint(original, path)
    assert load_checkpoint(path) == original


def test_checkpoint_schema_is_strict(tmp_path: Path) -> None:
    path = tmp_path / "bad.json"
    path.write_text("{}")
    with pytest.raises(ValueError, match="fields"):
        load_checkpoint(path)
