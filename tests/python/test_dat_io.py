from pathlib import Path

import pytest

from cfdx.dat_io import read_dat_restart


def test_read_dat_discovers_all_fields(tmp_path: Path) -> None:
    path = tmp_path / "state.dat"
    path.write_text(
        "CFDX-DAT 2\n"
        "cells 2\n"
        "iteration 42\n"
        "time 3.5\n"
        "field U 3\n"
        "1 2 3\n"
        "4 5 6\n"
        "field p 1\n"
        "10\n"
        "20\n"
        "field T 1\n"
        "300\n"
        "301\n"
        "field k 1\n"
        "0.1\n"
        "0.2\n"
        "field epsilon 1\n"
        "0.01\n"
        "0.02\n"
        "field omega 1\n"
        "1.0\n"
        "2.0\n",
        encoding="utf-8",
    )

    restart = read_dat_restart(path)

    assert restart.version == 2
    assert restart.cells == 2
    assert restart.iteration == 42
    assert restart.time == pytest.approx(3.5)
    assert set(restart.fields) == {"U", "p", "T", "k", "epsilon", "omega"}
    assert restart.fields["U"].dimension == 3
    assert restart.fields["U"].component(1) == [2.0, 5.0]
    assert restart.fields["T"].values == [300.0, 301.0]


def test_read_dat_accepts_current_v1_format(tmp_path: Path) -> None:
    path = tmp_path / "state.dat"
    path.write_text(
        "CFDX-DAT 1\n"
        "cells 1\n"
        "iteration 7\n"
        "time 1.25\n"
        "field U 3\n"
        "0.1 0.2 0.3\n"
        "field p 1\n"
        "12.0\n",
        encoding="utf-8",
    )

    restart = read_dat_restart(path)

    assert restart.version == 1
    assert set(restart.fields) == {"U", "p"}


@pytest.mark.parametrize(
    "payload",
    [
        "CFDX-DAT 2\ncells 1\niteration 0\ntime 0\n",
        "CFDX-DAT 2\ncells 1\niteration 0\ntime nan\nfield T 1\n1\n",
        "CFDX-DAT 2\ncells 1\niteration 0\ntime 0\nfield T 0\n1\n",
        "CFDX-DAT 2\ncells 1\niteration 0\ntime 0\nfield T 1\n1\nfield T 1\n2\n",
    ],
)
def test_read_dat_rejects_invalid_checkpoints(tmp_path: Path, payload: str) -> None:
    path = tmp_path / "invalid.dat"
    path.write_text(payload, encoding="utf-8")
    with pytest.raises(ValueError):
        read_dat_restart(path)


def test_hdf5_dat_roundtrip_all_fields(tmp_path: Path) -> None:
    from cfdx.dat_io import DatField, DatRestart, write_dat_hdf5

    restart = DatRestart(
        version=2,
        cells=2,
        iteration=12,
        time=4.5,
        fields={
            "U": DatField("U", 3, [1, 2, 3, 4, 5, 6]),
            "p": DatField("p", 1, [10, 11]),
            "T": DatField("T", 1, [300, 301]),
            "k": DatField("k", 1, [0.1, 0.2]),
            "omega": DatField("omega", 1, [2.0, 3.0]),
        },
    )
    path = tmp_path / "checkpoint.dat"
    write_dat_hdf5(path, restart)

    loaded = read_dat_restart(path)
    assert loaded == restart


def test_hdf5_checkpoint_persists_global_cell_ids_and_remaps() -> None:
    from cfdx.dat_io import DatField, DatRestart, remap_dat_restart, write_dat_hdf5

    restart = DatRestart(
        version=2,
        cells=4,
        iteration=9,
        time=1.5,
        fields={
            "p": DatField("p", 1, [10.0, 20.0, 30.0, 40.0]),
            "U": DatField("U", 3, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]),
        },
        cell_ids=(100, 7, 42, 900),
    )
    path = Path("checkpoint-ids.h5")
    write_dat_hdf5(path, restart)
    loaded = read_dat_restart(path)
    assert loaded.cell_ids == (100, 7, 42, 900)

    remapped = remap_dat_restart(loaded, (42, 100, 900))
    assert remapped.cells == 3
    assert remapped.cell_ids == (42, 100, 900)
    assert remapped.fields["p"].values == [30.0, 10.0, 40.0]
    assert remapped.fields["U"].values == [7, 8, 9, 1, 2, 3, 10, 11, 12]
    path.unlink()


def test_dat_checkpoint_rejects_duplicate_global_cell_ids(tmp_path: Path) -> None:
    from cfdx.dat_io import DatField, DatRestart, write_dat_hdf5

    restart = DatRestart(
        version=2,
        cells=2,
        iteration=0,
        time=0.0,
        fields={"p": DatField("p", 1, [1.0, 2.0])},
        cell_ids=(5, 5),
    )
    with pytest.raises(ValueError, match="duplicate cell ids"):
        write_dat_hdf5(tmp_path / "bad.h5", restart)
