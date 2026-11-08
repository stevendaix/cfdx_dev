"""Strict GUI/CLI hand-off to the C++ incompressible production solver.

The GUI HDF5 /case/config document and the native root case_setup_json
attribute are distinct schemas. They are intentionally *not* conflated.
A future translator must synchronise them before editable dual-schema cases
can be executed via the native solver.
"""

from __future__ import annotations

import json
from pathlib import Path

import h5py


def is_native_production_solver(executable: str | None) -> bool:
    if not executable:
        return False
    return Path(executable).name.lower() in {
        "cfdx_production_solver",
        "cfdx_production_solver.exe",
    }


def check_native_production_case(case_path: Path) -> None:
    """Fail before run when boundary settings cannot reach native runtime.

    An editable GUI document and a native case may coexist, but without a
    translation/synchronisation contract we cannot prove the GUI values are
    the values consumed by the C++ solver. Reject that ambiguous state.
    """
    path = Path(case_path)
    if not path.name.lower().endswith(".cfdx.h5"):
        raise ValueError("le solveur natif exige un cas .cfdx.h5")
    with h5py.File(path, "r") as h5:
        raw = h5.attrs.get("case_setup_json")
        if raw is None:
            raise ValueError(
                "case_setup_json absent : la sauvegarde GUI /case/config "
                "ne renseigne pas les conditions limites du solveur natif. "
                "Exporter un cas de production au schema natif avant Run."
            )
        if isinstance(raw, bytes):
            raw = raw.decode("utf-8")
        try:
            setup = json.loads(raw)
        except (TypeError, ValueError) as exc:
            raise ValueError("case_setup_json du solveur natif est invalide") from exc
        if not isinstance(setup, dict) or not isinstance(
            setup.get("boundary_conditions"), list
        ):
            raise TypeError("case_setup_json doit declarer boundary_conditions")
        if not setup["boundary_conditions"]:
            raise ValueError("case_setup_json ne contient aucune condition limite")
        if "case/config" in h5:
            raw_gui = h5["case/config"][()]
            if isinstance(raw_gui, bytes):
                raw_gui = raw_gui.decode("utf-8")
            try:
                gui_config = json.loads(raw_gui)
            except (TypeError, ValueError) as exc:
                raise ValueError("configuration GUI /case/config invalide") from exc
            if not isinstance(gui_config, dict):
                raise ValueError("configuration GUI /case/config non reconnue")
            if gui_config.get("boundaries"):
                raise ValueError(
                    "conditions limites GUI et natives coexistent sans "
                    "synchronisation prouvee : une modification dans le GUI "
                    "pourrait ne pas atteindre le solveur. Export natif requis."
                )


def native_production_command(
    executable: str,
    case_path: Path,
    *,
    output_dir: str | Path,
    max_iterations: int,
) -> list[str]:
    """Use the actual native C++ CLI, never a positional GUI case argument."""
    if not is_native_production_solver(executable):
        raise ValueError("executable different de cfdx_production_solver")
    if (
        isinstance(max_iterations, bool)
        or not isinstance(max_iterations, int)
        or max_iterations <= 0
    ):
        raise ValueError("max_iterations doit etre un entier strictement positif")
    if not str(output_dir).strip():
        raise ValueError("output_dir du solveur natif est vide")
    check_native_production_case(case_path)
    return [
        executable,
        "--mesh",
        str(case_path),
        "--output-dir",
        str(output_dir),
        "--iterations",
        str(max_iterations),
    ]
