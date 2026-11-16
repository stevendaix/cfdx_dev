"""Enable V87 pre-SST velocity screen on V86 BFS20540 trust case only.

No changes to geometry, BCs, viscosity, turbulence, numerical tolerances,
outer retry parameters, or final qualification gates. V87 is OFF by default.
"""

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import h5py
import numpy as np


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(preimage, destination):
    preimage = Path(preimage).resolve(strict=True)
    destination = Path(destination).absolute()
    if preimage == destination or destination.exists():
        raise ValueError("Output must be a NEW path; originals may not be overwritten")
    with h5py.File(preimage, "r") as f:
        if int(f.attrs.get("n_cells", 0)) != 20540:
            raise ValueError("Requires BFS20540 geometry")
        setup = json.loads(f.attrs["case_setup_json"])
        raw = setup["numerics"]["raw_settings"]
        if (
            raw.get("v86_bfs_physical_trust_region") != "enabled"
            or raw.get("v85_bounded_nonmonotone_outer_progress") != "enabled"
        ):
            raise ValueError("V87 requires V86 and V85 experimental options")
        if setup["numerics"]["max_iterations"] != 100:
            raise ValueError("Requires same 100-outer case")
    shutil.copy2(preimage, destination)
    try:
        with h5py.File(destination, "r+") as f:
            setup = json.loads(f.attrs["case_setup_json"])
            raw = setup["numerics"]["raw_settings"]
            if "v87_bfs_early_velocity_screen" in raw:
                raise ValueError("V87 option already present in preimage")
            raw["v87_bfs_early_velocity_screen"] = "enabled"
            f.attrs.modify("case_setup_json", json.dumps(setup, separators=(",", ":")))
        with h5py.File(preimage, "r") as a, h5py.File(destination, "r") as b:
            for key in a.attrs:
                if key == "case_setup_json":
                    continue
                if key not in b.attrs or not np.array_equal(a.attrs[key], b.attrs[key]):
                    raise RuntimeError("Unexpected modified metadata: " + key)

            def verify(name, obj):
                if isinstance(obj, h5py.Dataset) and (
                    name not in b or not np.array_equal(obj[...], b[name][...])
                ):
                    raise RuntimeError("Unexpected changed physical data: " + name)

            a.visititems(verify)
    except BaseException:
        destination.unlink(missing_ok=True)
        raise
    return {
        "input_sha256": sha(preimage),
        "output_sha256": sha(destination),
        "n_cells": 20540,
        "maximum_outer": 100,
        "change": "ONLY v87_bfs_early_velocity_screen=enabled",
        "physical_reference_qualified": False,
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    print(json.dumps(prepare(a.source, a.output), indent=2))
