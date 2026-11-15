"""Enable V86 BFS-only field trust region without changing FV equations or mesh.

The parameters 2*Uin and .03*Uin^2 bound *provisional* candidates and
trigger existing nonlinear rollback/retries. They do not certify NASA physics.
"""

import argparse
import hashlib
import json
import shutil
from pathlib import Path

import h5py
import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def prepare(preimage: Path, destination: Path):
    preimage = Path(preimage).resolve(strict=True)
    destination = Path(destination).absolute()
    if preimage == destination or destination.exists():
        raise ValueError("Output must be a new case")
    with h5py.File(preimage, "r") as h:
        if int(h.attrs.get("n_cells", 0)) != 20540:
            raise ValueError("BFS20540 required")
        s = json.loads(h.attrs["case_setup_json"])
        r = s["numerics"]["raw_settings"]
        if (
            s["numerics"]["max_iterations"] != 100
            or r.get("v85_bounded_nonmonotone_outer_progress") != "enabled"
        ):
            raise ValueError("Must start from 100-outer V85 provisional SIMPLE case")
        inlet = [
            bc for bc in s["boundary_conditions"] if bc.get("patch_name") == "inlet"
        ]
        if len(inlet) != 1 or inlet[0]["velocity_vector"] != [44.2, 0.0, 0.0]:
            raise ValueError(
                "Case-specific trust region requires exact TH03 inlet speed"
            )
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(preimage, destination)
    try:
        with h5py.File(destination, "r+") as h:
            s = json.loads(h.attrs["case_setup_json"])
            s["numerics"]["raw_settings"]["v86_bfs_physical_trust_region"] = "enabled"
            h.attrs.modify("case_setup_json", json.dumps(s, separators=(",", ":")))
        with h5py.File(preimage) as a, h5py.File(destination) as b:
            for attr in (
                "topology_hash",
                "mesh_hash",
                "geometry_hash",
                "boundary_patches",
            ):
                if a.attrs.get(attr) != b.attrs.get(attr):
                    raise RuntimeError("Changed mesh/provenance attr " + attr)

            def check(name, obj):
                if isinstance(obj, h5py.Dataset):
                    if name not in b or not np.array_equal(obj[...], b[name][...]):
                        raise RuntimeError("Altered physical data " + name)

            a.visititems(check)
    except BaseException:
        destination.unlink(missing_ok=True)
        raise
    return {
        "source_sha256": digest(preimage),
        "output_sha256": digest(destination),
        "n_cells": 20540,
        "max_outer": 100,
        "only_change": "v86_bfs_physical_trust_region=enabled",
        "geometry_qualification": "RECONSTRUCTED_SIX_BLOCK_ONLY",
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("source", type=Path)
    p.add_argument("output", type=Path)
    a = p.parse_args()
    print(json.dumps(prepare(a.source, a.output), indent=2))
