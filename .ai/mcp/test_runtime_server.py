from __future__ import annotations

import asyncio
import pathlib
import sys
import tempfile

import h5py
from mcp import Client

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from runtime_server import create_server


def make_artifacts(root: pathlib.Path) -> None:
    with h5py.File(root / "demo.cfdx.h5", "w") as h5:
        h5.attrs["format"] = "CFDX"
        h5.attrs["schema_version"] = 1
        h5.attrs["dimension"] = 3
        h5.attrs["precision"] = "float64"
        h5.attrs["case_revision"] = 4
        h5.create_dataset("case/config", data='{"name":"demo","physics":{"laminar":{"enabled":true}},"numerics":{"scheme":"simple"},"boundaries":{},"materials":{},"execution":{"policy":"AUTO","solver":null,"mpi_ranks":1,"deterministic":true,"restart_option":"--restart"}}')
        h5.create_dataset("points", data=[[0.0, 0.0, 0.0]])
        h5.create_dataset("fields/values", data=[[1.0]])

    with h5py.File(root / "demo.dat.h5", "w") as h5:
        h5.attrs["format"] = "CFDX-DAT"
        h5.attrs["version"] = 2
        h5.attrs["cells"] = 1
        h5.attrs["iteration"] = 42
        h5.attrs["time"] = 0.25
        h5.create_dataset("cell_ids", data=[7])
        h5.create_dataset("fields/p", data=[1.5])


async def exercise() -> None:
    with tempfile.TemporaryDirectory() as directory:
        root = pathlib.Path(directory)
        make_artifacts(root)
        client = Client(create_server(str(root)))
        async with client:
            listed = await client.list_tools()
            assert {tool.name for tool in listed.tools} == {"case.inspect", "case.validate", "checkpoint.inspect", "checkpoint.field.inspect", "checkpoint.compare"}
            for tool in listed.tools:
                annotations = tool.model_dump(by_alias=True).get("annotations", {})
                assert annotations["readOnlyHint"] is True
                assert annotations["openWorldHint"] is False

            case = await client.call_tool("case.inspect", {"case_path": "demo.cfdx.h5"})
            assert case.is_error is False
            data = case.structured_content
            assert data["ok"] is True
            assert data["artifact"] == "case"
            assert data["configuration"]["name"] == "demo"
            assert data["attributes"]["case_revision"] == 4
            assert data["runtime_state_present"] is False

            case_validation = await client.call_tool("case.validate", {"case_path": "demo.cfdx.h5"})
            assert case_validation.is_error is False
            assert case_validation.structured_content["ok"] is True
            assert case_validation.structured_content["artifact"] == "case-validation"
            assert case_validation.structured_content["checks"]["mesh"] is True
            assert case_validation.structured_content["checks"]["runtime_separated"] is True

            checkpoint = await client.call_tool("checkpoint.inspect", {"checkpoint_path": "demo.dat.h5"})
            assert checkpoint.is_error is False
            data = checkpoint.structured_content
            assert data["ok"] is True
            assert data["artifact"] == "checkpoint"
            assert data["attributes"]["iteration"] == 42
            assert data["attributes"]["time"] == 0.25
            assert data["metadata"] == {"version": 2, "cells": 1, "iteration": 42, "time": 0.25}
            assert data["cell_ids_present"] is True
            assert data["fields"][0]["name"] == "p"

            with h5py.File(root / "invalid.dat.h5", "w") as h5:
                h5.attrs["format"] = "CFDX-DAT"
                h5.attrs["version"] = 2
                h5.attrs["cells"] = 2
                h5.attrs["iteration"] = 1
                h5.attrs["time"] = 0.1
                h5.create_dataset("cell_ids", data=[1, 1])
                h5.create_dataset("fields/p", data=[1.0])
            invalid = await client.call_tool("checkpoint.inspect", {"checkpoint_path": "invalid.dat.h5"})
            assert invalid.is_error is False
            assert invalid.structured_content["ok"] is False
            assert "duplicate cell ids" in invalid.structured_content["errors"][0]

            scalar = await client.call_tool(
                "checkpoint.field.inspect",
                {"checkpoint_path": "demo.dat.h5", "field_name": "p"},
            )
            assert scalar.is_error is False
            data = scalar.structured_content
            assert data["ok"] is True
            assert data["dimension"] == 1
            assert data["cell_count"] == 1
            assert data["finite"] is True
            assert data["min"] == 1.5
            assert data["max"] == 1.5

            vector = await client.call_tool(
                "checkpoint.field.inspect",
                {"checkpoint_path": "demo.dat.h5", "field_name": "U", "component": 1},
            )
            assert vector.structured_content["ok"] is False

            with h5py.File(root / "demo.dat.h5", "a") as h5:
                h5.create_dataset("fields/U", data=[[1.0, -2.0, 3.0]])
            vector = await client.call_tool(
                "checkpoint.field.inspect",
                {"checkpoint_path": "demo.dat.h5", "field_name": "U", "component": 1},
            )
            assert vector.is_error is False
            data = vector.structured_content
            assert data["dimension"] == 3
            assert data["component"] == 1
            assert data["min"] == -2.0
            assert data["max"] == -2.0

            with h5py.File(root / "later.dat.h5", "w") as h5:
                h5.attrs["format"] = "CFDX-DAT"
                h5.create_group("fields").create_dataset("p", data=[2.0])
            comparison = await client.call_tool(
                "checkpoint.compare",
                {"left_checkpoint_path": "demo.dat.h5", "right_checkpoint_path": "later.dat.h5", "field_name": "p"},
            )
            assert comparison.is_error is False
            data = comparison.structured_content
            assert data["ok"] is True
            assert data["value_count"] == 1
            assert data["max_abs_difference"] == 0.5
            assert data["relative_l2_difference"] == 1.0 / 3.0

            for name, args in (
                ("case.inspect", {"case_path": "../outside.cfdx.h5"}),
                ("checkpoint.inspect", {"checkpoint_path": "/etc/passwd"}),
                ("checkpoint.field.inspect", {"checkpoint_path": "../outside.dat.h5", "field_name": "p"}),
                ("checkpoint.compare", {"left_checkpoint_path": "../a.dat.h5", "right_checkpoint_path": "later.dat.h5", "field_name": "p"}),
                ("checkpoint.inspect", {"checkpoint_path": "/etc/passwd"}),
                ("checkpoint.field.inspect", {"checkpoint_path": "../outside.dat.h5", "field_name": "p"}),
            ):
                result = await client.call_tool(name, args)
                assert result.is_error is False
                assert result.structured_content["ok"] is False


if __name__ == "__main__":
    asyncio.run(exercise())
