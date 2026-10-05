from __future__ import annotations

import asyncio
import pathlib
import tempfile

import h5py
from mcp import Client

import sys\nsys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))\n\nfrom runtime_server import create_server


def make_artifacts(root: pathlib.Path) -> None:
    case = root / "demo.cfdx.h5"
    with h5py.File(case, "w") as h5:
        h5.attrs["format"] = "CFDX"
        h5.attrs["schema_version"] = 1
        h5.attrs["dimension"] = 3
        h5.attrs["precision"] = "float64"
        h5.attrs["case_revision"] = 4
        h5.create_dataset(
            "case/config",
            data='{"name":"demo","physics":{"laminar":{"enabled":true}},"numerics":{"scheme":"simple"},"boundaries":{},"materials":{},"execution":{"policy":"AUTO","solver":null,"mpi_ranks":1,"deterministic":true,"restart_option":"--restart"}}',
        )
        h5.create_dataset("points", data=[[0.0, 0.0, 0.0]])
        h5.create_dataset("fields/values", data=[[1.0]])

    checkpoint = root / "demo.dat.h5"
    with h5py.File(checkpoint, "w") as h5:
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
            names = {tool.name for tool in listed.tools}
            assert names == {"case.inspect", "checkpoint.inspect"}
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

            checkpoint = await client.call_tool(
                "checkpoint.inspect", {"checkpoint_path": "demo.dat.h5"}
            )
            assert checkpoint.is_error is False
            data = checkpoint.structured_content
            assert data["ok"] is True
            assert data["artifact"] == "checkpoint"
            assert data["attributes"]["iteration"] == 42
            assert data["attributes"]["time"] == 0.25
            assert data["cell_ids_present"] is True
            assert data["fields"][0]["name"] == "p"

            for name, args in (
                ("case.inspect", {"case_path": "../outside.cfdx.h5"}),
                ("checkpoint.inspect", {"checkpoint_path": "/etc/passwd"}),
            ):
                result = await client.call_tool(name, args)
                assert result.is_error is False
                assert result.structured_content["ok"] is False


if __name__ == "__main__":
    asyncio.run(exercise())
