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
        h5.create_dataset(
            "case/config",
            data='{"name":"demo","physics":{"laminar":{"enabled":true}},"numerics":{"scheme":"simple"},"boundaries":{},"materials":{},"execution":{"policy":"AUTO","solver":null,"mpi_ranks":1,"deterministic":true,"restart_option":"--restart"}}',
        )
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
        (root / "execution.json").write_text(
            """{
              "format": "CFDX-EXECUTION",
              "schema_version": 1,
              "process_exit_code": 0,
              "converged": true,
              "iterations": 42,
              "convergence_status": 1,
              "convergence_reason": "converged",
              "artifacts": {"restart_dat": true, "convergence_json": true}
            }""",
            encoding="utf-8",
        )
        client = Client(create_server(str(root)))
        async with client:
            listed = await client.list_tools()
            assert {tool.name for tool in listed.tools} == {
                "case.inspect",
                "case.validate",
                "checkpoint.inspect",
                "checkpoint.field.inspect",
                "checkpoint.compare",
                "convergence.inspect",
                "execution.inspect",
            }
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

            case_validation = await client.call_tool(
                "case.validate",
                {"case_path": "demo.cfdx.h5"},
            )
            assert case_validation.is_error is False
            assert case_validation.structured_content["ok"] is False
            assert (
                "mesh data is missing"
                in case_validation.structured_content["errors"][0]
            )

            checkpoint = await client.call_tool(
                "checkpoint.inspect",
                {"checkpoint_path": "demo.dat.h5"},
            )
            assert checkpoint.is_error is False
            data = checkpoint.structured_content
            assert data["ok"] is True
            assert data["artifact"] == "checkpoint"
            assert data["attributes"]["iteration"] == 42
            assert data["attributes"]["time"] == 0.25
            assert data["metadata"] == {
                "version": 2,
                "cells": 1,
                "iteration": 42,
                "time": 0.25,
            }
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
            invalid = await client.call_tool(
                "checkpoint.inspect", {"checkpoint_path": "invalid.dat.h5"}
            )
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

            convergence = root / "convergence.json"
            convergence.write_text(
                """{
                  "format": "CFDX-CONVERGENCE",
                  "schema_version": 1,
                  "converged": true,
                  "iterations": 1,
                  "reference_momentum_residual": 2.0,
                  "convergence_status": 1,
                  "convergence_reason": "converged",
                  "history": [{
                    "iteration": 1,
                    "momentum_residual": 1.0,
                    "pressure_residual": 0.5,
                    "continuity_l1": 0.1,
                    "continuity_linf": 0.2,
                    "continuity_normalized": 0.3,
                    "momentum_equation_residual": 0.4,
                    "momentum_equation_residual_relative": 0.5,
                    "velocity_change_inf": 0.01,
                    "pressure_change_inf": 0.02,
                    "effective_alpha_u": 0.7,
                    "effective_alpha_p": 0.8,
                    "nonlinear_convergence_metric": 0.5,
                    "momentum_linear_iterations": 11,
                    "pressure_linear_iterations": 12,
                    "pressure_correctors_used": 2,
                    "linear_tolerance_used": 1e-8,
                    "corrected_flux_continuity_linf": 0.01,
                    "reconstructed_velocity_continuity_linf": 0.02,
                    "flux_velocity_mismatch_linf": 0.03,
                    "momentum_equation_residual_components": [0.1, 0.2, 0.3],
                    "momentum_equation_residual_internal": 0.4,
                    "momentum_equation_residual_boundary": 0.5,
                    "pressure_gradient_linf": 0.6,
                    "pressure_gradient_l2": 0.7,
                    "momentum_residual_cell": 4,
                    "momentum_residual_no_pressure": 0.8,
                    "momentum_pressure_contribution": 0.9,
                    "momentum_residual_patch": "wall",
                    "mass_boundary_flux": 1.0,
                    "mass_global_cell_balance": 1.1,
                    "mass_local_l1": 1.2,
                    "mass_local_linf": 1.3,
                    "mass_local_l2": 1.4,
                    "mass_normalized_imbalance": 1.5,
                    "mass_worst_cell": 5,
                    "mass_nonfinite_faces": 0,
                    "boundedness_nonfinite_velocity": 0,
                    "momentum_conservation_residual": [0.2, 0.3, 0.4],
                    "momentum_conservation_normalized": [0.5, 0.6, 0.7],
                    "momentum_conservation_worst_cell": [1, 2, 3]
                  }]
                }""",
                encoding="utf-8",
            )
            execution_result = await client.call_tool(
                "execution.inspect",
                {"execution_path": "execution.json"},
            )
            assert execution_result.is_error is False
            data = execution_result.structured_content
            assert data["ok"] is True
            assert data["artifact"] == "execution"
            assert data["format"] == "CFDX-EXECUTION"
            assert data["schema_version"] == 1
            assert data["process_exit_code"] == 0
            assert data["converged"] is True
            assert data["iterations"] == 42
            assert data["artifacts"] == {"restart_dat": True, "convergence_json": True}

            bad_execution = root / "bad" / "execution.json"
            bad_execution.parent.mkdir()
            bad_execution.write_text(
                (root / "execution.json").read_text(encoding="utf-8").replace(
                    '"schema_version": 1', '"schema_version": 99'
                ),
                encoding="utf-8",
            )
            bad_execution_result = await client.call_tool(
                "execution.inspect",
                {"execution_path": "bad/execution.json"},
            )
            assert bad_execution_result.is_error is False
            assert bad_execution_result.structured_content["ok"] is False
            assert "schema version" in (
                bad_execution_result.structured_content["errors"][0]
            )

            nonfinite_execution = root / "nonfinite" / "execution.json"
            nonfinite_execution.parent.mkdir()
            nonfinite_execution.write_text(
                (root / "execution.json").read_text(encoding="utf-8").replace(
                    '"iterations": 42', '"iterations": NaN'
                ),
                encoding="utf-8",
            )
            bad_nonfinite_execution = await client.call_tool(
                "execution.inspect",
                {"execution_path": "nonfinite/execution.json"},
            )
            assert bad_nonfinite_execution.is_error is False
            assert bad_nonfinite_execution.structured_content["ok"] is False
            assert "non-finite constant" in (
                bad_nonfinite_execution.structured_content["errors"][0]
            )

            bad_exit = root / "bad-exit" / "execution.json"
            bad_exit.parent.mkdir()
            bad_exit.write_text(
                (root / "execution.json").read_text(encoding="utf-8").replace(
                    '"process_exit_code": 0', '"process_exit_code": 3'
                ),
                encoding="utf-8",
            )
            bad_exit_result = await client.call_tool(
                "execution.inspect",
                {"execution_path": "bad-exit/execution.json"},
            )
            assert bad_exit_result.is_error is False
            assert bad_exit_result.structured_content["ok"] is False
            assert (
                "0, 1, or 2"
                in bad_exit_result.structured_content["errors"][0]
            )

            inconsistent = root / "inconsistent" / "execution.json"
            inconsistent.parent.mkdir()
            inconsistent.write_text(
                (root / "execution.json").read_text(encoding="utf-8").replace(
                    '"converged": true', '"converged": false'
                ),
                encoding="utf-8",
            )
            inconsistent_result = await client.call_tool(
                "execution.inspect",
                {"execution_path": "inconsistent/execution.json"},
            )
            assert inconsistent_result.is_error is False
            assert inconsistent_result.structured_content["ok"] is False
            assert (
                "exit code 0 requires converged=true"
                in inconsistent_result.structured_content["errors"][0]
            )

            convergence_result = await client.call_tool(
                "convergence.inspect",
                {"convergence_path": "convergence.json", "history_limit": 1},
            )
            assert convergence_result.is_error is False
            data = convergence_result.structured_content
            assert data["ok"] is True
            assert data["format"] == "CFDX-CONVERGENCE"
            assert data["schema_version"] == 1
            assert data["converged"] is True
            assert data["history_count"] == 1
            assert len(data["history"]) == 1
            assert data["history_truncated"] is False

            malformed = root / "bad.convergence.json"
            malformed.write_text(
                convergence.read_text(encoding="utf-8").replace(
                    '"schema_version": 1', '"schema_version": 99'
                ),
                encoding="utf-8",
            )
            bad_schema = await client.call_tool(
                "convergence.inspect",
                {"convergence_path": "bad.convergence.json"},
            )
            assert bad_schema.is_error is False
            assert bad_schema.structured_content["ok"] is False
            assert "schema version" in bad_schema.structured_content["errors"][0]

            nonfinite = root / "nonfinite.convergence.json"
            nonfinite.write_text(
                convergence.read_text(encoding="utf-8").replace(
                    '"momentum_residual": 1.0', '"momentum_residual": NaN'
                ),
                encoding="utf-8",
            )
            bad_nonfinite = await client.call_tool(
                "convergence.inspect",
                {"convergence_path": "nonfinite.convergence.json"},
            )
            assert bad_nonfinite.is_error is False
            assert bad_nonfinite.structured_content["ok"] is False
            assert "non-finite JSON constant" in bad_nonfinite.structured_content["errors"][0]

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
                (
                    "checkpoint.field.inspect",
                    {"checkpoint_path": "../outside.dat.h5", "field_name": "p"},
                ),
                (
                    "checkpoint.compare",
                    {
                        "left_checkpoint_path": "../a.dat.h5",
                        "right_checkpoint_path": "later.dat.h5",
                        "field_name": "p",
                    },
                ),
                ("convergence.inspect", {"convergence_path": "../outside/convergence.json"}),
                ("convergence.inspect", {"convergence_path": "/etc/convergence.json"}),
                ("execution.inspect", {"execution_path": "../outside/execution.json"}),
                ("execution.inspect", {"execution_path": "/etc/execution.json"}),
                ("checkpoint.inspect", {"checkpoint_path": "/etc/passwd"}),
                ("checkpoint.field.inspect", {"checkpoint_path": "../outside.dat.h5", "field_name": "p"}),
            ):
                result = await client.call_tool(name, args)
                assert result.is_error is False
                assert result.structured_content["ok"] is False


if __name__ == "__main__":
    asyncio.run(exercise())
