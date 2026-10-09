from __future__ import annotations

import asyncio
import os
import pathlib
import sys
import tempfile
from unittest.mock import patch

import h5py
from mcp import Client

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
sys.path.insert(
    0, str(pathlib.Path(__file__).resolve().parents[2] / "src" / "cfdx" / "python")
)

from cfdx.case import Case
from cfdx.case_io import save_case, validate_case_bundle
from cfdx.session import CFDXSession
from runtime_server import create_server


def make_valid_case(path: pathlib.Path) -> pathlib.Path:
    session = CFDXSession(
        case=Case(
            name="valid-case",
            physics={"laminar": {"enabled": True}},
            numerics={"scheme": "simple"},
        )
    )
    save_case(session, path)
    with h5py.File(path, "a") as h5:
        h5.create_dataset("points", data=[[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        h5.create_dataset("face_vertices", data=[0, 1])
        h5.create_dataset("face_offsets", data=[0, 2])
        h5.create_dataset("owner", data=[0])
        h5.create_dataset("neighbour", data=[])
        h5.create_dataset("cell_faces", data=[0])
        h5.create_dataset("cell_offsets", data=[0, 1])
        h5.require_group("fields").create_dataset("values", data=[0.0])
    assert all(validate_case_bundle(path).values())
    return path


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
        valid_case_path = make_valid_case(root / "valid.cfdx.h5")
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
                "execution.run",
                "case.configure",
            }
            for tool in listed.tools:
                annotations = tool.model_dump(by_alias=True).get("annotations", {})
                assert annotations["openWorldHint"] is False
                assert annotations["readOnlyHint"] is (
                    tool.name not in {"execution.run", "case.configure"}
                )

            disabled_execution = await client.call_tool(
                "execution.run", {"case_path": "demo.cfdx.h5"}
            )
            assert disabled_execution.is_error is False
            assert disabled_execution.structured_content["ok"] is False
            assert disabled_execution.structured_content["executed"] is False
            assert (
                "controlled execution is disabled"
                in (disabled_execution.structured_content["errors"][0])
            )

            invalid_timeout = await client.call_tool(
                "execution.run", {"case_path": "demo.cfdx.h5", "timeout": 0}
            )
            assert invalid_timeout.is_error is False
            assert invalid_timeout.structured_content["ok"] is False
            assert (
                "timeout must be greater than 0"
                in (invalid_timeout.structured_content["errors"][0])
            )

            old_write = os.environ.pop("CFDX_RUNTIME_ALLOW_WRITE", None)
            try:
                denied_write = await client.call_tool(
                    "case.configure",
                    {
                        "case_path": "valid.cfdx.h5",
                        "updates": {"numerics": {"scheme": "denied"}},
                    },
                )
                assert denied_write.is_error is False
                assert denied_write.structured_content["ok"] is False
                assert denied_write.structured_content["written"] is False
                assert (
                    "case configuration is disabled"
                    in (denied_write.structured_content["errors"][0])
                )
                os.environ["CFDX_RUNTIME_ALLOW_WRITE"] = "1"

                traversal_write = await client.call_tool(
                    "case.configure",
                    {
                        "case_path": "../outside.cfdx.h5",
                        "updates": {"name": "escape"},
                    },
                )
                assert traversal_write.structured_content["ok"] is False

                outside_case = root.parent / (root.name + "-outside.cfdx.h5")
                outside_case.write_bytes(valid_case_path.read_bytes())
                escape_link = root / "escape.cfdx.h5"
                escape_link.symlink_to(outside_case)
                symlink_write = await client.call_tool(
                    "case.configure",
                    {
                        "case_path": "escape.cfdx.h5",
                        "updates": {"name": "escape"},
                    },
                )
                assert symlink_write.structured_content["ok"] is False
                assert (
                    "inside the configured runtime root"
                    in (symlink_write.structured_content["errors"][0])
                )
                outside_case.unlink(missing_ok=True)

                invalid_write = await client.call_tool(
                    "case.configure",
                    {
                        "case_path": "valid.cfdx.h5",
                        "updates": {"numerics": ["not", "an", "object"]},
                    },
                )
                assert invalid_write.structured_content["ok"] is False

                configured = await client.call_tool(
                    "case.configure",
                    {
                        "case_path": "valid.cfdx.h5",
                        "updates": {
                            "name": "updated-case",
                            "numerics": {"scheme": "coupled"},
                        },
                    },
                )
                assert configured.is_error is False
                configured_data = configured.structured_content
                assert configured_data["ok"] is True
                assert configured_data["written"] is True
                assert configured_data["updated_sections"] == ["name", "numerics"]
                assert configured_data["case_revision"] == 1
                assert configured_data["numerics_revision"] == 1
                with h5py.File(valid_case_path, "r") as h5:
                    config = h5["case/config"][()]
                    if isinstance(config, bytes):
                        config = config.decode("utf-8")
                    assert '"updated-case"' in config
                    assert "points" in h5
                    assert "face_vertices" in h5
                    assert "fields/values" in h5
                assert all(validate_case_bundle(valid_case_path).values())

                before_failed_write = valid_case_path.read_bytes()
                with patch(
                    "runtime_server.os.replace",
                    side_effect=OSError("injected atomic replace failure"),
                ):
                    failed_write = await client.call_tool(
                        "case.configure",
                        {
                            "case_path": "valid.cfdx.h5",
                            "updates": {"name": "must-not-commit"},
                        },
                    )
                assert failed_write.structured_content["ok"] is False
                assert failed_write.structured_content["written"] is False
                assert valid_case_path.read_bytes() == before_failed_write
                assert not list(root.glob(".valid.*.cfdx.h5"))
            finally:
                if old_write is None:
                    os.environ.pop("CFDX_RUNTIME_ALLOW_WRITE", None)
                else:
                    os.environ["CFDX_RUNTIME_ALLOW_WRITE"] = old_write

            solver = root / "fake-solver"
            solver.write_text(
                "#!/usr/bin/env python3\n"
                "print('O' * 20000)\n"
                "import sys\n"
                "print('E' * 20000, file=sys.stderr)\n",
                encoding="utf-8",
            )
            solver.chmod(0o755)
            old_enable = os.environ.get("CFDX_RUNTIME_ALLOW_EXECUTE")
            old_solver = os.environ.get("CFDX_RUNTIME_SOLVER")
            os.environ["CFDX_RUNTIME_ALLOW_EXECUTE"] = "1"
            os.environ["CFDX_RUNTIME_SOLVER"] = str(solver)
            try:
                successful_execution = await client.call_tool(
                    "execution.run",
                    {"case_path": "demo.cfdx.h5", "timeout": 10},
                )
                assert successful_execution.is_error is False
                execution_data = successful_execution.structured_content
                assert execution_data["ok"] is True
                assert execution_data["executed"] is True
                assert execution_data["timed_out"] is False
                assert execution_data["returncode"] == 0
                assert len(execution_data["stdout"]) <= 16384
                assert len(execution_data["stderr"]) <= 16384
                assert execution_data["stdout"].rstrip().endswith("O" * 100)
                assert execution_data["stderr"].rstrip().endswith("E" * 100)

                outside_case = root.parent / (root.name + "-execution-outside.cfdx.h5")
                outside_case.write_bytes((root / "demo.cfdx.h5").read_bytes())
                execution_escape = root / "escape-execution.cfdx.h5"
                execution_escape.symlink_to(outside_case)
                solver_marker = root / "solver-was-launched.marker"
                solver.write_text(
                    "#!/usr/bin/env python3\\n"
                    "from pathlib import Path\\n"
                    f"Path({str(solver_marker)!r}).write_text('launched')\\n",
                    encoding="utf-8",
                )
                solver.chmod(0o755)
                rejected_execution = await client.call_tool(
                    "execution.run",
                    {"case_path": execution_escape.name, "timeout": 10},
                )
                rejected_data = rejected_execution.structured_content
                assert rejected_data["ok"] is False
                assert rejected_data["executed"] is False
                assert "inside the configured runtime root" in rejected_data["errors"][0]
                assert not solver_marker.exists()
                execution_escape.unlink()
                outside_case.unlink()

                solver.write_text(
                    "#!/usr/bin/env python3\n"
                    "import signal, time\n"
                    "signal.signal(signal.SIGTERM, signal.SIG_IGN)\n"
                    "while True: time.sleep(0.1)\n",
                    encoding="utf-8",
                )
                solver.chmod(0o755)
                timed_execution = await client.call_tool(
                    "execution.run",
                    {"case_path": "demo.cfdx.h5", "timeout": 0.1},
                )
                assert timed_execution.is_error is False
                timeout_data = timed_execution.structured_content
                assert timeout_data["ok"] is False
                assert timeout_data["executed"] is True
                assert timeout_data["timed_out"] is True
                assert timeout_data["returncode"] is not None
                assert "exceeded timeout" in timeout_data["errors"][0]
            finally:
                if old_enable is None:
                    os.environ.pop("CFDX_RUNTIME_ALLOW_EXECUTE", None)
                else:
                    os.environ["CFDX_RUNTIME_ALLOW_EXECUTE"] = old_enable
                if old_solver is None:
                    os.environ.pop("CFDX_RUNTIME_SOLVER", None)
                else:
                    os.environ["CFDX_RUNTIME_SOLVER"] = old_solver

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
                (root / "execution.json")
                .read_text(encoding="utf-8")
                .replace('"schema_version": 1', '"schema_version": 99'),
                encoding="utf-8",
            )
            bad_execution_result = await client.call_tool(
                "execution.inspect",
                {"execution_path": "bad/execution.json"},
            )
            assert bad_execution_result.is_error is False
            assert bad_execution_result.structured_content["ok"] is False
            assert (
                "schema version"
                in (bad_execution_result.structured_content["errors"][0])
            )

            nonfinite_execution = root / "nonfinite" / "execution.json"
            nonfinite_execution.parent.mkdir()
            nonfinite_execution.write_text(
                (root / "execution.json")
                .read_text(encoding="utf-8")
                .replace('"iterations": 42', '"iterations": NaN'),
                encoding="utf-8",
            )
            bad_nonfinite_execution = await client.call_tool(
                "execution.inspect",
                {"execution_path": "nonfinite/execution.json"},
            )
            assert bad_nonfinite_execution.is_error is False
            assert bad_nonfinite_execution.structured_content["ok"] is False
            assert (
                "non-finite constant"
                in (bad_nonfinite_execution.structured_content["errors"][0])
            )

            bad_exit = root / "bad-exit" / "execution.json"
            bad_exit.parent.mkdir()
            bad_exit.write_text(
                (root / "execution.json")
                .read_text(encoding="utf-8")
                .replace('"process_exit_code": 0', '"process_exit_code": 3'),
                encoding="utf-8",
            )
            bad_exit_result = await client.call_tool(
                "execution.inspect",
                {"execution_path": "bad-exit/execution.json"},
            )
            assert bad_exit_result.is_error is False
            assert bad_exit_result.structured_content["ok"] is False
            assert "0, 1, or 2" in bad_exit_result.structured_content["errors"][0]

            inconsistent = root / "inconsistent" / "execution.json"
            inconsistent.parent.mkdir()
            inconsistent.write_text(
                (root / "execution.json")
                .read_text(encoding="utf-8")
                .replace('"converged": true', '"converged": false'),
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
            assert (
                "non-finite JSON constant"
                in bad_nonfinite.structured_content["errors"][0]
            )

            with h5py.File(root / "later.dat.h5", "w") as h5:
                h5.attrs["format"] = "CFDX-DAT"
                h5.create_group("fields").create_dataset("p", data=[2.0])
            comparison = await client.call_tool(
                "checkpoint.compare",
                {
                    "left_checkpoint_path": "demo.dat.h5",
                    "right_checkpoint_path": "later.dat.h5",
                    "field_name": "p",
                },
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
                (
                    "convergence.inspect",
                    {"convergence_path": "../outside/convergence.json"},
                ),
                ("convergence.inspect", {"convergence_path": "/etc/convergence.json"}),
                ("execution.inspect", {"execution_path": "../outside/execution.json"}),
                ("execution.inspect", {"execution_path": "/etc/execution.json"}),
                ("execution.run", {"case_path": "../outside.cfdx.h5"}),
                (
                    "case.configure",
                    {"case_path": "../outside.cfdx.h5", "updates": {"name": "bad"}},
                ),
                ("checkpoint.inspect", {"checkpoint_path": "/etc/passwd"}),
                (
                    "checkpoint.field.inspect",
                    {"checkpoint_path": "../outside.dat.h5", "field_name": "p"},
                ),
            ):
                result = await client.call_tool(name, args)
                assert result.is_error is False
                assert result.structured_content["ok"] is False


if __name__ == "__main__":
    asyncio.run(exercise())
