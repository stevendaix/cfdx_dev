"""Read-only CFDX Runtime MCP server.

This first Runtime MCP surface inspects existing CFDX case/checkpoint artifacts.
It never creates, modifies, executes, or deletes cases or solver state.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from typing import Any

import h5py
from mcp.server import MCPServer
from mcp.types import ToolAnnotations


def _root(configured: str | None) -> Path:
    return Path(configured or os.environ.get("CFDX_RUNTIME_ROOT", ".")).resolve()


def _safe_path(root: Path, value: str) -> Path:
    if not value or Path(value).is_absolute():
        raise ValueError("path must be a non-empty root-relative path")
    candidate = (root / value).resolve()
    try:
        candidate.relative_to(root)
    except ValueError as exc:
        raise ValueError("path must remain inside the configured runtime root") from exc
    return candidate


def _scalar(value: Any) -> Any:
    if hasattr(value, "item"):
        value = value.item()
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _attributes(h5: h5py.File) -> dict[str, Any]:
    return {str(key): _scalar(value) for key, value in h5.attrs.items()}


def _datasets(group: h5py.Group, prefix: str = "") -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for name, item in group.items():
        path = f"{prefix}/{name}" if prefix else name
        if isinstance(item, h5py.Dataset):
            rows.append({"path": path, "shape": list(item.shape), "dtype": str(item.dtype), "size": int(item.size)})
        elif isinstance(item, h5py.Group):
            rows.extend(_datasets(item, path))
    return rows


def _inspect_case(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"ok": False, "errors": [f"case file does not exist: {path.name}"]}
    if not path.name.lower().endswith(".cfdx.h5"):
        return {"ok": False, "errors": ["case must use the canonical .cfdx.h5 extension"]}
    try:
        with h5py.File(path, "r") as h5:
            if _scalar(h5.attrs.get("format")) != "CFDX":
                return {"ok": False, "errors": ["not a CFDX case artifact"]}
            raw_config = h5["case/config"][()]
            if isinstance(raw_config, bytes):
                raw_config = raw_config.decode("utf-8")
            config = json.loads(raw_config)
            return {"ok": True, "artifact": "case", "path": path.name, "attributes": _attributes(h5), "configuration": config, "datasets": _datasets(h5), "runtime_state_present": "runtime" in h5}
    except (OSError, KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX case artifact: {exc}"]}


def _inspect_checkpoint(path: Path) -> dict[str, Any]:
    if not path.is_file():
        return {"ok": False, "errors": [f"checkpoint file does not exist: {path.name}"]}
    if not path.name.lower().endswith(".dat.h5"):
        return {"ok": False, "errors": ["checkpoint must use the canonical .dat.h5 extension"]}
    try:
        with h5py.File(path, "r") as h5:
            if _scalar(h5.attrs.get("format")) != "CFDX-DAT":
                return {"ok": False, "errors": ["not a CFDX DAT HDF5 checkpoint"]}
            if "fields" not in h5:
                return {"ok": False, "errors": ["checkpoint has no fields group"]}
            fields = [{"name": name, "shape": list(dataset.shape), "dtype": str(dataset.dtype), "size": int(dataset.size)} for name, dataset in h5["fields"].items()]
            return {"ok": True, "artifact": "checkpoint", "path": path.name, "attributes": _attributes(h5), "cell_ids_present": "cell_ids" in h5, "fields": fields}
    except (OSError, KeyError, TypeError, ValueError) as exc:
        return {"ok": False, "errors": [f"invalid CFDX DAT checkpoint: {exc}"]}


def create_server(root: str | None = None) -> MCPServer:
    runtime_root = _root(root)
    server = MCPServer("CFDX Runtime MCP", instructions="Read-only CFDX runtime artifact inspection. These tools inspect case/checkpoint files only; they never execute CFDX, mutate artifacts, or claim solver success.")
    annotations = ToolAnnotations(read_only_hint=True, open_world_hint=False)

    @server.tool(name="case.inspect", title="Inspect CFDX case", annotations=annotations)
    def case_inspect(case_path: str) -> dict[str, Any]:
        """Inspect a root-relative .cfdx.h5 case without loading solver state."""
        try:
            return _inspect_case(_safe_path(runtime_root, case_path))
        except ValueError as exc:
            return {"ok": False, "errors": [str(exc)]}

    @server.tool(name="checkpoint.inspect", title="Inspect CFDX checkpoint", annotations=annotations)
    def checkpoint_inspect(checkpoint_path: str) -> dict[str, Any]:
        """Inspect a root-relative .dat.h5 checkpoint without executing CFDX."""
        try:
            return _inspect_checkpoint(_safe_path(runtime_root, checkpoint_path))
        except ValueError as exc:
            return {"ok": False, "errors": [str(exc)]}

    return server


mcp = create_server()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=None)
    parser.add_argument("--transport", choices=("stdio", "streamable-http"), default="stdio")
    args = parser.parse_args()
    create_server(str(args.root) if args.root is not None else None).run(transport=args.transport)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
