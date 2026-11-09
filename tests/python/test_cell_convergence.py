"""Fail-closed hotspot trace contracts, non-finite/mismatched/scoped reads."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from cfdx.cell_convergence import (
    CellConvergenceError,
    cell_trace_metric_name,
    read_cell_trace,
    sampled_cell_scalar_map,
)


def fixture_trace(tmp_path: Path) -> Path:
    source = tmp_path / "cell_convergence.jsonl"
    header = {
        "kind": "header",
        "format": "CFDX-CELL-CONVERGENCE",
        "schema_version": 1,
        "mesh_source": "/work/mesh.cfdx.h5",
        "n_cells": 6,
        "n_faces": 20,
        "every": 1,
        "top_k": 3,
        "reference_relative_tolerance": 1e-6,
    }
    cell = {
        "cell_id": 2,
        "center": [-1.0, 0.1, 6],
        "momentum_abs": [0.2, 0, 0.02],
        "momentum_scaled": 0.1,
        "mass_imbalance_abs": 0.004,
        "nondecrease_streak": 3,
    }
    frame = {
        "kind": "frame",
        "iteration": 4,
        "assessed": 6,
        "nonfinite": 0,
        "momentum_above_reference": 5,
        "persistent_momentum_cells": 2,
        "operator_comparable": True,
        "rhs_scale": [2, 1, 1],
        "cells": [cell],
    }
    end = {"kind": "end", "iterations": 4, "last_frame_iteration": 4}
    source.write_text("".join(json.dumps(row) + "\n" for row in (header, frame, end)))
    return source


def test_correct_spatial_hotspot_and_no_qualification_claim(tmp_path):
    data = read_cell_trace(fixture_trace(tmp_path))
    assert data.complete
    assert data.accepted_iterations == 4
    assert len(data.frames) == 1
    sample = data.latest.cells[0]
    assert sample.center == (-1.0, 0.1, 6.0)
    assert sample.stagnant
    assert sample.cell_id == 2
    assert sample.momentum_scaled == pytest.approx(0.1)
    assert data.latest.ranked("persistence") == (sample,)
    assert not data.to_dict()["solver_convergence_certificate"]
    assert "top-K" in data.to_dict()["spatial_coverage"]
    assert len(data.sha256) == 64


def test_live_file_only_discards_uncommitted_final_fragment(tmp_path):
    source = fixture_trace(tmp_path)
    payload = source.read_bytes().splitlines(keepends=True)
    source.write_bytes(payload[0] + payload[1] + b'{"kind":"end","iteration')
    trace = read_cell_trace(source)
    assert not trace.complete
    assert trace.latest.iteration == 4


def test_header_only_live_file(tmp_path):
    source = fixture_trace(tmp_path)
    source.write_bytes(source.read_bytes().splitlines(keepends=True)[0])
    assert read_cell_trace(source).latest is None


@pytest.mark.parametrize(
    "mutate,match",
    [
        (lambda h, f, e: h.update({"format": "WRONG"}), "format|supported"),
        (lambda h, f, e: h.update({"mesh_source": ""}), "mesh source"),
        (lambda h, f, e: h.update({"n_cells": True}), "n_cells"),
        (lambda h, f, e: h.update({"top_k": 0}), "top_k"),
        (lambda h, f, e: f.update({"assessed": 7}), "counts"),
        (lambda h, f, e: f.update({"nonfinite": 1}), "counts"),
        (lambda h, f, e: f.update({"operator_comparable": "true"}), "boolean"),
        (lambda h, f, e: f.update({"rhs_scale": [1, 0, 1]}), "positive"),
        (
            lambda h, f, e: f.update({"cells": [{**f["cells"][0], "cell_id": 6}]}),
            "out-of-range",
        ),
        (
            lambda h, f, e: f.update({"cells": [f["cells"][0], f["cells"][0]]}),
            "duplicate",
        ),
        (
            lambda h, f, e: f["cells"][0].update({"center": [0, float("nan"), 0]}),
            "non-finite",
        ),
        (
            lambda h, f, e: f["cells"][0].update({"momentum_scaled": 0.5}),
            "does not match",
        ),
        (
            lambda h, f, e: f["cells"][0].update({"nondecrease_streak": 7}),
            "persistence",
        ),
        (
            lambda h, f, e: f["cells"][0].update({"mass_imbalance_abs": -1}),
            "nonnegative",
        ),
        (lambda h, f, e: e.update({"iterations": 3}), "disagrees"),
        (lambda h, f, e: e.update({"last_frame_iteration": 3}), "disagrees"),
    ],
)
def test_invalid_trace_fail_closed(tmp_path, mutate, match):
    source = fixture_trace(tmp_path)
    h, f, e = [json.loads(line) for line in source.read_text().splitlines()]
    mutate(h, f, e)
    source.write_text("\n".join(json.dumps(row) for row in (h, f, e)) + "\n")
    with pytest.raises(CellConvergenceError, match=match):
        read_cell_trace(source)


def test_duplicate_key_rejected(tmp_path):
    source = fixture_trace(tmp_path)
    source.write_text(
        source.read_text().replace('"cell_id": 2', '"cell_id": 2, "cell_id": 2')
    )
    with pytest.raises(CellConvergenceError, match="duplicate"):
        read_cell_trace(source)


def test_file_budgets_and_symlinks_fail_closed(tmp_path):
    source = fixture_trace(tmp_path)
    source.write_bytes(source.read_bytes() + b" " * 2048)
    with pytest.raises(CellConvergenceError, match="budget"):
        read_cell_trace(source, max_bytes=1024)
    alias = tmp_path / "alias.jsonl"
    alias.symlink_to(source)
    with pytest.raises(CellConvergenceError, match="open trace safely"):
        read_cell_trace(alias)


def test_incomplete_or_out_of_order_frames_not_certified(tmp_path):
    source = fixture_trace(tmp_path)
    lines = source.read_text().splitlines()
    source.write_text(
        lines[0] + "\n" + lines[1] + "\n" + lines[1] + "\n" + lines[2] + "\n"
    )
    with pytest.raises(CellConvergenceError, match="increase"):
        read_cell_trace(source)


def test_no_corruption_if_last_partial_line_is_not_ascii(tmp_path):
    source = fixture_trace(tmp_path)
    lines = source.read_bytes().splitlines(keepends=True)
    source.write_bytes(b"".join(lines[:2]) + b"\xff\xfe")
    assert not read_cell_trace(source).complete


def test_jsonl_cli_smoke(tmp_path, capsys):
    from cfdx.cell_convergence import main

    assert main([str(fixture_trace(tmp_path)), "--metric", "persistence"]) == 0
    text = capsys.readouterr().out
    assert "cell        2" in text
    assert "no convergence certificate" in text


def test_ranked_limits(tmp_path):
    frame = read_cell_trace(fixture_trace(tmp_path)).latest
    with pytest.raises(ValueError, match="momentum"):
        frame.ranked("other")
    with pytest.raises(ValueError, match="limit"):
        frame.ranked(limit=0)


def test_no_accepted_iteration_complete_but_not_converged(tmp_path):
    source = fixture_trace(tmp_path)
    header = source.read_text().splitlines()[0]
    source.write_text(
        header
        + "\n"
        + json.dumps({"kind": "end", "iterations": 0, "last_frame_iteration": 0})
        + "\n"
    )
    trace = read_cell_trace(source)
    assert trace.complete
    assert trace.accepted_iterations == 0
    assert trace.latest is None
    assert not trace.to_dict()["solver_convergence_certificate"]


def test_end_marker_allows_last_periodic_frame_before_last_iteration(tmp_path):
    source = fixture_trace(tmp_path)
    h, f, e = [json.loads(line) for line in source.read_text().splitlines()]
    e["iterations"] = 7
    source.write_text("\n".join(json.dumps(row) for row in (h, f, e)) + "\n")
    trace = read_cell_trace(source)
    assert trace.complete
    assert trace.accepted_iterations == 7
    assert trace.latest.iteration == 4


def test_reject_end_with_zero_accepted_and_positive_frame(tmp_path):
    source = fixture_trace(tmp_path)
    lines = source.read_text().splitlines()
    e = json.loads(lines[-1])
    e["iterations"] = 0
    source.write_text("\n".join(lines[:-1] + [json.dumps(e)]) + "\n")
    with pytest.raises(CellConvergenceError, match="end marker"):
        read_cell_trace(source)


def test_observed_history_does_not_interpolate_missing_samples(tmp_path):
    trace = read_cell_trace(fixture_trace(tmp_path))
    history = trace.observed_history(2)
    assert len(history) == 1
    assert history[0][0] == 4
    assert history[0][1].cell_id == 2
    assert trace.observed_history(3) == ()  # unobserved is NOT converged
    with pytest.raises(ValueError, match="cell ID"):
        trace.observed_history(trace.n_cells)
    with pytest.raises(ValueError, match="cell ID"):
        trace.observed_history(True)


def test_sampled_sparse_cell_map_preserves_unknown_cells(tmp_path):
    trace = read_cell_trace(fixture_trace(tmp_path))
    title, values = sampled_cell_scalar_map(trace, trace.latest, "momentum")
    assert title == "Momentum / RHS"
    assert len(values) == trace.n_cells
    assert values[2] == pytest.approx(0.1)
    assert str(values[0]) == "nan"
    title, values = sampled_cell_scalar_map(trace, trace.latest, "persistence")
    assert title == "Non-decrease streak"
    assert values[2] == 3.0
    with pytest.raises(ValueError, match="metric"):
        cell_trace_metric_name("wrong")
