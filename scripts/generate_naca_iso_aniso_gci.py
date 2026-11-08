"""Build matched NACA O-grid sequences with different *measured* wall anisotropy.

Each sequence is geometrically nested radially by subdivision of a common
24-layer base O-grid, while its surface nodes are independently refined.
The two alternatives use cosine clustering or approximately uniform arc-length
surface sampling of the same analytical NACA0012 contour. Both have identical
farfield, radial power, radial count and first-layer scale.

Neither naming convention nor this geometric gate proves CFD/grid convergence.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import math
from pathlib import Path

import numpy as np

SOURCE = Path(__file__).resolve().with_name("generate_2d_validation_meshes.py")


def load_generator():
    spec = importlib.util.spec_from_file_location("cfdx_naca_meshgen_gci", SOURCE)
    assert spec and spec.loader
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def arc_length_surface(n: int, native, cluster: float = 0.0):
    """Same analytic NACA contour, almost uniform arc-length panel spacing.

    Keep LE and both finite-thickness TE endpoints explicit; the reference
    cosine discretization is far denser than any exported wall grid.
    """
    dense_n = max(8192, 32 * n)
    samples = np.asarray(native(dense_n)[: dense_n + 1], dtype=np.float64)
    step = np.linalg.norm(np.diff(samples, axis=0), axis=1)
    arc = np.concatenate(([0.0], np.cumsum(step)))
    t = np.linspace(0.0, 1.0, n + 1)
    if cluster > 0:
        # Smooth LE clustering with a NONZERO endpoint derivative. Unlike
        # cosine clustering, the smallest wall segment refines as O(1/N).
        t = 1.0 - np.expm1(cluster * (1.0 - t)) / np.expm1(cluster)
    targets = float(arc[-1]) * t
    upper = np.column_stack(
        (
            np.interp(targets, arc, samples[:, 0]),
            np.interp(targets, arc, samples[:, 1]),
        )
    )
    # TE upper -> LE -> TE lower, with a single LE node.
    return [tuple(p) for p in upper] + [(float(x), -float(y)) for x, y in upper[-2::-1]]


def wall_aspect_metrics(points, n: int, n_radial: int) -> dict:
    count = 2 * n + 1
    xy = np.asarray(points, dtype=np.float64)[: 2 * count, :2]
    xy = xy.reshape(2, count, 2)
    chordwise = np.linalg.norm(np.roll(xy[0], -1, axis=0) - xy[0], axis=1)
    wallnormal = np.linalg.norm(xy[1] - xy[0], axis=1)
    if not np.all(np.isfinite(chordwise)) or np.min(chordwise) <= 0:
        raise ValueError("invalid NACA airfoil wall segment")
    ratio = np.maximum(chordwise / wallnormal, wallnormal / chordwise)
    return {
        "wall_edge_aspect_ratio_p50": float(np.percentile(ratio, 50)),
        "wall_edge_aspect_ratio_p95": float(np.percentile(ratio, 95)),
        "wall_edge_aspect_ratio_max": float(np.max(ratio)),
        "minimum_tangential_wall_edge": float(np.min(chordwise)),
        "maximum_tangential_wall_edge": float(np.max(chordwise)),
    }


def shear_interior_rings(mesh, n: int, n_radial: int, strength: float = 0.01):
    """Controlled volume-mesh skew at EXACTLY identical wall/farfield nodes.

    Tangential displacement is proportional to sin(pi * eta), so it vanishes
    on the physical boundaries. This is a geometry A/B perturbation; it does
    not solve CFD or change the PDE/physics. The first-cell spacing may still
    shift slightly and must be checked by the family-comparison gate.
    """
    nodes, cells, boundary = mesh
    m = 2 * n + 1
    points = np.asarray(nodes, dtype=np.float64)
    if len(points) % ((n_radial + 1) * m):
        raise ValueError("invalid O-grid node layout for shearing")
    span_levels = len(points) // ((n_radial + 1) * m)
    rings = points.reshape(span_levels, n_radial + 1, m, 3).copy()
    for j in range(1, n_radial):
        shift = strength * m * math.sin(math.pi * j / n_radial)
        base = math.floor(shift)
        fraction = shift - base
        for z in range(span_levels):
            ring = rings[z, j].copy()
            rings[z, j] = (1 - fraction) * np.roll(
                ring, -base, axis=0
            ) + fraction * np.roll(ring, -base - 1, axis=0)
    return rings.reshape(-1, 3).tolist(), cells, boundary


def contour_polygon_error(surface) -> dict:
    """Midpoint/interior vertical contour interpolation error, not force error.

    The finite-thickness TE closing edge is deliberately excluded: it is
    a separate physical edge, not one of the NACA upper/lower wall curves.
    """
    pts = np.asarray(surface, dtype=np.float64)

    def thickness(x):
        x = np.clip(x, 0.0, 1.0)
        return 0.6 * (
            0.2969 * np.sqrt(x)
            - 0.1260 * x
            - 0.3516 * x * x
            + 0.2843 * x**3
            - 0.1015 * x**4
        )

    errors = []
    for i in range(len(pts) - 1):
        one, two = pts[i], pts[i + 1]
        for frac in (0.25, 0.5, 0.75):
            approx = (1 - frac) * one + frac * two
            sign = 1.0 if one[1] + two[1] >= 0.0 else -1.0
            errors.append(abs(approx[1] - sign * thickness(approx[0])))
    return {
        "max_vertical_contour_interpolation_error": float(max(errors)),
        "rms_vertical_contour_interpolation_error": float(
            np.sqrt(np.mean(np.square(errors)))
        ),
        "note": (
            "vertical polygon-versus-analytic NACA error, "
            "not a normal-distance or force error"
        ),
    }


def build_families(
    out: Path,
    powers: float = 2.0,
    levels: tuple[int, ...] = (64, 128, 256),
    make_hdf5: bool = True,
    only_family: str | None = None,
    shear_strength: float = 0.01,
) -> dict:
    if powers < 1 or powers > 3:
        raise ValueError("invalid radial exponent")
    if not math.isfinite(shear_strength) or not 0.0 <= shear_strength <= 0.03:
        raise ValueError("interior shear strength must be between 0 and 0.03")
    gen = load_generator()
    native = gen.naca0012_surface
    out.mkdir(parents=True, exist_ok=True)
    result = {}
    original_grid = gen.naca_o_grid
    clusters = {
        "controlled_anisotropic": 3.2,
        "controlled_anisotropic_mild": 2.0,
    }
    parameterizations = {
        "quasi_isotropic": "arc_length",
        "sheared_isosurface": "arc_length",
        "controlled_anisotropic": "finite_slope_LE_clustering_3p2",
        "controlled_anisotropic_mild": "finite_slope_LE_clustering_2p0",
        "anisotropic": "cosine",
    }
    for family in (
        "quasi_isotropic",
        "sheared_isosurface",
        "controlled_anisotropic",
        "controlled_anisotropic_mild",
        "anisotropic",
    ):
        if only_family is not None and family != only_family:
            continue
        gen.naca_o_grid = original_grid
        if family == "anisotropic":
            gen.naca0012_surface = native
        else:
            cluster = clusters.get(family, 0.0)
            gen.naca0012_surface = lambda n, native=native, cluster=cluster: (
                arc_length_surface(n, native, cluster=cluster)
            )
        if family == "sheared_isosurface":
            gen.naca_o_grid = lambda n, nr, *a, **kw: shear_interior_rings(
                original_grid(n, nr, *a, **kw), n, nr, strength=shear_strength
            )
        rows = []
        (out / family).mkdir(exist_ok=True)
        for n in levels:
            if n % 64 != 0 or (24 * n) % 64 != 0:
                raise ValueError("only N64/N128/N256 nested levels supported")
            nrad = 24 * n // 64
            points, cells, boundaries = gen.naca_o_grid(
                n, nrad, 20.0, radial_power=powers, base_radial_layers=24
            )
            row = gen.naca_o_grid_quality(
                n, nrad, 20.0, radial_power=powers, base_radial_layers=24
            )
            row.update(wall_aspect_metrics(points, n, nrad))
            row.update(contour_polygon_error(gen.naca0012_surface(n)))
            row["interior_shear_strength"] = (
                shear_strength if family == "sheared_isosurface" else 0.0
            )
            row["surface_parameterization"] = parameterizations[family]
            if row["inverted_or_degenerate_cells"]:
                raise ValueError(f"{family} N{n} has inverted cells")
            rows.append(row)
            if make_hdf5:
                gen.write_mesh(
                    out / family / f"naca0012_{family}_n{n}.h5",
                    points,
                    cells,
                    boundaries,
                )
        (out / family / "quality.json").write_text(
            json.dumps(rows, indent=2) + "\n", encoding="utf-8"
        )
        result[family] = rows
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--radial-power", type=float, default=2.0)
    p.add_argument("--shear-strength", type=float, default=0.01)
    p.add_argument("--quality-only", action="store_true")
    p.add_argument(
        "--only-family",
        choices=(
            "quasi_isotropic",
            "sheared_isosurface",
            "controlled_anisotropic",
            "controlled_anisotropic_mild",
            "anisotropic",
        ),
    )
    args = p.parse_args()
    result = build_families(
        args.output_dir,
        args.radial_power,
        make_hdf5=not args.quality_only,
        only_family=args.only_family,
        shear_strength=args.shear_strength,
    )
    for family, rows in result.items():
        for r in rows:
            print(
                f"{family} N{r['n_surface']} "
                f"cells={r['cell_count']} "
                f"wall_AR_median={r['wall_edge_aspect_ratio_p50']:.2f} "
                f"wall_AR_max={r['wall_edge_aspect_ratio_max']:.2f} "
                f"nonorthog_max={r['maximum_nonorthogonality_degrees']:.1f}"
            )


if __name__ == "__main__":
    main()
