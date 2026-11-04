"""Independent, fail-closed TH-03 input/geometry contract audit.

This does NOT validate a CFD solution. It distinguishes the geometrical
OpenFOAM-style surrogate from the frozen Driver-Seegmiller experimental case,
checks declared SI units, inlet turbulence, boundary semantics and refinement.
No tolerance is modified and no simulation is declared converged here.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import h5py
import numpy as np

# Reuse the independent geometric/topological checker, not the case generator's
# metadata-only counts. This rejects inverted faces and false EMPTY patches.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from audit_th03_openfoam_mesh import audit as audit_mesh, read_hdf5

REFERENCE = Path(__file__).resolve().parents[1] / (
    'docs/validation/14-benchmarks/09-thermohydraulics/TH-03.json')
NOMINAL_INLET_VISCOSITY_RATIO = 0.009


def _bc(setup: dict, patch: str) -> dict:
    matches = [x for x in setup['boundary_conditions'] if x['patch_name'] == patch]
    if len(matches) != 1:
        raise ValueError(f'TH-03: expected exactly one {patch} boundary')
    return matches[0]


def _scalar(bc: dict, field: str) -> float:
    matches = [x for x in bc.get('sst_scalar_conditions', []) if x['field'] == field]
    if len(matches) != 1 or matches[0]['type'] != 'fixed_value':
        raise ValueError(f'TH-03: requires explicit inlet {field} fixed_value')
    v = float(matches[0]['value'])
    if not math.isfinite(v) or v <= 0:
        raise ValueError(f'TH-03: invalid inlet {field}')
    return v


def _close(a: float, b: float, relative: float = 1e-8) -> bool:
    return math.isfinite(a) and math.isclose(a, b, rel_tol=relative, abs_tol=1e-12)


def audit(path: Path) -> dict:
    mesh_quality = audit_mesh(read_hdf5(path))
    reference = json.loads(REFERENCE.read_text(encoding='utf-8'))
    phys = reference['physics']
    h = float(phys['geometry']['step_height_m'])
    flow = phys['flow']
    with h5py.File(path, 'r') as f:
        points = np.asarray(f['points'], dtype=float)
        cell_count = len(f['cell_offsets']) - 1
        if points.ndim != 2 or points.shape[1] != 3 or not np.isfinite(points).all():
            raise ValueError('TH-03: invalid finite 3D mesh coordinates')
        setup = json.loads(f.attrs['case_setup_json'])
        count = int(f.attrs['n_cells'])
    if count != cell_count:
        raise ValueError('TH-03: cell count metadata mismatch')
    if setup['units'] != 'SI' or setup['turbulence_model'] != 'k_omega_sst':
        raise ValueError('TH-03: expected SI, k-omega SST case')
    if len(setup['materials']) != 1:
        raise ValueError('TH-03: expected one homogeneous fluid material')
    material = setup['materials'][0]
    rho = float(material['density'])
    mu = float(material['dynamic_viscosity'])
    if not (rho > 0 and mu > 0 and math.isfinite(rho) and math.isfinite(mu)):
        raise ValueError('TH-03: invalid material density / dynamic viscosity')
    nu = mu / rho
    inlet = _bc(setup, 'inlet')
    if inlet['type'] != 'inlet' or inlet['value_type'] != 'fixed':
        raise ValueError('TH-03: requires velocity Dirichlet inlet')
    ux = float(inlet['velocity_vector'][0])
    k = _scalar(inlet, 'k')
    omega = _scalar(inlet, 'omega')
    if ux <= 0 or not math.isfinite(ux):
        raise ValueError('TH-03: invalid reference inlet streamwise velocity')
    reynolds = ux * h / nu
    ratio = k / (nu * omega)
    intensity = math.sqrt(2.0 * k / 3.0) / ux
    omega_initial = float(setup['initial_condition']['scalar_fields']['omega'])
    raw = setup['numerics']['raw_settings']
    pressure_tol = float(raw.get('pressure_linear_tolerance', 'nan'))
    if not (math.isfinite(pressure_tol) and pressure_tol > 0):
        raise ValueError('TH-03: pressure tolerance must be finite, positive')
    pmin = points.min(axis=0) / h
    pmax = points.max(axis=0) / h
    outlet = points[np.isclose(points[:, 0] / h, 50, atol=1e-9)]
    step = points[np.isclose(points[:, 0] / h, 0, atol=1e-9)]
    if not len(outlet) or not len(step):
        raise ValueError('TH-03: missing outlet or step mesh vertices')
    floor_exit = float(outlet[:, 1].min() / h)
    step_height_over_h = float((step[:, 1].min()) / h)
    geobounds = (
        np.allclose(pmin, [-130, 0, 0], atol=1e-6) and
        np.allclose(pmax, [50, 9, 1], atol=1e-6)
    )
    # Estimate ONLY using the upstream measured Cf; true y+ requires the
    # converged local CFDX wall traction, especially around reattachment.
    y_at_step = np.sort(np.unique(np.round(step[:, 1] / h, 12)))
    positive = y_at_step[y_at_step > 1e-10]
    near_step_y = float(positive.min() * h) if len(positive) else float('nan')
    cf_reference = float(reference['physics']['inlet_boundary_layer_at_x_over_h_minus4']['skin_friction_coefficient'])
    friction_velocity = ux * math.sqrt(cf_reference / 2.)
    rough_yplus = 0.5 * near_step_y * friction_velocity / nu
    wall = _bc(setup, 'lowerWall')
    wall_kinds = {x['field']: x['type'] for x in wall['sst_scalar_conditions']}
    wall_resolved = (wall['value_type'] == 'wall_no_slip' and
                     wall_kinds == {'k': 'wall_k_zero', 'omega': 'wall_omega_sst'})
    current_convection = [x['configuration_key'] for x in
                          setup['numerics']['selection']['entries'] if x['family'] == 'convection']
    first_order = 'numerics.convection.upwind' in current_convection
    experimental_re = float(flow['reynolds_step'])
    # The frozen CMP32 experiment is a straight recovery duct with measured
    # x/h=-4 approach profiles. The reconstructed tutorial has straight
    # EXTERNAL walls: its internal block interface rises to y/h=3.125 at exit.
    # Do not confuse the internal partition with a sloping physical wall.
    checks = {
        'six_block_domain_bounding_box': bool(geobounds),
        'inlet_nut_over_nu_matches_documented_openfoam_0p009': _close(ratio, NOMINAL_INLET_VISCOSITY_RATIO, 1e-3),
        'inlet_turbulence_intensity_matches_nominal_openfoam_0p061percent': _close(intensity, 0.00061, 0.01),
        'experimental_Re_h_37500': _close(reynolds, experimental_re, 1e-3),
        'straight_recovery_floor_at_outlet': abs(floor_exit) < 1e-8,
        'measured_x_over_h_minus4_profile_declared': False,
        'wall_model_matches_openfoam_wall_functions': False,
        'second_order_convection_selected': not first_order,
    }
    return {
        'case_file': str(path),
        'cells': int(cell_count),
        'dimensions': {'bounds_min_over_h': pmin.tolist(), 'bounds_max_over_h': pmax.tolist(),
                       'outlet_lower_wall_y_over_h': floor_exit,
                       'wall_height_at_step_over_h': step_height_over_h},
        'physics': {'step_height_m': h, 'inlet_velocity_m_s': ux,
                    'density_kg_m3': rho, 'kinematic_viscosity_m2_s': nu,
                    'reynolds_h': reynolds, 'experimental_reynolds_h': experimental_re,
                    'inlet_k_m2_s2': k, 'inlet_omega_s_inverse': omega,
                    'initial_omega_s_inverse': omega_initial,
                    'inlet_nut_over_nu': ratio,
                    'inlet_turbulence_intensity': intensity,
                    'nominal_openfoam_inlet_nut_over_nu': NOMINAL_INLET_VISCOSITY_RATIO,
                    'estimated_first_cell_yplus_from_upstream_experimental_cf': rough_yplus,
                    'yplus_note': 'diagnostic estimate, NOT converged-CFD wall y+',
                    'wall_resolved_sst': wall_resolved,
                    'pressure_linear_tolerance': pressure_tol},
        'checks': checks,
        'mesh_status': ('SURROGATE_GEOMETRY_ONLY' if geobounds and
                        mesh_quality['status'] == 'PASS_GEOMETRY_ONLY' else 'FAIL_GEOMETRY'),
        'mesh_quality': mesh_quality,
        'experimental_comparability': 'NOT_COMPARABLE',
        'physical_qualification': 'NOT_QUALIFIED',
        'warnings': [
            'OpenFOAM-style reconstructed mesh is not native blockMesh polyMesh.',
            'NASA nominal Re_h=36000 is not the frozen experimental Re_h=37500.',
            'Experimental measured approach profiles not supplied; uniform inlet is not equivalent.',
            'Exterior walls are straight, but exact OpenFOAM blockMesh grading is not verified.',
            'Wall-resolved SST and OpenFOAM wall functions are different models.',
            'No convergence field, mass budget, profile, Cf, Cp or three-grid GCI is certified here.',
        ],
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('case', type=Path)
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    result = audit(args.case)
    serialized = json.dumps(result, indent=2, allow_nan=False) + '\n'
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(serialized, encoding='utf-8')
    print(serialized)
    # A preflight is useful even when it cannot qualify the physics.
    # 2 means geometry is broken, not merely 'not yet qualified'.
    return 2 if result['mesh_status'] == 'FAIL_GEOMETRY' else 0


if __name__ == '__main__':
    raise SystemExit(main())
