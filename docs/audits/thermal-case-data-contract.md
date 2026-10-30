# Thermal-only case data contract (A1)

**Issue:** #1150  
**Contract version:** `THERMAL_CASE_SCHEMA_VERSION = 1`  
**Status:** typed in-memory contract and validation rules defined; persistence adapter not yet implemented.

## 1. Design boundary

The thermal case is a solver-independent physics description. It must not require a momentum equation, pressure-velocity coupling, turbulence model, CFD runtime state, or an initialized velocity field.

The implementation reuses `ThermophysicalProperties` from `thermophysical_models.h`. It does not create a second conductivity / heat-capacity / density law system. The thermal material validator intentionally requires only (k(T)), (c_p(T)), and (ho(T,p)); viscosity and turbulent Prandtl are not prerequisites for a thermal-only fluid.

The data model is in `src/cfdx/physics/thermal_case.h`:

- `ThermalCase`: version, time mode, collections and whole-case validation.
- `ThermalMaterial`: shared thermophysical properties and declared validity range.
- `ThermalRegion`: solid or meshed fluid, with explicit transport hypothesis.
- `ThermalNode`: homogeneous 0D fluid volume, with mass derived from (ho(T,p)V).
- `ThermalInterface`: explicit two-region interface and contact-resistance model.
- `ThermalLink`: conductance connection between a region patch and/or a lumped node.

This is a data contract, not a solver. It does not assemble equations or claim that the new physical capabilities are implemented.

## 2. Physical semantics

### Regions

| Kind | Model | Contract |
|---|---|---|
| solid | n/a | A mesh reference and material are mandatory; fluid transport fields are forbidden |
| fluid | `stagnant_diffusion` | Thermal diffusion only; no velocity, pressure, momentum or mass-flux field is needed |
| fluid | `prescribed_mass_flux` | Requires the name of an authoritative face mass-flux field; it is not reconstructed from cell velocity |

Convection/natural circulation is not implied by `stagnant_diffusion`. A prescribed mass flux is metadata for a later adapter/solver stage; A1 does not read or validate the actual field against a mesh.

### 0D nodes

A node stores volume, material, initial absolute temperature and pressure. Its mass is derived from the selected material density at the initial state. A1 represents closed lumped thermal capacitance. Mass-flow network transport and enthalpy advection are deliberately deferred to A5; there is no implicit mass flow in the current model.

### Interfaces and links

- A conjugate-conduction interface must connect one solid region and one fluid region and identify both patches.
- Contact resistance is expressed in (m^2 K/W), non-negative, and cannot be silently combined with the zero-resistance conjugate model.
- A thermal link uses a positive conductance in (W/K). Its sign convention is (Q_{source\to target}=G(T_{source}-T_{target})).
- Region endpoints identify a patch; lumped-node endpoints must not carry a patch.
- A1 does not yet define radiation view factors, convection correlations, non-conformal mapping, or mass-flow links.

## 3. Units and numeric policy

All values in the typed contract are SI:

| Quantity | Unit |
|---|---|
| temperature | K (absolute) |
| pressure | Pa (absolute) |
| volume | m³ |
| density | kg/m³ |
| conductivity | W/(m·K) |
| specific heat | J/(kg·K) |
| contact resistance | m²·K/W |
| conductance | W/K |
| time step | s |
| mass flux field | kg/(m²·s), oriented positive from face owner to neighbour / outward at boundary as defined by the mesh flux convention |

The contract stores canonical SI values. Unit conversion belongs in input adapters; the core validator never guesses units or silently rescales values. NaN, infinity, non-positive absolute temperature/pressure/volume, invalid property ranges, missing references and unsupported schema versions are rejected.

For a steady case, `time_step` must be exactly zero. For a transient case it must be finite and strictly positive. This avoids ambiguous defaults at the solver boundary.

## 4. HDF5 mapping target

The case remains the canonical `.cfdx.h5` setup artifact; transient solution state belongs in `.dat.h5`. A1 defines the intended thermal subtree without changing the global CFDX case schema version:

```text
/physics/thermal/
    schema_version                 uint32 = 1
    time_mode                      "steady" | "transient"
    time_step                      float64 [s]
    materials/<id>/
        name
        reference_temperature      float64 [K]
        reference_pressure         float64 [Pa]
        temperature_range          float64[2] [K]
        properties/
            conductivity/          existing scalar-property model encoding
            heat_capacity/          existing scalar-property model encoding
            density/                existing density-model encoding
    regions/<id>/
        name
        kind                       "solid" | "fluid"
        mesh_id
        material_id
        initial_temperature        float64 [K]
        fluid_model                "stagnant_diffusion" | "prescribed_mass_flux"
        mass_flux_field            optional field reference
    nodes/<id>/
        name, material_id
        volume                     float64 [m3]
        initial_temperature        float64 [K]
        pressure                   float64 [Pa]
    interfaces/<id>/
        name, region_a, patch_a, region_b, patch_b
        model                      "conjugate_conduction" | "contact_resistance"
        contact_resistance         float64 [m2.K/W]
    links/<id>/
        name, source, target
        endpoint kinds and patch references
        model                      "conductance"
        conductance                float64 [W/K]
```

Stable IDs are references, not array positions. Names are human-readable labels and must be unique within each collection. Required attributes/datasets and their shapes/types are validated by the future HDF5 adapter. Missing mandatory values, unknown enum strings and unsupported versions must fail closed. No HDF5 reader/writer is claimed by A1.

## 5. Validation performed by A1

`ThermalCase::validate()` checks:

1. supported schema version and time-mode/time-step consistency;
2. non-empty case, material collection and physical entities;
3. unique material, region, node, interface and link names;
4. finite positive (k,c_p,ho) at reference and declared temperature-range endpoints;
5. region material references and initial temperature range;
6. solid/fluid transport compatibility and mandatory prescribed face-flux reference;
7. positive 0D volume and valid initial thermodynamic state;
8. interface region references, patch identifiers and solid-fluid pairing for conjugate conduction;
9. link endpoint types/references, patch semantics and positive finite conductance.

A1 validation is intentionally structural and thermophysical. Mesh existence, patch existence, face-flux dimensions/conservation, field orientation and HDF5 round-trip checks require the later mesh/HDF5 adapter and are not represented as passed here.

## 6. Deferred work

- A2: typed case adapter/CLI and stationary/transient autonomous driver.
- A3: mesh-aware conformal CHT assembly and equal/opposite interface heat-flux audit.
- A4: closed 0D node energy equation and independent analytical transient oracle.
- A5: mass-flow network and conservative enthalpy transport.
- A6: imported prescribed face-flux validation and advection-diffusion.
- A7: S2S integration and separately qualified extensions.

## 7. Test target

`test_thermal_case_contracts` covers valid mixed entities, viscosity independence, transient-step rules, prescribed-flux requirements, solid/fluid incompatibility, unknown material references, invalid lumped volume, CHT pairing, link validation and duplicate material rejection. Test registration is added to CMake; execution status must be read from CI on the PR head.
