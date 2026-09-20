# Boundary Patches — `.msh` PhysicalNames → Mesh BoundaryPatches

## Pattern — Gmsh `.msh` v2.2 + PhysicalNames

The `.msh` file (`mms_square_16x16.msh`) contains `$PhysicalNames` and sometimes `$PhysicalGroups`. The `gmsh_importer.cpp` must:

1. Read `.msh` content via `std::ifstream`
2. Call `parsePhysicalNames(file_content)` → `std::map<uint32_t, std::string>`
3. Call `parsePhysicalGroups(file_content)` → `std::map<uint32_t, std::vector<uint32_t>>`
4. Build boundary patches from the detected names (`left`, `right`, `bottom`, `top`, `wall`)

```cpp
std::map<std::string, std::vector<std::size_t>> boundary_patches_map;
for (const auto& [tag, name] : physical_names) {
    std::string clean_name = name;
    clean_name.erase(std::remove(clean_name.begin(), clean_name.end(), '"'), clean_name.end());
    clean_name.erase(std::remove(clean_name.begin(), clean_name.end(), ' '), clean_name.end());
    if (clean_name == "left" || clean_name == "right" ||
        clean_name == "bottom" || clean_name == "top" || clean_name == "wall") {
        boundary_patches_map[clean_name] = {};  // filled later with face_ids
    }
}
```

## Pattern — Geometric Face Distribution (2D Meshes)

For each face with `mesh.ownership().neighbour(f) == -1` (boundary face), calculate the geometric center:

```cpp
const auto off = mesh.faces().offsets_data()[f];
const auto n = mesh.faces().offsets_data()[f + 1] - off;
double cx = 0.0, cy = 0.0;
for (std::size_t v = 0; v < n; ++v) {
    std::size_t v_idx = mesh.faces().vertices_data()[off + v];
    cx += pts_ref.x(static_cast<std::size_t>(v_idx));
    cy += pts_ref.y(static_cast<std::size_t>(v_idx));
}
cx /= n; cy /= n;
```

Then distribute by position (`eps` threshold):
- `cx < eps` → `left`
- `cx > 1.0 - eps` → `right`
- `cy < eps` → `bottom`
- `cy > 1.0 - eps` → `top`

## Key Rules

- `mesh.boundary().add_patch(name, PatchType::WALL)` creates patches; `.name` is the string attribute (not `.name()` method).
- `mesh.boundary().patch(p).name` is accessed as `.name` (public member), `.size()` is `.n_faces()` or `.n_boundary_faces()` (method call).
- The `.msh` file must contain `$PhysicalNames` for boundary identification; without it (`mms_meshes` generated without PhysicalNames), `mesh.boundary().n_patches()` returns 0.
- `PatchType::WALL` (not `BOUNDARY`) is the correct enum value for boundary conditions.
- `mesh.boundary()` is rebuilt with `mesh.boundary() = patches;` after constructing `cfdx::core::BoundaryPatches`.

## Pitfalls

- **Pitfall**: `mesh.boundary().add_patch("left", PatchType::BOUNDARY)` — `BOUNDARY` is not a valid enum value; the correct value is `WALL`.
- **Pitfall**: Calling `.add_patch()` with duplicate names (`bottom` added 4 times) throws `BoundaryPatches: duplicate patch name 'bottom'`. Each unique patch name must be added exactly once, with all boundary face IDs collected per name.
- **Pitfall**: Using `eps = 1e-6` (absolute) is too strict for geometric face center comparison; `eps = 0.5 * h` (where `h = 1/sqrt(n_cells)`) is too loose. The correct value is `eps = 0.01 - 0.02` for square domains.
- **Pitfall**: Not reading the full `.msh` file before parsing `PhysicalNames`. The parser must seek back to the file start (`file.seekg(0)`) or open a new stream to read `$PhysicalNames` after reading `$Elements`.
