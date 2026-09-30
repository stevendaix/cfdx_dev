// M0.7-T01 — Gauss gradient
//
// Spécification CFDX v0.7 §28 :
//   Première implémentation : Gauss linear.
//
// Pour un champ scalaire φ (Field<double, CELL>) :
//   ∇φ_c = (1/V_c) Σ_f φ_f Sf_f
//
//   - φ_f est interpolé géométriquement à partir des cellules owner/neighbour.
//     φ_f = (1-w) * φ_owner + w * φ_neighbour
//     w = |Cf - C_owner| / |C_neighbour - C_owner|  (pondération géométrique)
//   - Pour les faces de frontière, on utilise φ_f = φ_owner (zero-gradient).
//
// Le résultat est un Field<double, CELL> de dimension 3 (le gradient est un vecteur
// 3D par cellule). Le stockage est SoA (x, y, z contigus).

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/ownership.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/geometry/geometry_cache.h"
#include "cfdx/core/fvm/least_squares_gradient.h"
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>
#include <algorithm>

namespace cfdx {
namespace core {

// Calcule le poids d'interpolation géométrique pour une face.
// w = |Cf - C_owner| / |C_neighbour - C_owner|
// Pour face frontière : w = 0 (zero gradient)
inline double geometric_interpolation_weight(
    const Vec3& face_centre,
    const Vec3& owner_centre,
    const Vec3* neighbour_centre)
{
    if (!neighbour_centre) return 0.0;  // Boundary face

    const Vec3 d_owner = face_centre - owner_centre;
    const Vec3 d_neigh = *neighbour_centre - owner_centre;
    const double dist_owner = d_owner.mag();
    const double dist_total = d_neigh.mag();

    if (!std::isfinite(dist_total) || dist_total < 1e-15)
        throw std::runtime_error("geometric_interpolation_weight: degenerate cell-centre distance");
    return std::min(1.0, std::max(0.0, dist_owner / dist_total));
}

// Calcule le gradient d'un champ scalaire cellulaire par la méthode Gauss.
//
inline Field<double, Location::CELL> compute_gradient_gauss(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    const GeometryCache& geometry)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();

    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("compute_gradient_gauss: invalid geometry cache");
    if (cell_field.size() != n_cells) {
        throw std::runtime_error("compute_gradient_gauss: field size != n_cells");
    }
    if (cell_field.dimension() != 1) {
        throw std::runtime_error("compute_gradient_gauss: field must be scalar (dim=1)");
    }

    Field<double, Location::CELL> grad(n_cells, cell_field.name() + "_grad", cell_field.metadata().unit + "/m", 3);

    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* cell_values = cell_field.component_data(0);

    std::vector<double> face_field_values(n_faces, 0.0);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        if (owner >= n_cells)
            throw std::runtime_error("compute_gradient_gauss: owner index out of range");
        const std::int64_t neighbour = own.neighbour(f);
        if (neighbour >= 0) {
            const std::size_t nb = static_cast<std::size_t>(neighbour);
            if (nb >= n_cells)
                throw std::runtime_error("compute_gradient_gauss: neighbour index out of range");
            const double w = geometric_interpolation_weight(
                geometry.face_centres[f], geometry.cell_centres[owner],
                &geometry.cell_centres[nb]);
            face_field_values[f] =
                (1.0 - w) * cell_values[owner] + w * cell_values[nb];
        } else {
            face_field_values[f] = cell_values[owner];
        }
    }

    double* gx = grad.component_data(0);
    double* gy = grad.component_data(1);
    double* gz = grad.component_data(2);
    for (std::size_t c = 0; c < n_cells; ++c) {
        Vec3 sum;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const Vec3 Sf_cell = (own.owner(f) == c)
                ? geometry.face_Sf[f]
                : geometry.face_Sf[f] * (-1.0);
            sum = sum + Sf_cell * face_field_values[f];
        }
        const double volume = geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error("compute_gradient_gauss: non-positive cell volume");
        const double inv_volume = 1.0 / volume;
        gx[c] = sum.x * inv_volume;
        gy[c] = sum.y * inv_volume;
        gz[c] = sum.z * inv_volume;
    }
    return grad;
}

inline Field<double, Location::CELL> compute_gradient_gauss(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return compute_gradient_gauss(cell_field, mesh, geometry);
}


inline Field<double, Location::CELL> compute_gradient_least_squares(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh)
{
    const std::size_t n_cells = mesh.n_cells();
    if (cell_field.size() != n_cells)
        throw std::runtime_error("compute_gradient_least_squares: field size != n_cells");
    if (cell_field.dimension() != 1)
        throw std::runtime_error("compute_gradient_least_squares: field must be scalar (dim=1)");

    Field<double, Location::CELL> grad(
        n_cells, cell_field.name() + "_grad_ls", cell_field.metadata().unit + "/m", 3);

    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const auto& own = mesh.ownership();
    const auto geometry = make_geometry_cache(mesh);
    const double* values = cell_field.component_data(0);

    for (std::size_t c = 0; c < n_cells; ++c) {
        std::vector<Vec3> neighbours;
        std::vector<double> neighbour_values;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const std::size_t owner = own.owner(f);
            const auto nraw = own.neighbour(f);
            // The adjacent cell can be either the stored owner or the stored
            // neighbour depending on which side of the face cell c lies. Using
            // only faces owned by c drops every higher-index neighbour and
            // degrades the stencil to a one-sided, first-order operator.
            std::size_t nb = n_cells;
            if (owner == c) {
                if (nraw >= 0) nb = static_cast<std::size_t>(nraw);
            } else {
                nb = owner;
            }
            if (nb >= n_cells || nb == c)
                continue;
            neighbours.push_back(geometry.cell_centres[nb]);
            neighbour_values.push_back(values[nb]);
        }

        const Vec3 g = least_squares_gradient(
            geometry.cell_centres[c], values[c], neighbours, neighbour_values);
        grad(c, 0) = g.x;
        grad(c, 1) = g.y;
        grad(c, 2) = g.z;
    }
    return grad;
}

// Vertex-based (secondary) Green-Gauss gradient.
//
// Cell values are interpolated to the mesh vertices with inverse-distance
// weighting, vertex values are averaged onto each face, and the gradient is
// then the Green-Gauss sum:
//   phi_v   = sum_c (phi_c / |C_c - P_v|) / sum_c (1/|C_c - P_v|)
//   phi_f   = average of phi_v over the face vertices
//   grad_c  = (1/V_c) sum_f phi_f Sf_c
//
// This is the classic "Green-Gauss node based" variant. It is linear-exact on
// affine meshes (the weighted vertex average is exact there) but does not, in
// general, reproduce quadratic fields to machine precision, so its
// verification requirements differ from the cell-based scheme.
inline Field<double, Location::CELL> compute_gradient_gauss_vertex(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    const GeometryCache& geometry)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();
    const std::size_t n_points = mesh.n_points();

    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("compute_gradient_gauss_vertex: invalid geometry cache");
    if (cell_field.size() != n_cells)
        throw std::runtime_error("compute_gradient_gauss_vertex: field size != n_cells");
    if (cell_field.dimension() != 1)
        throw std::runtime_error("compute_gradient_gauss_vertex: field must be scalar (dim=1)");

    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const auto* face_vertices = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* cell_values = cell_field.component_data(0);

    // Unique vertices per cell (a hexahedral cell touches each vertex through
    // several faces, so dedupe before weighting).
    std::vector<std::vector<std::size_t>> cell_vertices(n_cells);
    {
        std::vector<std::size_t> scratch;
        std::vector<char> seen(n_points, 0);
        for (std::size_t c = 0; c < n_cells; ++c) {
            scratch.clear();
            for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
                const std::size_t f = cell_faces[k];
                for (Offset o = face_offsets[f]; o < face_offsets[f + 1]; ++o) {
                    const std::size_t v = face_vertices[o];
                    if (!seen[v]) { seen[v] = 1; scratch.push_back(v); }
                }
            }
            cell_vertices[c] = scratch;
            for (const std::size_t v : scratch) seen[v] = 0;
        }
    }

    // Inverse-distance vertex averaging from the neighbouring cells.
    std::vector<std::vector<std::pair<std::size_t, double>>> vertex_cells(n_points);
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Vec3 cc = geometry.cell_centres[c];
        for (const std::size_t v : cell_vertices[c]) {
            const double dx = mesh.points().x(v) - cc.x;
            const double dy = mesh.points().y(v) - cc.y;
            const double dz = mesh.points().z(v) - cc.z;
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (!(r > 1e-14) || !std::isfinite(r))
                throw std::runtime_error("compute_gradient_gauss_vertex: degenerate vertex-cell distance");
            vertex_cells[v].emplace_back(c, 1.0 / r);
        }
    }

    std::vector<double> vertex_values(n_points, 0.0);
    std::vector<double> vertex_weights(n_points, 0.0);
    for (std::size_t v = 0; v < n_points; ++v) {
        double acc = 0.0;
        for (const auto& [c, w] : vertex_cells[v]) {
            acc += w * cell_values[c];
            vertex_weights[v] += w;
        }
        if (vertex_weights[v] > 0.0) vertex_values[v] = acc / vertex_weights[v];
    }

    // Average vertex values onto each face.
    std::vector<double> face_values(n_faces, 0.0);
    for (std::size_t f = 0; f < n_faces; ++f) {
        double acc = 0.0;
        std::size_t count = 0;
        for (Offset o = face_offsets[f]; o < face_offsets[f + 1]; ++o, ++count)
            acc += vertex_values[face_vertices[o]];
        face_values[f] = count ? acc / static_cast<double>(count) : 0.0;
    }

    // Green-Gauss assembly, identical to the cell-based operator.
    Field<double, Location::CELL> grad(
        n_cells, cell_field.name() + "_grad_vg", cell_field.metadata().unit + "/m", 3);
    double* gx = grad.component_data(0);
    double* gy = grad.component_data(1);
    double* gz = grad.component_data(2);
    for (std::size_t c = 0; c < n_cells; ++c) {
        Vec3 sum;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const Vec3 Sf = (own.owner(f) == c)
                ? geometry.face_Sf[f]
                : geometry.face_Sf[f] * (-1.0);
            sum = sum + Sf * face_values[f];
        }
        const double volume = geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error("compute_gradient_gauss_vertex: non-positive cell volume");
        const double inv_volume = 1.0 / volume;
        gx[c] = sum.x * inv_volume;
        gy[c] = sum.y * inv_volume;
        gz[c] = sum.z * inv_volume;
    }
    return grad;
}

inline Field<double, Location::CELL> compute_gradient_gauss_vertex(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return compute_gradient_gauss_vertex(cell_field, mesh, geometry);
}

// Point-linear (skew-corrected) Green-Gauss gradient.
//
// The two-point Gauss interpolation evaluates the face value on the
// owner->neighbour centre line, which is NOT the face centroid on a generally
// shaped polyhedron (e.g. a tetrahedron): the resulting O(1) interpolation
// bias makes the operator inconsistent for linear fields (see the #461
// investigation note in NUMERICAL_METHOD_CONTRACTS.md).
//
// The point-linear scheme (Greenshields & Weller, Notes on CFD, sec. 3.15)
// builds phi_f from:
//   - vertex values, interpolated from the adjacent cells with inverse
//     distance weights;
//   - phi_fp, the value at the FACE POINT: the intersection of the centre line
//     with the face plane, obtained by linear interpolation along the line;
//   - a signed-area-weighted average over the triangles (v_j, v_{j+1}, fp).
// For a linear field every ingredient is exact, so the triangle-average equals
// phi at the face centroid and the gradient is linear-exact on any
// planar-faced mesh; on smooth fields the error is second order.
inline Field<double, Location::CELL> compute_gradient_gauss_point(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh,
    const GeometryCache& geometry)
{
    const std::size_t n_cells = mesh.n_cells();
    const std::size_t n_faces = mesh.n_faces();
    const std::size_t n_points = mesh.n_points();

    if (!is_valid(geometry, mesh))
        throw std::invalid_argument("compute_gradient_gauss_point: invalid geometry cache");
    if (cell_field.size() != n_cells)
        throw std::runtime_error("compute_gradient_gauss_point: field size != n_cells");
    if (cell_field.dimension() != 1)
        throw std::runtime_error("compute_gradient_gauss_point: field must be scalar (dim=1)");

    const auto* cell_faces = mesh.cells().faces_data();
    const auto* cell_offsets = mesh.cells().offsets_data();
    const auto* face_vertices = mesh.faces().vertices_data();
    const auto* face_offsets = mesh.faces().offsets_data();
    const FaceOwnership& own = mesh.ownership();
    const double* cell_values = cell_field.component_data(0);

    // Unique vertices per cell, then inverse-distance vertex averaging.
    std::vector<std::vector<std::size_t>> cell_vertices(n_cells);
    {
        std::vector<std::size_t> scratch;
        std::vector<char> seen(n_points, 0);
        for (std::size_t c = 0; c < n_cells; ++c) {
            scratch.clear();
            for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
                const std::size_t f = cell_faces[k];
                for (Offset o = face_offsets[f]; o < face_offsets[f + 1]; ++o) {
                    const std::size_t v = face_vertices[o];
                    if (!seen[v]) { seen[v] = 1; scratch.push_back(v); }
                }
            }
            cell_vertices[c] = scratch;
            for (const std::size_t v : scratch) seen[v] = 0;
        }
    }
    std::vector<double> vertex_values(n_points, 0.0);
    std::vector<double> vertex_weights(n_points, 0.0);
    for (std::size_t c = 0; c < n_cells; ++c) {
        const Vec3 cc = geometry.cell_centres[c];
        for (const std::size_t v : cell_vertices[c]) {
            const double dx = mesh.points().x(v) - cc.x;
            const double dy = mesh.points().y(v) - cc.y;
            const double dz = mesh.points().z(v) - cc.z;
            const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (!(r > 1e-14) || !std::isfinite(r))
                throw std::runtime_error("compute_gradient_gauss_point: degenerate vertex-cell distance");
            const double w = 1.0 / r;
            vertex_values[v] += w * cell_values[c];
            vertex_weights[v] += w;
        }
    }
    for (std::size_t v = 0; v < n_points; ++v)
        if (vertex_weights[v] > 0.0) vertex_values[v] /= vertex_weights[v];

    // Point-linear face reconstruction.
    std::vector<double> face_values(n_faces, 0.0);
    for (std::size_t f = 0; f < n_faces; ++f) {
        const std::size_t owner = own.owner(f);
        const std::int64_t nraw = own.neighbour(f);
        if (nraw < 0) {   // zero-gradient boundary policy
            face_values[f] = cell_values[owner];
            continue;
        }
        const std::size_t nb = static_cast<std::size_t>(nraw);
        const Vec3 P = geometry.cell_centres[owner];
        const Vec3 N = geometry.cell_centres[nb];
        const Vec3 Cf = geometry.face_centres[f];
        const Vec3 Sf = geometry.face_Sf[f];
        const double sn = Sf.mag();
        if (!(sn > 1e-14) || !std::isfinite(sn))
            throw std::runtime_error("compute_gradient_gauss_point: degenerate face normal");
        const Vec3 normal = Sf * (1.0 / sn);
        const Vec3 d = N - P;
        const double denom = normal.dot(d);

        const Offset o0 = face_offsets[f];
        const Offset on = face_offsets[f + 1];
        const std::size_t nv = static_cast<std::size_t>(on - o0);

        double t = 0.0;
        bool on_line = std::abs(denom) > 1e-14;
        if (on_line) {
            t = normal.dot(Cf - P) / denom;
            on_line = std::isfinite(t);
        }
        if (!on_line || nv < 3) {
            // Fallback: plain average of the inverse-distance vertex values.
            double acc = 0.0;
            for (Offset o = o0; o < on; ++o) acc += vertex_values[face_vertices[o]];
            face_values[f] = acc / static_cast<double>(nv);
            continue;
        }

        const double phi_fp = cell_values[owner] + t * (cell_values[nb] - cell_values[owner]);
        const Vec3 Xfp = P + d * t;

        double area = 0.0;
        double sum = 0.0;
        for (std::size_t j = 0; j < nv; ++j) {
            const std::size_t vj = face_vertices[o0 + j];
            const std::size_t vj1 = face_vertices[o0 + ((j + 1) % nv)];
            const Vec3 Xj{mesh.points().x(vj), mesh.points().y(vj), mesh.points().z(vj)};
            const Vec3 Xj1{mesh.points().x(vj1), mesh.points().y(vj1), mesh.points().z(vj1)};
            const double sa = 0.5 * ((Xj1 - Xj).cross(Xfp - Xj)).dot(normal);   // signed area
            if (std::abs(sa) <= 1e-30 * sn * sn) continue;
            area += sa;
            sum += sa * (vertex_values[vj] + vertex_values[vj1] + phi_fp) / 3.0;
        }
        face_values[f] = area != 0.0
            ? sum / area
            : 0.5 * (cell_values[owner] + cell_values[nb]);
    }

    // Green-Gauss assembly, identical to the cell-based operator.
    Field<double, Location::CELL> grad(
        n_cells, cell_field.name() + "_grad_ggp", cell_field.metadata().unit + "/m", 3);
    double* gx = grad.component_data(0);
    double* gy = grad.component_data(1);
    double* gz = grad.component_data(2);
    for (std::size_t c = 0; c < n_cells; ++c) {
        Vec3 sum;
        for (Offset k = cell_offsets[c]; k < cell_offsets[c + 1]; ++k) {
            const std::size_t f = cell_faces[k];
            const Vec3 Sf = (own.owner(f) == c)
                ? geometry.face_Sf[f]
                : geometry.face_Sf[f] * (-1.0);
            sum = sum + Sf * face_values[f];
        }
        const double volume = geometry.cell_volumes[c];
        if (!(volume > 0.0) || !std::isfinite(volume))
            throw std::runtime_error("compute_gradient_gauss_point: non-positive cell volume");
        gx[c] = sum.x / volume;
        gy[c] = sum.y / volume;
        gz[c] = sum.z / volume;
    }
    return grad;
}

inline Field<double, Location::CELL> compute_gradient_gauss_point(
    const Field<double, Location::CELL>& cell_field,
    const Mesh& mesh)
{
    const GeometryCache geometry = make_geometry_cache(mesh);
    return compute_gradient_gauss_point(cell_field, mesh, geometry);
}

}  // namespace core
}  // namespace cfdx
