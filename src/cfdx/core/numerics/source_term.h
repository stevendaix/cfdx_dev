// M0.7-T06 — Source Term Linearization
//
// Spécification CFDX v0.7 §32 :
//   Terme source linéarisé : S(φ) = Su + Sp * φ
//   - Su : partie explicite (independent of φ)
//   - Sp : coefficient implicite (≤ 0 pour stabilité)
//   - Discrétisation : aP * φP = Σ aN * φN + Su + Sp * φP
//     → (aP - Sp) * φP = Σ aN * φN + Su
//
//   Utilisé pour linéariser termes sources non-linéaires :
//   - Réactions chimiques
//   - Rayonnement
//   - Sources d'énergie
//   - Termes de dissipation

#pragma once

#include "cfdx/core/field/field.h"
#include "cfdx/core/mesh/mesh.h"
#include "cfdx/core/mesh/index_types.h"
#include "cfdx/core/geometry/face_geometry.h"
#include "cfdx/core/geometry/cell_geometry.h"
#include <cstddef>
#include <functional>
#include <vector>

namespace cfdx {
namespace core {

// Structure contenant les coefficients de linéarisation du terme source.
// Su : partie explicite (source constante)
// Sp : coefficient implicite (multiplie par φ, doit être ≤ 0)
struct SourceTerm {
    Field<double, Location::CELL> Su;  // Explicit source
    Field<double, Location::CELL> Sp;  // Implicit coefficient (≤ 0)

    SourceTerm() = default;

    SourceTerm(std::size_t n_cells, const std::string& name = "source")
        : Su(n_cells, name + "_Su", "kg/m3/s", 1),
          Sp(n_cells, name + "_Sp", "1/s", 1)
    {
        Su.fill(0.0);
        Sp.fill(0.0);
    }

    // Create from explicit source only (Sp = 0)
    static SourceTerm explicit_only(std::size_t n_cells,
                                     const std::string& name,
                                     const Field<double, Location::CELL>& su) {
        SourceTerm st(n_cells, name);
        st.Su = su;
        return st;
    }

    // Create from implicit coefficient only (Su = 0)
    static SourceTerm implicit_only(std::size_t n_cells,
                                     const std::string& name,
                                     const Field<double, Location::CELL>& sp) {
        SourceTerm st(n_cells, name);
        st.Sp = sp;
        return st;
    }

    // Apply linearization to diagonal and RHS of linear system.
    // For: aP * φP = Σ aN * φN + Su + Sp * φP
    // Returns: new_diag = diag - Sp, new_rhs = rhs + Su
    inline void apply_to_system(
        std::vector<double>& diag,
        std::vector<double>& rhs) const
    {
        const std::size_t n = Su.size();
        const double* su_data = Su.component_data(0);
        const double* sp_data = Sp.component_data(0);

        for (std::size_t c = 0; c < n; ++c) {
            diag[c] -= sp_data[c];  // Sp ≤ 0, so -Sp ≥ 0 (increases diagonal dominance)
            rhs[c] += su_data[c];
        }
    }

    // Apply to single cell
    inline void apply_to_cell(std::size_t c, double& diag, double& rhs) const {
        const double* su_data = Su.component_data(0);
        const double* sp_data = Sp.component_data(0);
        diag -= sp_data[c];
        rhs += su_data[c];
    }
};

// Linearized source term for common cases.
// Returns SourceTerm with Su and Sp computed from φ.
using LinearizationFunc = std::function<SourceTerm(const Field<double, Location::CELL>&, const Mesh&)>;

// Constant source: S = Su_const (Sp = 0)
inline SourceTerm make_constant_source(std::size_t n_cells, double su_value, const std::string& name = "const_source") {
    SourceTerm st(n_cells, name);
    st.Su.fill(su_value);
    return st;
}

// Linear source: S = k * φ  (Sp = k, Su = 0)
// Note: k should be ≤ 0 for stability (implicit treatment)
inline SourceTerm make_linear_source(const Field<double, Location::CELL>& phi, double k, const std::string& name = "linear_source") {
    SourceTerm st(phi.size(), name);
    (void)phi; // size() used
    double* sp_data = st.Sp.component_data(0);

    for (std::size_t c = 0; c < phi.size(); ++c) {
        sp_data[c] = k;
    }
    return st;
}

// Decay source: S = -λ * φ  (Sp = -λ ≤ 0, Su = 0)
// Always stable (implicit)
inline SourceTerm make_decay_source(const Field<double, Location::CELL>& phi, double lambda, const std::string& name = "decay_source") {
    return make_linear_source(phi, -lambda, name);
}

// Linearized non-linear source: S(φ) = f(φ)
// Linearization around φ0: S(φ) ≈ f(φ0) + f'(φ0) * (φ - φ0)
// = [f(φ0) - f'(φ0)*φ0] + f'(φ0)*φ
// So: Su = f(φ0) - f'(φ0)*φ0, Sp = f'(φ0)  (clamped to ≤ 0)
inline SourceTerm linearize_source(const Field<double, Location::CELL>& phi,
                                    std::function<double(double)> f,
                                    std::function<double(double)> df,
                                    const std::string& name = "lin_source")
{
    SourceTerm st(phi.size(), name);
    const double* phi_data = phi.component_data(0);
    double* su_data = st.Su.component_data(0);
    double* sp_data = st.Sp.component_data(0);

    for (std::size_t c = 0; c < phi.size(); ++c) {
        const double phi_c = phi_data[c];
        const double f_val = f(phi_c);
        const double df_val = df(phi_c);
        const double sp_clamped = std::min(0.0, df_val);  // Ensure Sp ≤ 0

        su_data[c] = f_val - sp_clamped * phi_c;
        sp_data[c] = sp_clamped;
    }
    return st;
}

// Exponential source: S = A * exp(-B * φ)
// Linearization: Su = A*exp(-B*φ0) + A*B*φ0*exp(-B*φ0), Sp = -A*B*exp(-B*φ0)
inline SourceTerm make_exponential_source(const Field<double, Location::CELL>& phi,
                                           double A, double B,
                                           const std::string& name = "exp_source")
{
    return linearize_source(phi,
        [A, B](double phi_val) { return A * std::exp(-B * phi_val); },
        [A, B](double phi_val) { return -A * B * std::exp(-B * phi_val); },
        name);
}

// Power law source: S = A * φ^n
// Linearization: Su = A*φ0^n - n*A*φ0^n, Sp = n*A*φ0^(n-1)
inline SourceTerm make_power_source(const Field<double, Location::CELL>& phi,
                                     double A, double n,
                                     const std::string& name = "power_source")
{
    return linearize_source(phi,
        [A, n](double phi_val) { return A * std::pow(phi_val, n); },
        [A, n](double phi_val) { return n * A * std::pow(phi_val, n - 1.0); },
        name);
}

// Source term with spatial variation: S = f(x, y, z, φ)
// User provides Su_func(x,y,z) and Sp_func(x,y,z,φ)
inline SourceTerm make_spatial_source(const Mesh& mesh,
                                       std::function<double(const Vec3&)> su_func,
                                       std::function<double(const Vec3&, double)> sp_func,
                                       const std::string& name = "spatial_source")
{
    const std::size_t n_cells = mesh.n_cells();
    SourceTerm st(n_cells, name);

    const CellConnectivity& cells = mesh.cells();
    const auto* cell_faces = cells.faces_data();
    const auto* cell_offsets = cells.offsets_data();
    (void)cell_faces; // used in compute_cell_geometry

    // Compute cell centres for spatial functions
    std::vector<Vec3> cell_centres(n_cells);

    const FaceConnectivity& faces = mesh.faces();
    const auto* face_verts = faces.vertices_data();
    const auto* face_offsets = faces.offsets_data();
    const PointCloud& pts = mesh.points();
    (void)face_verts; // used in compute_face_geometry
    (void)pts;        // used in compute_face_geometry

    std::vector<Vec3> face_centres(mesh.n_faces());
    std::vector<Vec3> face_Sf(mesh.n_faces());

    for (std::size_t f = 0; f < mesh.n_faces(); ++f) {
        const Offset off = face_offsets[f];
        const Offset n = face_offsets[f + 1] - off;
        (void)n; // used in compute_face_geometry
        const FaceGeometry fg = compute_face_geometry(
            pts.x_data(), pts.y_data(), pts.z_data(), face_verts, off, n);
        face_centres[f] = fg.centre;
        face_Sf[f] = fg.Sf;
    }

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Offset off = cell_offsets[c];
        const Offset n = cell_offsets[c + 1] - off;
        (void)n; // used in compute_cell_geometry
        const CellGeometry cg = compute_cell_geometry_oriented(
            face_centres.data(), face_Sf.data(), cell_faces + off, n,
            static_cast<CellIndex>(c), mesh.ownership());
        cell_centres[c] = cg.centre;
    }

    double* su_data = st.Su.component_data(0);
    double* sp_data = st.Sp.component_data(0);

    for (std::size_t c = 0; c < n_cells; ++c) {
        const Vec3& cc = cell_centres[c];
        su_data[c] = su_func(cc);
        sp_data[c] = std::min(0.0, sp_func(cc, 0.0));  // Evaluate at φ=0, clamp Sp
    }

    return st;
}

// Helper: update Sp based on current φ (for iterative linearization)
// Su is recomputed as: Su = f(φ) - Sp * φ
// Sp is clamped to ≤ 0
inline void update_linearization(const Field<double, Location::CELL>& phi,
                                  SourceTerm& st,
                                  std::function<double(double)> f,
                                  std::function<double(double)> df)
{
    const double* phi_data = phi.component_data(0);
    double* su_data = st.Su.component_data(0);
    double* sp_data = st.Sp.component_data(0);

    for (std::size_t c = 0; c < phi.size(); ++c) {
        const double phi_c = phi_data[c];
        const double f_val = f(phi_c);
        const double df_val = df(phi_c);
        const double sp_clamped = std::min(0.0, df_val);

        su_data[c] = f_val - sp_clamped * phi_c;
        sp_data[c] = sp_clamped;
    }
}

// Compute source contribution to RHS vector (Su only, no Sp treatment)
// Useful for explicit source treatment
inline Field<double, Location::CELL> compute_source_rhs(const SourceTerm& st) {
    return st.Su;  // Copy of explicit source
}

}  // namespace core
}  // namespace cfdx