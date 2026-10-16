// ghia_pcd_audit_dump.cpp — diagnostic-only audit tool (#914).
//
// Branch: ci/long-campaign-diagnostics (never merged, PR #903). This file
// lives under .github/workflows/ so that pushing it never re-triggers the
// multi-hour long-campaign workflow (paths-ignore) and only wakes the
// ghia-pcd-audit workflow (paths filter).
//
// Purpose (issue #914, audit Fp-vs-Kp of the nested PCD action): dump the
// EXACT production coupled matrix of the ghia 32x32 Re=100 PCD case
// (tests/validation/test_ghia_cavity.cpp --pcd, post-#908 controls) at
// nonlinear iteration 2 — the first iteration with a non-zero Rhie-Chow mass
// flux, the regime where the Fp/Kp ordering actually matters — together with
// the PCD pressure operators (Mp/Kp/Fp) assembled from the same state, in
// round-trip bit-exact %a (C99 hexadecimal float) format. A FP64 twin can
// then replay the outer FGMRES solve with the action -Fp^{-1} vs -Kp^{-1}
// and count external iterations on the production operators.
//
// Method (pattern #872/#898 evidence):
//   Phase 0: cavity mesh 32x32 + controls, byte-for-byte the --pcd test
//            (Re=100, rho=1, nu=0.01, alpha_u=0.7, alpha_p=0.3, coupled
//            budget 5000, coupled tolerance 1e-10, FGMRES restart 4096,
//            pcd fp budget 200 / restart 512 — the #908 mirror —, schur_model
//            PCD, upwind, bounded convection, gauge cell 528).
//   Phase 1: reassemble the iteration-1 momentum equations from the exact
//            zero state (mass flux from U=0, grad p from p=0, zero body
//            force) — identical inputs, deterministic assembly, so this is
//            bit-identical to what production assembles before its first
//            coupled solve.
//   Phase 1b: audit-assemble the ITERATION-1 coupled system (verbatim
//            assembly body extracted from steady_incompressible_solver.h)
//            and validate against the production run-38809298638 log line
//            (M_l2=2.54234, C_l2=2.35617, schur_diag_min=0.0536722).
//   Phase 2: run exactly one production nonlinear iteration
//            (solve_steady_incompressible, callback stop after iteration 1,
//            ~20 s on the CI runner). U/p then hold the relaxed, gauged
//            post-iteration-1 state, and the production code path prints
//            its own COUPLED_MATRIX_SUMMARY for iteration 1 into the log.
//   Phase 3: replicate the production loop body to rebuild the iteration-2
//            assembly state: coupled_rAU = 1/ex1.diagonal, Rhie-Chow mass
//            flux from (U1,p1,rAU), grad p1, convected fields from U1,
//            momentum equations ex2/ey2/ez2.
//   Phase 4: audit-assemble the ITERATION-2 coupled system and validate
//            against the production iteration-2 log line (M_l2=2.85644,
//            C_l2=2.1722, schur_diag_min=0.0465948, nnz=76991).
//   Phase 5: assemble the PCD pressure operators from the same state
//            (production passes U_old=U1 to assemble_pcd_pressure_operators).
//   Phase 6: dump everything in %a format (chunked CSR so every file stays
//            below the ~1 MB GitHub contents-API read limit) + FNV-1a
//            checksums + meta.
//
// Compile with the production flags (CMakeLists: CXX_EXTENSIONS OFF,
// CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG"): the iteration-1 coupled solve is
// accepted at relative 9.76e-11 against a 1e-10 gate, so any FP-rounding
// divergence (e.g. -ffp-contract or -O3 differences) would invalidate the
// trajectory. g++ -std=c++17 -O2 -DNDEBUG matches the Release build.

#include "cfdx/physics/steady_incompressible_solver.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace cfdx::core;
using namespace cfdx::physics;

namespace {

// ---------------------------------------------------------------------------
// Cavity mesh — verbatim copy of make_cavity_mesh() from
// tests/validation/test_ghia_cavity.cpp (32x32 qualification grid).
// ---------------------------------------------------------------------------
Mesh make_cavity_mesh(std::size_t nx, std::size_t ny)
{
    if (nx < 4 || ny < 4) throw std::invalid_argument("cavity mesh is too small");
    Mesh mesh;
    mesh.points().resize((nx + 1) * (ny + 1) * 2);
    const auto point_id = [nx](std::size_t i, std::size_t j, std::size_t k) {
        return (j * (nx + 1) + i) * 2 + k;
    };
    for (std::size_t j=0;j<=ny;++j) for (std::size_t i=0;i<=nx;++i) {
        const double x=static_cast<double>(i)/static_cast<double>(nx);
        const double y=static_cast<double>(j)/static_cast<double>(ny);
        mesh.points().set(point_id(i,j,0),x,y,0.0);
        mesh.points().set(point_id(i,j,1),x,y,1.0);
    }

    std::map<std::vector<std::size_t>,std::size_t> face_map;
    std::vector<std::vector<std::size_t>> cell_faces(nx*ny);
    // Ownership is sized after face creation because shared faces are deduplicated.
    auto add_or_get_face=[&](std::initializer_list<std::size_t> vertices,std::size_t cell) {
        std::vector<std::size_t> sorted(vertices);
        std::sort(sorted.begin(),sorted.end());
        const auto found=face_map.find(sorted);
        if(found!=face_map.end()) {
            mesh.ownership().set_neighbour(found->second,static_cast<int>(cell));
            return found->second;
        }
        const std::size_t face=mesh.faces().n_faces();
        mesh.faces().push_face(vertices);
        mesh.ownership().resize(mesh.faces().n_faces());
        face_map.emplace(std::move(sorted),face);
        mesh.ownership().set_owner(face,cell);
        mesh.ownership().set_neighbour(face,FaceOwnership::BOUNDARY);
        return face;
    };

    for(std::size_t j=0;j<ny;++j) for(std::size_t i=0;i<nx;++i) {
        const std::size_t cell=j*nx+i;
        const std::size_t a=point_id(i,j,0), b=point_id(i+1,j,0);
        const std::size_t c=point_id(i+1,j+1,0), d=point_id(i,j+1,0);
        const std::size_t e=point_id(i,j,1), f=point_id(i+1,j,1);
        const std::size_t g=point_id(i+1,j+1,1), h=point_id(i,j+1,1);
        cell_faces[cell]={
            add_or_get_face({a,d,c,b},cell), add_or_get_face({e,f,g,h},cell),
            add_or_get_face({a,b,f,e},cell), add_or_get_face({d,h,g,c},cell),
            add_or_get_face({a,e,h,d},cell), add_or_get_face({b,c,g,f},cell)};
    }
    for(const auto& faces:cell_faces) mesh.cells().push_cell(faces);

    Patch bottom; bottom.name="bottom"; bottom.type=PatchType::WALL;
    Patch top; top.name="top"; top.type=PatchType::WALL;
    Patch left; left.name="left"; left.type=PatchType::WALL;
    Patch right; right.name="right"; right.type=PatchType::WALL;
    Patch front; front.name="front"; front.type=PatchType::WALL;
    Patch back; back.name="back"; back.type=PatchType::WALL;

    for(std::size_t f=0;f<mesh.n_faces();++f) {
        if(mesh.ownership().neighbour(f)>=0) continue;
        const auto& all_vertices=mesh.faces().vertices();
        const auto begin=all_vertices.begin()+static_cast<std::ptrdiff_t>(mesh.faces().face_offset(f));
        const auto end=begin+static_cast<std::ptrdiff_t>(mesh.faces().face_size(f));
        double x=0,y=0,z=0;
        for(auto it=begin;it!=end;++it) { x+=mesh.points().x(*it); y+=mesh.points().y(*it); z+=mesh.points().z(*it); }
        const double nverts=static_cast<double>(mesh.faces().face_size(f));
        x/=nverts; y/=nverts; z/=nverts;
        constexpr double tol=1e-12;
        if(std::abs(y)<tol) bottom.face_ids.push_back(f);
        else if(std::abs(y-1)<tol) top.face_ids.push_back(f);
        else if(std::abs(x)<tol) left.face_ids.push_back(f);
        else if(std::abs(x-1)<tol) right.face_ids.push_back(f);
        else if(std::abs(z)<tol) front.face_ids.push_back(f);
        else if(std::abs(z-1)<tol) back.face_ids.push_back(f);
        else throw std::runtime_error("cavity: unclassified boundary face");
    }
    mesh.boundary().add_patch(bottom); mesh.boundary().add_patch(top);
    mesh.boundary().add_patch(left); mesh.boundary().add_patch(right);
    mesh.boundary().add_patch(front); mesh.boundary().add_patch(back);
    return mesh;
}

// ---------------------------------------------------------------------------
// Round-trip bit-exact dump helpers (%a C99 hexadecimal floats).
// ---------------------------------------------------------------------------

std::FILE* audit_open(const std::string& path)
{
    std::FILE* f = std::fopen(path.c_str(), "w");
    if (!f) throw std::runtime_error("audit dump: cannot open " + path);
    return f;
}

std::uint64_t fnv1a_init() { return 1469598103934665603ULL; }

std::uint64_t fnv1a_bytes(std::uint64_t h, const void* data, std::size_t size)
{
    const auto* p = static_cast<const unsigned char*>(data);
    for (std::size_t i = 0; i < size; ++i) {
        h ^= p[i];
        h *= 1099511628211ULL;
    }
    return h;
}

// %a round-trips every finite double exactly (C99). Values are written one
// per line; the twin verifies its parse against the FNV-1a checksums.
void dump_vector_file(const std::string& path, const std::string& name,
                      const double* values, std::size_t n)
{
    std::FILE* f = audit_open(path);
    std::fprintf(f, "CFDX_VEC name=%s size=%llu\n",
                 name.c_str(), static_cast<unsigned long long>(n));
    for (std::size_t i = 0; i < n; ++i)
        std::fprintf(f, "%a\n", values[i]);
    std::fclose(f);
}

template <cfdx::core::Location L>
void dump_field_file(const std::string& path, const std::string& name,
                     const cfdx::core::Field<double, L>& field)
{
    std::FILE* f = audit_open(path);
    std::fprintf(f, "CFDX_FIELD name=%s components=%llu size=%llu\n",
                 name.c_str(),
                 static_cast<unsigned long long>(field.dimension()),
                 static_cast<unsigned long long>(field.size()));
    for (std::size_t d = 0; d < field.dimension(); ++d) {
        std::fprintf(f, "COMP %llu\n", static_cast<unsigned long long>(d));
        const double* v = field.component_data(d);
        for (std::size_t i = 0; i < field.size(); ++i)
            std::fprintf(f, "%a\n", v[i]);
    }
    std::fclose(f);
}

// Chunked CSR: every chunk file stays far below the ~1 MB GitHub contents
// API read limit. A chunk carries the global row_offsets slice for its row
// range plus the columns and values of those rows, so the concatenation of
// all chunks (in order) is the full CSR triple.
void dump_csr_chunks(const std::string& prefix, const std::string& name,
                     const SparseMatrix& m, std::size_t rows_per_chunk)
{
    const std::size_t nr = m.n_rows();
    const std::size_t nc = m.n_cols();
    const std::size_t nnz = m.nnz();
    const auto* ro = m.row_offsets_data();
    const auto* co = m.columns_data();
    const double* va = m.values_data();
    const std::size_t chunks = (nr + rows_per_chunk - 1) / rows_per_chunk;
    for (std::size_t c = 0; c < chunks; ++c) {
        const std::size_t r0 = c * rows_per_chunk;
        const std::size_t r1 = std::min(nr, r0 + rows_per_chunk);
        char idx[32];
        std::snprintf(idx, sizeof(idx), "%04llu", static_cast<unsigned long long>(c));
        std::FILE* f = audit_open(prefix + "." + idx);
        std::fprintf(f, "CFDX_CSR_CHUNK name=%s chunk=%llu chunks=%llu row0=%llu row_end=%llu n_rows=%llu n_cols=%llu nnz=%llu\n",
                     name.c_str(),
                     static_cast<unsigned long long>(c),
                     static_cast<unsigned long long>(chunks),
                     static_cast<unsigned long long>(r0),
                     static_cast<unsigned long long>(r1),
                     static_cast<unsigned long long>(nr),
                     static_cast<unsigned long long>(nc),
                     static_cast<unsigned long long>(nnz));
        std::fprintf(f, "OFFSETS");
        for (std::size_t i = r0; i <= r1; ++i)
            std::fprintf(f, " %llu", static_cast<unsigned long long>(ro[i]));
        std::fprintf(f, "\nCOLUMNS");
        for (std::size_t k = ro[r0]; k < ro[r1]; ++k)
            std::fprintf(f, " %llu", static_cast<unsigned long long>(co[k]));
        std::fprintf(f, "\nVALUES %llu\n",
                     static_cast<unsigned long long>(ro[r1] - ro[r0]));
        for (std::size_t k = ro[r0]; k < ro[r1]; ++k)
            std::fprintf(f, "%a\n", va[k]);
        std::fclose(f);
    }
}

// ---------------------------------------------------------------------------
// Summary recomputation (same formulas as the COUPLED_MATRIX_SUMMARY block
// of the production assembly body) so the tool can validate its own dumps
// against the production CI log values without parsing stderr.
// ---------------------------------------------------------------------------

struct AuditBlockNorms {
    double m_l2 = 0.0, g_l2 = 0.0, d_l2 = 0.0, c_l2 = 0.0;
    std::size_t nnz = 0;
};

AuditBlockNorms audit_block_norms(const SparseMatrix& A, std::size_t nv)
{
    AuditBlockNorms out;
    double m2 = 0.0, g2 = 0.0, d2 = 0.0, c2 = 0.0;
    const auto* ro = A.row_offsets_data();
    const auto* co = A.columns_data();
    const auto* va = A.values_data();
    const std::size_t n = A.n_rows();
    for (std::size_t r = 0; r < n; ++r) {
        for (std::size_t k = ro[r]; k < ro[r + 1]; ++k) {
            const double value = va[k];
            const std::size_t col = co[k];
            out.nnz++;
            if (r < nv && col < nv) m2 += value * value;
            else if (r < nv && col >= nv) g2 += value * value;
            else if (r >= nv && col < nv) d2 += value * value;
            else c2 += value * value;
        }
    }
    out.m_l2 = std::sqrt(m2);
    out.g_l2 = std::sqrt(g2);
    out.d_l2 = std::sqrt(d2);
    out.c_l2 = std::sqrt(c2);
    return out;
}

// Effective algebraic Schur diagonal minimum — verbatim logic of the summary
// block (row-wise D Auu^-1 G accumulation plus the pressure diagonal).
double audit_schur_effective_diag_min(const SparseMatrix& A, std::size_t nv)
{
    const std::size_t n = A.n_rows();
    const std::size_t np = n - nv;
    const auto* ro = A.row_offsets_data();
    const auto* co = A.columns_data();
    const auto* va = A.values_data();
    double schur_diag_min = std::numeric_limits<double>::infinity();
    for (std::size_t c = 0; c < np; ++c) {
        const std::size_t pr = nv + c;
        double effective = 0.0;
        for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k) {
            const std::size_t u = co[k];
            if (u >= nv) continue;
            const double d = va[k];
            double inv_m = 0.0;
            bool found = false;
            for (std::size_t mk = ro[u]; mk < ro[u + 1]; ++mk) {
                if (co[mk] == u) {
                    const double diag = va[mk];
                    if (std::isfinite(diag) && std::abs(diag) > 1e-30) {
                        inv_m = 1.0 / diag;
                        found = true;
                    }
                    break;
                }
            }
            if (!found) continue;
            for (std::size_t gk = ro[u]; gk < ro[u + 1]; ++gk) {
                const std::size_t q = co[gk];
                if (q == pr)
                    effective -= d * inv_m * va[gk];
            }
        }
        for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k)
            if (co[k] == pr) effective += va[k];
        schur_diag_min = std::min(schur_diag_min, effective);
    }
    return schur_diag_min;
}

bool audit_close(double measured, double expected)
{
    const double scale = std::max(1.0, std::abs(expected));
    return std::abs(measured - expected) <= 1e-4 * scale;
}

// ---------------------------------------------------------------------------
// Audit wrapper around the VERBATIM production assembly body of
// solve_coupled_momentum_continuity (steady_incompressible_solver.h).
// The body below is byte-identical to the production function between its
// opening and the COUPLED_MATRIX_SUMMARY block; the only additions are the
// parameter list and the return statement.
// ---------------------------------------------------------------------------

struct CoupledAuditSystem {
    SparseMatrix A;
    Vector b;
    Vector x0;
    std::vector<double> schur_diagonal;
};

CoupledAuditSystem audit_assemble_coupled(
    const Mesh& mesh,
    const FvGeometry& geometry,
    const ScalarEquation& ex,
    const ScalarEquation& ey,
    const ScalarEquation& ez,
    const Field<double, Location::CELL>& U_old,
    const Field<double, Location::CELL>& p_old,
    const VelocityBoundaryConditions& velocity_bcs,
    const ScalarBoundaryConditions& pressure_bcs,
    double rho,
    std::size_t reference_cell,
    double reference_value,
    const IncompressibleSolverControls& controls,
    const DiagnosticsControls& diagnostics)
{
    using namespace cfdx::core;
    const std::size_t nc = mesh.n_cells();
    const std::size_t nv = 3 * nc;
    const std::size_t np = nc;
    const std::size_t n = nv + np;

    SparseMatrix A(n, n);
    Vector b(n, 0.0);

    // The momentum equations received the current pressure explicitly during
    // segregated assembly. Undo that explicit contribution here: the coupled
    // matrix contains the pressure gradient implicitly instead.
    const auto grad_p_old =
        gauss_gradient_with_boundary(p_old, mesh, geometry, pressure_bcs);
    auto add_momentum_block = [&](const ScalarEquation& eq, std::size_t component) {
        const std::size_t row_base = component * nc;
        const auto* ro = eq.matrix.row_offsets_data();
        const auto* co = eq.matrix.columns_data();
        const auto* va = eq.matrix.values_data();
        const double* old_component = U_old.component_data(component);
        for (std::size_t r = 0; r < nc; ++r) {
            const double diagonal = eq.diagonal[r];
            if (!(diagonal > 0.0) || !std::isfinite(diagonal))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: invalid momentum diagonal at cell " +
                    std::to_string(r) + " component " + std::to_string(component));

            bool has_diagonal = false;
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t column = component * nc + co[k];
                A.push_back(row_base + r, column, va[k]);
                has_diagonal = has_diagonal || (column == row_base + r);
            }
            // The scalar equation keeps its authoritative diagonal separately
            // from the CSR entries. If a transport row reaches the coupled
            // assembler without an explicit diagonal entry, restore that
            // coefficient from the authoritative metadata instead of creating
            // an empty coupled momentum row.
            if (!has_diagonal)
                A.push_back(row_base + r, row_base + r, diagonal);
            b(row_base + r) =
                eq.rhs(r) + grad_p_old.component_data(component)[r] * geometry.cell_volumes[r];
        }
    };
    add_momentum_block(ex, 0);
    add_momentum_block(ey, 1);
    add_momentum_block(ez, 2);

    // Move the current explicit pressure term from b back into the linear
    // system. The FV Gauss gradient is linear in the pressure field, so its
    // coefficients can be assembled face-by-face.
    for (std::size_t c = 0; c < nc; ++c) {
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                // Same owner/neighbour rule as the continuity block below:
                // if c is the neighbour, nr points back to c itself.
                const int other_cell = owner
                    ? nr
                    : static_cast<int>(mesh.ownership().owner(f));
                if (other_cell < 0 ||
                    static_cast<std::size_t>(other_cell) >= nc ||
                    static_cast<std::size_t>(other_cell) == c)
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face topology on face " +
                        std::to_string(f) + " cell " + std::to_string(c));
                const std::size_t ncell = static_cast<std::size_t>(other_cell);
                // Match gauss_gradient_with_boundary() exactly: the
                // cell-centred pressure gradient uses a distance-weighted
                // face value, not an unconditional arithmetic average.
                // On a skewed mesh dc and dn need not be equal. Using 0.5 here
                // therefore assembled a different pressure-gradient operator
                // from the one used by the segregated momentum equation and
                // from the independent final momentum audit.
                const double dc =
                    (geometry.face_centres[f] - geometry.cell_centres[c]).mag();
                const double dn =
                    (geometry.cell_centres[ncell] - geometry.face_centres[f]).mag();
                if (!(dc > 0.0) || !(dn > 0.0) || !std::isfinite(dc) ||
                    !std::isfinite(dn))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face-centre distances on face " +
                        std::to_string(f));
                const double weight_c = dn / (dc + dn);
                const double weight_n = dc / (dc + dn);
                A.push_back(c, 3*nc + c, Sf.x * weight_c);
                A.push_back(c, 3*nc + ncell, Sf.x * weight_n);
                A.push_back(nc + c, 3*nc + c, Sf.y * weight_c);
                A.push_back(nc + c, 3*nc + ncell, Sf.y * weight_n);
                A.push_back(2*nc + c, 3*nc + c, Sf.z * weight_c);
                A.push_back(2*nc + c, 3*nc + ncell, Sf.z * weight_n);
            } else {
                const std::size_t patch = geometry.face_patch[f];
                const ScalarBoundaryCondition* bc = nullptr;
                ScalarBoundaryCondition local_bc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto it = pressure_bcs.find(name);
                    if (it != pressure_bcs.end()) {
                        local_bc = it->second;
                        bc = &local_bc;
                    }
                }
                if (bc && bc->type == ScalarBoundaryType::FIXED_VALUE) {
                    // Known boundary pressure is moved to the momentum RHS.
                    b(c) -= Sf.x * bc->value;
                    b(nc + c) -= Sf.y * bc->value;
                    b(2*nc + c) -= Sf.z * bc->value;
                } else {
                    // The transport-layer Gauss gradient treats zero-gradient
                    // pressure as owner-cell extrapolation. Retain that exact
                    // operator here so the coupled and segregated residuals
                    // remain comparable.
                    A.push_back(c, 3*nc + c, Sf.x);
                    A.push_back(nc + c, 3*nc + c, Sf.y);
                    A.push_back(2*nc + c, 3*nc + c, Sf.z);
                }
            }
        }
    }

    // Continuity is assembled from the same face velocity interpolation used
    // by the pressure-based solver. The pressure Schur diagonal below is the
    // exact diagonal of the same momentum-weighted (Rhie-Chow) face-flux
    // derivative and is passed explicitly to the coupled preconditioner.
    std::vector<double> coupled_pressure_schur_diagonal(nc, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        const std::size_t row = nv + c;
        const Offset off = mesh.cells().offsets_data()[c];
        const Offset count = mesh.cells().offsets_data()[c + 1] - off;
        for (Offset k = 0; k < count; ++k) {
            const std::size_t f = mesh.cells().faces_data()[off + k];
            const bool owner = mesh.ownership().owner(f) == c;
            const double sign = owner ? 1.0 : -1.0;
            const Vec3 Sf = geometry.face_area_vectors[f] * sign;
            const auto nr = mesh.ownership().neighbour(f);

            if (nr >= 0) {
                // neighbour(f) is the opposite cell only on the owner side.
                // When c is the neighbour, neighbour(f) == c; using it directly
                // creates a self-coupled face (zero centre distance) and corrupts
                // the coupled pressure Schur coefficient.
                const int other_cell = owner
                    ? nr
                    : static_cast<int>(mesh.ownership().owner(f));
                if (other_cell < 0 ||
                    static_cast<std::size_t>(other_cell) >= nc ||
                    static_cast<std::size_t>(other_cell) == c)
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face topology on face " +
                        std::to_string(f) + " cell " + std::to_string(c));
                const std::size_t ncell = static_cast<std::size_t>(other_cell);
                const double half_rho = 0.5 * rho;
                A.push_back(row, c, half_rho * Sf.x);
                A.push_back(row, ncell, half_rho * Sf.x);
                A.push_back(row, nc + c, half_rho * Sf.y);
                A.push_back(row, nc + ncell, half_rho * Sf.y);
                A.push_back(row, 2*nc + c, half_rho * Sf.z);
                A.push_back(row, 2*nc + ncell, half_rho * Sf.z);

                const Vec3 rawSf = geometry.face_area_vectors[f];
                const double area = rawSf.mag();
                const double d = (geometry.cell_centres[ncell] -
                                  geometry.cell_centres[c]).mag();
                const double nx = rawSf.x / area;
                const double ny = rawSf.y / area;
                const double nz = rawSf.z / area;
                const double aox = ex.diagonal[c], anx = ex.diagonal[ncell];
                const double aoy = ey.diagonal[c], any = ey.diagonal[ncell];
                const double aoz = ez.diagonal[c], anz = ez.diagonal[ncell];
                if (!(aox > 0.0) || !(anx > 0.0) ||
                    !(aoy > 0.0) || !(any > 0.0) ||
                    !(aoz > 0.0) || !(anz > 0.0) ||
                    !std::isfinite(aox) || !std::isfinite(anx) ||
                    !std::isfinite(aoy) || !std::isfinite(any) ||
                    !std::isfinite(aoz) || !std::isfinite(anz))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid momentum diagonal on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                // The coupled Schur complement uses the same dAU = V/A_P
                // pressure response as the segregated Rhie-Chow operator.
                const double rAUx_f = 0.5 * (
                    geometry.cell_volumes[c] / aox +
                    geometry.cell_volumes[ncell] / anx);
                const double rAUy_f = 0.5 * (
                    geometry.cell_volumes[c] / aoy +
                    geometry.cell_volumes[ncell] / any);
                const double rAUz_f = 0.5 * (
                    geometry.cell_volumes[c] / aoz +
                    geometry.cell_volumes[ncell] / anz);
                const double rfn = rAUx_f*nx*nx + rAUy_f*ny*ny + rAUz_f*nz*nz;
                if (!(rfn > 0.0) || !std::isfinite(rfn))
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid face rAU on face " +
                        std::to_string(f));
                const double D = rho * rfn * area / d;
                if (!(d > 0.0) || !(area > 0.0) ||
                    !std::isfinite(d) || !std::isfinite(area) ||
                    !(D > 0.0) || !std::isfinite(D)) {
                    if (diagnostics.coupled_matrix_summary) {
                        std::cerr << "COUPLED_FACE_DEBUG face=" << f
                                  << " cells=" << c << "/" << ncell
                                  << " area=" << area
                                  << " distance=" << d
                                  << " a=(" << aox << "," << anx
                                  << ";" << aoy << "," << any
                                  << ";" << aoz << "," << anz << ")"
                                  << " rAU=(" << rAUx_f << "," << rAUy_f
                                  << "," << rAUz_f << ")"
                                  << " rfn=" << rfn
                                  << " D=" << D << "\n";
                    }
                    throw std::runtime_error(
                        "solve_coupled_momentum_continuity: invalid pressure coefficient on face " +
                        std::to_string(f) + " cells " + std::to_string(c) + "/" +
                        std::to_string(ncell));
                }
                // phi_p = D (p_P - p_N). Only this orthogonal derivative is
                // implicit in the coupled matrix. The non-orthogonal remainder
                // is a deferred correction evaluated from p_old below, exactly
                // as in the segregated Rhie-Chow operator.
                A.push_back(row, nv + c, D);
                A.push_back(row, nv + ncell, -D);
                coupled_pressure_schur_diagonal[c] += D;

                const Vec3 dvec = geometry.cell_centres[ncell] -
                                  geometry.cell_centres[c];
                const Vec3 e{dvec.x/d, dvec.y/d, dvec.z/d};
                const double orthogonal_area = rawSf.dot(e);
                const Vec3 Snon{
                    rawSf.x - orthogonal_area*e.x,
                    rawSf.y - orthogonal_area*e.y,
                    rawSf.z - orthogonal_area*e.z};
                const Vec3 gp{
                    0.5*(grad_p_old.component_data(0)[c] +
                         grad_p_old.component_data(0)[ncell]),
                    0.5*(grad_p_old.component_data(1)[c] +
                         grad_p_old.component_data(1)[ncell]),
                    0.5*(grad_p_old.component_data(2)[c] +
                         grad_p_old.component_data(2)[ncell])};
                const double nonorth_flux = rho * rfn * gp.dot(Snon);
                // The continuity equation is div(phi)=0 and the pressure
                // correction is written as phi = phi_U - D grad(p).
                // Therefore the explicit non-orthogonal pressure flux is
                // subtracted from the current continuity residual.
                // nonorth_flux is evaluated with the raw owner->neighbour
                // face orientation, so the neighbour-side visit must flip
                // its sign: like the implicit pressure part, which enters
                // every row as D*(p_self - p_other), the deferred flux
                // must enter this cell's outward divergence with the
                // face sign.
                b(row) += sign * nonorth_flux;
            } else {
                const std::size_t patch = geometry.face_patch[f];
                if (patch < mesh.boundary().n_patches() &&
                    mesh.boundary().patch(patch).type == PatchType::EMPTY)
                    continue;
                const VelocityBoundaryCondition* vbc = nullptr;
                VelocityBoundaryCondition local_vbc;
                const ScalarBoundaryCondition* pbc = nullptr;
                ScalarBoundaryCondition local_pbc;
                if (patch < mesh.boundary().n_patches()) {
                    const auto& name = mesh.boundary().patch(patch).name;
                    const auto vit = velocity_bcs.find(name);
                    if (vit != velocity_bcs.end()) {
                        local_vbc = vit->second;
                        vbc = &local_vbc;
                    }
                    const auto pit = pressure_bcs.find(name);
                    if (pit != pressure_bcs.end()) {
                        local_pbc = pit->second;
                        pbc = &local_pbc;
                    }
                }

                if (vbc && vbc->type == VelocityBoundaryCondition::Type::FIXED_VALUE) {
                    b(row) -= rho * vbc->value.dot(Sf);
                } else {
                    A.push_back(row, c, rho * Sf.x);
                    A.push_back(row, nc + c, rho * Sf.y);
                    A.push_back(row, 2*nc + c, rho * Sf.z);
                }

                if (pbc && pbc->type == ScalarBoundaryType::FIXED_VALUE) {
                    const Vec3 rawSf = geometry.face_area_vectors[f];
                    const double area = rawSf.mag();
                    const double d = (geometry.face_centres[f] -
                                      geometry.cell_centres[c]).mag();
                    const double nx = rawSf.x / area;
                    const double ny = rawSf.y / area;
                    const double nz = rawSf.z / area;
                    const double a_x = ex.diagonal[c];
                    const double a_y = ey.diagonal[c];
                    const double a_z = ez.diagonal[c];
                    if (!(a_x > 0.0) || !(a_y > 0.0) || !(a_z > 0.0) ||
                        !std::isfinite(a_x) || !std::isfinite(a_y) || !std::isfinite(a_z))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary momentum diagonal on face " +
                            std::to_string(f) + " cell " + std::to_string(c));
                    const double rAUx = geometry.cell_volumes[c] / a_x;
                    const double rAUy = geometry.cell_volumes[c] / a_y;
                    const double rAUz = geometry.cell_volumes[c] / a_z;
                    const double rfn = rAUx*nx*nx + rAUy*ny*ny + rAUz*nz*nz;
                    if (!(rfn > 0.0) || !std::isfinite(rfn))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary face rAU on face " +
                            std::to_string(f));
                    const double D = rho * rfn * area / d;
                    if (!(D > 0.0) || !std::isfinite(D))
                        throw std::runtime_error(
                            "solve_coupled_momentum_continuity: invalid boundary pressure coefficient on face " +
                            std::to_string(f));
                    A.push_back(row, nv + c, D);
                    coupled_pressure_schur_diagonal[c] += D;
                    b(row) += D * pbc->value;
                }
            }
        }
    }

    // A pure-Neumann pressure field has one gauge null mode. Eliminate the
    // reference pressure column and replace the reference continuity row by
    // p_ref = reference_value. This keeps the coupled system deterministic.
    bool has_fixed_pressure = false;
    for (const auto& [name, bc] : pressure_bcs) {
        (void)name;
        has_fixed_pressure = has_fixed_pressure ||
            bc.type == ScalarBoundaryType::FIXED_VALUE;
    }

    if (!has_fixed_pressure) {
        // Finalize the assembled COO matrix before reading its CSR arrays.
        // Without this explicit finalize(), row_offsets_/values_/columns_ still
        // describe an empty CSR shell, so gauge elimination would silently
        // discard every assembled momentum and continuity coefficient.
        A.finalize();

        // Rebuild the matrix while eliminating the known gauge column.
        SparseMatrix reduced(n, n);
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const auto* va = A.values_data();
        for (std::size_t r = 0; r < n; ++r) {
            if (r == nv + reference_cell) {
                reduced.push_back(r, nv + reference_cell, 1.0);
                b(r) = reference_value;
                continue;
            }
            for (std::uint32_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t col = co[k];
                if (col == nv + reference_cell) {
                    b(r) -= va[k] * reference_value;
                } else {
                    reduced.push_back(r, col, va[k]);
                }
            }
        }
        reduced.finalize();
        A = std::move(reduced);
    } else {
        A.finalize();
    }

    Vector x(n, 0.0);
    for (std::size_t c = 0; c < nc; ++c) {
        x(c) = U_old.component_data(0)[c];
        x(nc + c) = U_old.component_data(1)[c];
        x(2*nc + c) = U_old.component_data(2)[c];
        x(nv + c) = p_old(c);
    }

    // PR #832 fix: zero the gauge row of the initial guess before GMRES.
    // The gauge row is an identity row with rhs = reference_value (default 0.0).
    // If x0(gauge) = p_old(reference_cell) is non-zero, the right-preconditioned
    // GMRES may struggle to correct it because the preconditioner preserves
    // the gauge component. Setting x0(gauge) = 0.0 gives the solver the correct
    // starting point and avoids the stagnation observed in PR #830.
    x(nv + reference_cell) = controls.pressure_reference_value;

    // Reject non-finite coefficients before entering GMRES. A bad local
    // coefficient must be diagnosed at assembly time rather than surfacing as
    // an opaque Krylov failure.
    {
        const auto* values = A.values_data();
        for (std::size_t k = 0; k < A.nnz(); ++k) {
            if (!std::isfinite(values[k]))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite matrix coefficient at nnz " +
                    std::to_string(k));
        }
        for (std::size_t i = 0; i < n; ++i) {
            if (!std::isfinite(b(i)))
                throw std::runtime_error(
                    "solve_coupled_momentum_continuity: non-finite RHS at row " +
                    std::to_string(i));
        }
    }

    // Optional algebraic microscope. It is deliberately computed from
    // the exact post-gauge matrix passed to GMRES, so it exposes the actual
    // M/G/D/C blocks rather than reconstructed proxy operators.
    if (diagnostics.coupled_matrix_summary) {
        double m2 = 0.0, g2 = 0.0, d2 = 0.0, c2 = 0.0;
        std::size_t zero_rows = 0;
        std::vector<std::size_t> row_nnz(n, 0);
        std::vector<std::size_t> col_nnz(n, 0);
        const auto* ro = A.row_offsets_data();
        const auto* co = A.columns_data();
        const auto* va = A.values_data();

        for (std::size_t r = 0; r < n; ++r) {
            row_nnz[r] = ro[r + 1] - ro[r];
            if (row_nnz[r] == 0) ++zero_rows;
            for (std::size_t k = ro[r]; k < ro[r + 1]; ++k) {
                const std::size_t col = co[k];
                const double value = va[k];
                col_nnz[col]++;
                if (r < nv && col < nv)
                    m2 += value * value;
                else if (r < nv && col >= nv)
                    g2 += value * value;
                else if (r >= nv && col < nv)
                    d2 += value * value;
                else
                    c2 += value * value;
            }
        }
        std::size_t zero_cols = 0;
        for (const auto count : col_nnz)
            if (count == 0) ++zero_cols;

        double schur_diag_min = std::numeric_limits<double>::infinity();
        double schur_diag_max = 0.0;
        double schur_mismatch_max = 0.0;
        std::size_t schur_near_zero = 0;
        for (std::size_t c = 0; c < nc; ++c) {
            const std::size_t pr = nv + c;
            double effective = 0.0;
            for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k) {
                const std::size_t u = co[k];
                if (u >= nv) continue;
                const double d = va[k];
                double inv_m = 0.0;
                bool found = false;
                for (std::size_t mk = ro[u]; mk < ro[u + 1]; ++mk) {
                    if (co[mk] == u) {
                        const double diag = va[mk];
                        if (std::isfinite(diag) && std::abs(diag) > 1e-30) {
                            inv_m = 1.0 / diag;
                            found = true;
                        }
                        break;
                    }
                }
                if (!found) continue;
                for (std::size_t gk = ro[u]; gk < ro[u + 1]; ++gk) {
                    const std::size_t q = co[gk];
                    if (q == pr)
                        effective -= d * inv_m * va[gk];
                }
            }
            for (std::size_t k = ro[pr]; k < ro[pr + 1]; ++k)
                if (co[k] == pr) effective += va[k];

            schur_diag_min = std::min(schur_diag_min, effective);
            schur_diag_max = std::max(schur_diag_max, std::abs(effective));
            if (std::abs(effective) <=
                1e-12 * std::max(1.0, schur_diag_max))
                ++schur_near_zero;
            if (c < coupled_pressure_schur_diagonal.size() &&
                std::isfinite(coupled_pressure_schur_diagonal[c]))
                schur_mismatch_max = std::max(
                    schur_mismatch_max,
                    std::abs(effective - coupled_pressure_schur_diagonal[c]));
        }

        std::cerr << "COUPLED_MATRIX_SUMMARY n=" << n
                  << " nnz=" << A.nnz()
                  << " M_l2=" << std::sqrt(m2)
                  << " G_l2=" << std::sqrt(g2)
                  << " D_l2=" << std::sqrt(d2)
                  << " C_l2=" << std::sqrt(c2)
                  << " zero_rows=" << zero_rows
                  << " zero_cols=" << zero_cols
                  << " schur_diag_min=" << schur_diag_min
                  << " schur_diag_absmax=" << schur_diag_max
                  << " schur_near_zero=" << schur_near_zero
                  << " rhie_chow_schur_vs_algebraic_diag_delta=" << schur_mismatch_max
                  << " gauge_row=" << (nv + reference_cell)
                  << "\n";
    }
    // Production continuation (steady_incompressible_solver.h, right after
    // the COUPLED_MATRIX_SUMMARY block): the gauge entry of the Rhie-Chow
    // Schur diagonal hint is set to one before the linear plan is selected.
    // The hint is only consumed by the summary's mismatch metric in
    // production (it is not passed to the preconditioner); it is returned
    // and dumped here for completeness, with the same override applied.
    coupled_pressure_schur_diagonal[reference_cell] = 1.0;

    return CoupledAuditSystem{std::move(A), std::move(b), std::move(x),
                              std::move(coupled_pressure_schur_diagonal)};
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::string out_dir = argc > 1 ? argv[1] : ".";
        const auto t0 = std::chrono::steady_clock::now();
        const auto elapsed = [&t0]() {
            return std::chrono::duration<double>(
                std::chrono::steady_clock::now() - t0).count();
        };

        constexpr std::size_t nx = 32, ny = 32;
        const double reynolds = 100.0;
        std::cerr << "AUDIT_START case=ghia_pcd_32x32 reynolds=" << reynolds
                  << " out_dir=" << out_dir << "\n";

        // ---- Phase 0: mesh, boundary conditions, controls (--pcd, 32x32). ----
        Mesh mesh = make_cavity_mesh(nx, ny);
        const FvGeometry geometry = build_fv_geometry(mesh);
        const std::size_t nc = mesh.n_cells();

        VelocityBoundaryConditions ubc;
        ubc["bottom"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ubc["left"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ubc["right"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ubc["front"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ubc["back"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{0,0,0}};
        ubc["top"]={VelocityBoundaryCondition::Type::FIXED_VALUE,{1,0,0}};

        ScalarBoundaryConditions pbc;
        for (const char* name : {"bottom","top","left","right","front","back"})
            pbc[name] = {ScalarBoundaryType::ZERO_GRADIENT, 0.0, 0.0};

        // Controls byte-for-byte from tests/validation/test_ghia_cavity.cpp
        // (--pcd path), except convergence.max_iterations which is replaced
        // by an iteration_output_callback stop after iteration 1 (the audit
        // needs exactly the post-iteration-1 state; the coupled solve of
        // iteration 1 itself is untouched by this change).
        IncompressibleSolverControls controls;
        controls.algorithm = PressureVelocityAlgorithm::COUPLED;
        controls.density = 1.0;
        controls.kinematic_viscosity = 1.0 / reynolds;
        controls.linear_max_iterations = 10000;
        controls.linear_tolerance = 1e-8;
        controls.pressure_reference_cell = (ny / 2) * nx + (nx / 2);
        controls.pressure_reference_value = 0.0;
        controls.coupling.alpha_u = 0.7;
        controls.coupling.alpha_p = 0.3;
        controls.coupling.coupled_max_iterations = 5000;
        controls.coupling.coupled_linear_tolerance = 1e-10;
        controls.coupled_linear_solver.krylov = KrylovModel::FGMRES;
        controls.coupled_linear_solver.preconditioner = PreconditionerModel::PCD;
        // #898/#894 restart rule (test comment, verbatim policy).
        const std::size_t coupled_unknowns = 4 * nx * ny;
        const std::size_t coupled_budget = controls.coupling.coupled_max_iterations;
        controls.coupled_linear_solver.gmres_restart = static_cast<int>(
            coupled_unknowns <= 4096
                ? std::min<std::size_t>(coupled_unknowns,
                                        std::max<std::size_t>(1024, coupled_budget))
                : std::min<std::size_t>(coupled_unknowns, 512));
        // #908/#870 nested Fp restart rule (test comment, verbatim policy).
        const std::size_t pressure_unknowns = nx * ny;
        const std::size_t pcd_fp_budget = 200;
        controls.pcd_convection_diffusion_max_iterations = pcd_fp_budget;
        controls.pcd_convection_diffusion_linear_solver.gmres_restart =
            static_cast<int>(std::min<std::size_t>(
                pressure_unknowns, std::max<std::size_t>(512, pcd_fp_budget)));
        // schur_model: the CouplingControls default is already PCD; set it
        // explicitly so the audit never depends on that default silently.
        controls.coupling.schur_model = CoupledSchurModel::PCD;
        controls.use_bounded_convection = true;
        controls.convection_scheme = ConvectionScheme::UPWIND;
        controls.diagnostics.coupled_matrix_summary = true;
        controls.diagnostics.iteration_trace = true;
        controls.diagnostics.iteration_trace_frequency = 50;
        controls.convergence.max_iterations = 2;
        controls.convergence.relative_tolerance = 1e-8;
        controls.convergence.continuity_tolerance = 1e-8;

        const std::size_t reference_cell = controls.pressure_reference_cell;
        const double reference_value = controls.pressure_reference_value;
        const double mu_eff = controls.density *
            (controls.kinematic_viscosity + controls.turbulent_viscosity);

        // ---- Phase 1: iteration-1 momentum equations from the zero state. ----
        // Exact replica of the production first assembly: mass flux from
        // U=0, pressure gradient from p=0, zero body force, convected fields
        // from U=0, upwind/bounded — deterministic inputs, so this is
        // bit-identical to the equations production assembles before its
        // first coupled solve.
        Field<double, Location::CELL> U0(nc, "U0", "m/s", 3); U0.fill(0.0);
        Field<double, Location::CELL> p0(nc, "p0", "Pa", 1); p0.fill(0.0);

        const auto mass_flux_0 = make_mass_flux(mesh, geometry, U0, controls.density, ubc);
        const auto grad_p0 = gauss_gradient_with_boundary(p0, mesh, geometry, pbc);

        Field<double, Location::CELL> body_x(nc, "body_x", "N/m3", 1);
        Field<double, Location::CELL> body_y(nc, "body_y", "N/m3", 1);
        Field<double, Location::CELL> body_z(nc, "body_z", "N/m3", 1);
        for (std::size_t c = 0; c < nc; ++c) {
            body_x(c) = controls.density * controls.body_force.x;
            body_y(c) = controls.density * controls.body_force.y;
            body_z(c) = controls.density * controls.body_force.z;
        }

        Field<double, Location::CELL> u_x0(nc, "Ux0", "m/s", 1);
        Field<double, Location::CELL> u_y0(nc, "Uy0", "m/s", 1);
        Field<double, Location::CELL> u_z0(nc, "Uz0", "m/s", 1);
        for (std::size_t c = 0; c < nc; ++c) {
            u_x0(c) = U0.component_data(0)[c];
            u_y0(c) = U0.component_data(1)[c];
            u_z0(c) = U0.component_data(2)[c];
        }
        ScalarBoundaryConditions ubc_x, ubc_y, ubc_z;
        for (const auto& [name, bc] : ubc) {
            ubc_x[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.x, 0.0};
            ubc_y[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.y, 0.0};
            ubc_z[name] = {bc.type == VelocityBoundaryCondition::Type::FIXED_VALUE
                               ? ScalarBoundaryType::FIXED_VALUE : ScalarBoundaryType::ZERO_GRADIENT,
                           bc.value.z, 0.0};
        }

        const auto ex1 = assemble_momentum_component(
            mesh, geometry, mass_flux_0, grad_p0, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme, &u_x0);
        const auto ey1 = assemble_momentum_component(
            mesh, geometry, mass_flux_0, grad_p0, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme, &u_y0);
        const auto ez1 = assemble_momentum_component(
            mesh, geometry, mass_flux_0, grad_p0, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme, &u_z0);
        std::cerr << "AUDIT_PHASE1 iteration1_momentum_equations_assembled t_s="
                  << elapsed() << "\n";

        // ---- Phase 1b: audit assembly of the ITERATION-1 coupled system. ----
        // Expected values: production run 38809298638,
        // diag-evidence/logs/test_ghia_cavity_n8_pcd.log, first
        // COUPLED_MATRIX_SUMMARY line.
        const auto audit1 = audit_assemble_coupled(
            mesh, geometry, ex1, ey1, ez1, U0, p0, ubc, pbc, controls.density,
            reference_cell, reference_value, controls, controls.diagnostics);

        // ---- Phase 2: exactly one production nonlinear iteration (~20 s). ----
        // The real solve_steady_incompressible runs with the exact --pcd
        // controls; the callback stops the loop after iteration 1. U and p
        // then hold the relaxed, gauged post-iteration-1 state (the state
        // production uses as U_old/p_old for iteration 2), and the production
        // code path has printed its own COUPLED_MATRIX_SUMMARY for
        // iteration 1 into the log.
        Field<double, Location::CELL> U = U0;
        Field<double, Location::CELL> p = p0;
        controls.iteration_output_callback =
            [](std::size_t iter, double, const Mesh&,
               const Field<double, Location::CELL>&,
               const Field<double, Location::CELL>&) -> bool {
                std::cerr << "AUDIT_CALLBACK iteration=" << iter
                          << " stop_after_iteration=1\n";
                return false;
            };
        const auto solve = solve_steady_incompressible(mesh, U, p, ubc, pbc, controls);
        std::cerr << "AUDIT_PHASE2 production_nonlinear_iterations=" << solve.iterations
                  << " t_s=" << elapsed() << "\n";
        if (solve.iterations != 1)
            throw std::runtime_error(
                "audit: production solve did not stop after exactly one iteration");
        if (!solve.history.empty()) {
            const auto& it = solve.history.back();
            std::fprintf(stderr,
                "AUDIT_PHASE2_STATE continuity_normalized=%.6g velocity_change_inf=%.6g "
                "pressure_change_inf=%.6g effective_alpha_u=%.6g effective_alpha_p=%.6g\n",
                it.continuity_normalized, it.velocity_change_inf,
                it.pressure_change_inf, it.effective_alpha_u, it.effective_alpha_p);
        }

        // ---- Phase 3: iteration-2 assembly state (production loop replica). ----
        // coupled_rAU comes from the iteration-1 momentum diagonals (the
        // equations production assembled at the start of iteration 1), the
        // flux is the authoritative Rhie-Chow flux production builds after
        // the coupled update + Picard relaxation + gauge policy.
        std::array<std::vector<double>, 3> coupled_rAU;
        for (std::size_t d = 0; d < 3; ++d) coupled_rAU[d].resize(nc);
        for (std::size_t c = 0; c < nc; ++c) {
            if (!(ex1.diagonal[c] > 0.0) || !std::isfinite(ex1.diagonal[c]) ||
                !(ey1.diagonal[c] > 0.0) || !std::isfinite(ey1.diagonal[c]) ||
                !(ez1.diagonal[c] > 0.0) || !std::isfinite(ez1.diagonal[c]))
                throw std::runtime_error("audit: invalid iteration-1 momentum diagonal");
            coupled_rAU[0][c] = 1.0 / ex1.diagonal[c];
            coupled_rAU[1][c] = 1.0 / ey1.diagonal[c];
            coupled_rAU[2][c] = 1.0 / ez1.diagonal[c];
        }
        const auto mass_flux_1 = make_rhie_chow_mass_flux(
            mesh, geometry, U, p, coupled_rAU, controls.density, ubc, pbc);
        const auto grad_p1 = gauss_gradient_with_boundary(p, mesh, geometry, pbc);

        Field<double, Location::CELL> u_x1(nc, "Ux1", "m/s", 1);
        Field<double, Location::CELL> u_y1(nc, "Uy1", "m/s", 1);
        Field<double, Location::CELL> u_z1(nc, "Uz1", "m/s", 1);
        for (std::size_t c = 0; c < nc; ++c) {
            u_x1(c) = U.component_data(0)[c];
            u_y1(c) = U.component_data(1)[c];
            u_z1(c) = U.component_data(2)[c];
        }
        const auto ex2 = assemble_momentum_component(
            mesh, geometry, mass_flux_1, grad_p1, body_x, mu_eff, ubc_x, 0,
            controls.use_bounded_convection, controls.convection_scheme, &u_x1);
        const auto ey2 = assemble_momentum_component(
            mesh, geometry, mass_flux_1, grad_p1, body_y, mu_eff, ubc_y, 1,
            controls.use_bounded_convection, controls.convection_scheme, &u_y1);
        const auto ez2 = assemble_momentum_component(
            mesh, geometry, mass_flux_1, grad_p1, body_z, mu_eff, ubc_z, 2,
            controls.use_bounded_convection, controls.convection_scheme, &u_z1);
        std::cerr << "AUDIT_PHASE3 iteration2_momentum_equations_assembled t_s="
                  << elapsed() << "\n";

        // ---- Phase 4: audit assembly of the ITERATION-2 coupled system. ----
        // Expected values: production run 38809298638, second
        // COUPLED_MATRIX_SUMMARY line (the audit target state).
        const auto audit2 = audit_assemble_coupled(
            mesh, geometry, ex2, ey2, ez2, U, p, ubc, pbc, controls.density,
            reference_cell, reference_value, controls, controls.diagnostics);

        // ---- Phase 5: PCD pressure operators from the same state. ----
        // Production passes U_old (= U1 here) to assemble_pcd_pressure_operators.
        const auto pcd_ops = assemble_pcd_pressure_operators(
            mesh, geometry, U, controls.density, controls.kinematic_viscosity,
            ubc, pbc, reference_cell, ConvectionScheme::UPWIND);
        std::cerr << "AUDIT_PHASE5 pcd_operators_assembled"
                  << " Mp_nnz=" << pcd_ops.mass.nnz()
                  << " Kp_nnz=" << pcd_ops.laplacian.nnz()
                  << " Fp_nnz=" << pcd_ops.convection_diffusion.nnz()
                  << " t_s=" << elapsed() << "\n";

        // ---- Validation against the production CI log (run 38809298638). ----
        const std::size_t nv = 3 * nc;
        bool all_ok = true;
        const auto validate = [&](const char* label, const SparseMatrix& A,
                                  double eM, double eG, double eD, double eC,
                                  double eSchur, std::size_t eNnz) {
            const auto norms = audit_block_norms(A, nv);
            const double smin = audit_schur_effective_diag_min(A, nv);
            const bool okM = audit_close(norms.m_l2, eM);
            const bool okG = audit_close(norms.g_l2, eG);
            const bool okD = audit_close(norms.d_l2, eD);
            const bool okC = audit_close(norms.c_l2, eC);
            const bool okS = audit_close(smin, eSchur);
            const bool okN = norms.nnz == eNnz;
            std::fprintf(stderr,
                "AUDIT_%s_SUMMARY n=%llu nnz=%llu M_l2=%.6g G_l2=%.6g D_l2=%.6g C_l2=%.6g "
                "schur_diag_min=%.6g gauge_row=%llu\n",
                label,
                static_cast<unsigned long long>(A.n_rows()),
                static_cast<unsigned long long>(norms.nnz),
                norms.m_l2, norms.g_l2, norms.d_l2, norms.c_l2, smin,
                static_cast<unsigned long long>(nv + reference_cell));
            std::fprintf(stderr,
                "AUDIT_%s_EXPECTED nnz=%llu M_l2=%.6g G_l2=%.6g D_l2=%.6g C_l2=%.6g "
                "schur_diag_min=%.6g\n",
                label, static_cast<unsigned long long>(eNnz),
                eM, eG, eD, eC, eSchur);
            std::fprintf(stderr,
                "AUDIT_%s_MATCH M_l2=%s G_l2=%s D_l2=%s C_l2=%s schur_diag_min=%s nnz=%s\n",
                label, okM ? "ok" : "MISMATCH", okG ? "ok" : "MISMATCH",
                okD ? "ok" : "MISMATCH", okC ? "ok" : "MISMATCH",
                okS ? "ok" : "MISMATCH", okN ? "ok" : "MISMATCH");
            const bool ok = okM && okG && okD && okC && okS && okN;
            std::fprintf(stderr, "AUDIT_%s_MATCH overall=%s\n",
                         label, ok ? "PASS" : "FAIL");
            all_ok = all_ok && ok;
        };
        validate("ITER1", audit1.A, 2.54234, 1.43614, 1.39124, 2.35617,
                 0.0536722, 76991);
        validate("ITER2", audit2.A, 2.85644, 1.43614, 1.39124, 2.1722,
                 0.0465948, 76991);
        // Production iteration-1 nonlinear state (INCOMPRESSIBLE_ITER line).
        if (!solve.history.empty()) {
            const bool ok_cont = audit_close(solve.history.back().continuity_normalized,
                                             0.00310014);
            std::fprintf(stderr, "AUDIT_ITER1_NONLINEAR continuity_normalized=%.6g "
                                 "expected=0.00310014 match=%s\n",
                         solve.history.back().continuity_normalized,
                         ok_cont ? "ok" : "MISMATCH");
            all_ok = all_ok && ok_cont;
        }

        // ---- Phase 6: bit-exact dumps + checksums + meta. ----
        const auto to_vector = [](const Vector& v) {
            std::vector<double> out(v.size());
            for (std::size_t i = 0; i < v.size(); ++i) out[i] = v(i);
            return out;
        };
        const auto b1 = to_vector(audit1.b);
        const auto x01 = to_vector(audit1.x0);
        const auto b2 = to_vector(audit2.b);
        const auto x02 = to_vector(audit2.x0);
        const auto schur2 = audit2.schur_diagonal;
        const auto schur1 = audit1.schur_diagonal;

        dump_csr_chunks(out_dir + "/A_iter1.csr", "A_iter1", audit1.A, 1024);
        dump_csr_chunks(out_dir + "/A_iter2.csr", "A_iter2", audit2.A, 1024);
        dump_csr_chunks(out_dir + "/Mp.csr", "Mp", pcd_ops.mass, 1024);
        dump_csr_chunks(out_dir + "/Kp.csr", "Kp", pcd_ops.laplacian, 1024);
        dump_csr_chunks(out_dir + "/Fp.csr", "Fp", pcd_ops.convection_diffusion, 1024);
        dump_vector_file(out_dir + "/b_iter1.vec", "b_iter1", b1.data(), b1.size());
        dump_vector_file(out_dir + "/x0_iter1.vec", "x0_iter1", x01.data(), x01.size());
        dump_vector_file(out_dir + "/b_iter2.vec", "b_iter2", b2.data(), b2.size());
        dump_vector_file(out_dir + "/x0_iter2.vec", "x0_iter2", x02.data(), x02.size());
        dump_vector_file(out_dir + "/schur_diag_iter1.vec", "schur_diag_iter1",
                         schur1.data(), schur1.size());
        dump_vector_file(out_dir + "/schur_diag_iter2.vec", "schur_diag_iter2",
                         schur2.data(), schur2.size());
        dump_field_file(out_dir + "/U_iter2.field", "U_iter2", U);
        dump_field_file(out_dir + "/p_iter2.field", "p_iter2", p);
        dump_field_file(out_dir + "/mass_flux_iter2.field", "mass_flux_iter2",
                        mass_flux_1);

        const auto checksum = [](const char* name, const double* v, std::size_t n) {
            std::uint64_t h = fnv1a_init();
            for (std::size_t i = 0; i < n; ++i) {
                std::uint64_t bits = 0;
                std::memcpy(&bits, &v[i], sizeof(bits));
                h = fnv1a_bytes(h, &bits, sizeof(bits));
            }
            std::fprintf(stderr, "AUDIT_CHECKSUM name=%s count=%llu fnv1a=%016llx\n",
                         name, static_cast<unsigned long long>(n),
                         static_cast<unsigned long long>(h));
            return h;
        };
        const auto matrix_checksum = [](const char* name, const SparseMatrix& m) {
            std::uint64_t h = fnv1a_init();
            const auto* va = m.values_data();
            for (std::size_t k = 0; k < m.nnz(); ++k) {
                std::uint64_t bits = 0;
                std::memcpy(&bits, &va[k], sizeof(bits));
                h = fnv1a_bytes(h, &bits, sizeof(bits));
            }
            std::fprintf(stderr, "AUDIT_CHECKSUM name=%s count=%llu fnv1a=%016llx\n",
                         name, static_cast<unsigned long long>(m.nnz()),
                         static_cast<unsigned long long>(h));
            return h;
        };
        std::vector<std::string> meta;
        meta.push_back("# Ghia 32x32 PCD coupled-matrix audit dump (#914)");
        meta.push_back("production_reference=run_38809298638_diag-evidence_logs_test_ghia_cavity_n8_pcd");
        meta.push_back("case=ghia_cavity_32x32_re100_pcd");
        meta.push_back("n_cells=" + std::to_string(nc));
        meta.push_back("n_faces=" + std::to_string(mesh.n_faces()));
        meta.push_back("n=" + std::to_string(4 * nc));
        meta.push_back("nnz=" + std::to_string(audit2.A.nnz()));
        meta.push_back("reference_cell=" + std::to_string(reference_cell));
        meta.push_back("gauge_row=" + std::to_string(3 * nc + reference_cell));
        meta.push_back("density=" + std::to_string(controls.density));
        meta.push_back("kinematic_viscosity=" + std::to_string(controls.kinematic_viscosity));
        meta.push_back("alpha_u=" + std::to_string(controls.coupling.alpha_u));
        meta.push_back("alpha_p=" + std::to_string(controls.coupling.alpha_p));
        meta.push_back("coupled_max_iterations=" + std::to_string(controls.coupling.coupled_max_iterations));
        meta.push_back("coupled_linear_tolerance=1e-10");
        meta.push_back("coupled_krylov=FGMRES");
        meta.push_back("coupled_restart=" + std::to_string(controls.coupled_linear_solver.gmres_restart));
        meta.push_back("schur_model=pcd");
        meta.push_back("pcd_fp_budget=" + std::to_string(controls.pcd_convection_diffusion_max_iterations));
        meta.push_back("pcd_fp_restart=" + std::to_string(controls.pcd_convection_diffusion_linear_solver.gmres_restart));
        meta.push_back("pcd_fp_tolerance=1e-8");
        meta.push_back("pcd_kp_krylov=GMRES");
        meta.push_back("pcd_kp_preconditioner=NativeAMG");
        meta.push_back("pcd_kp_budget=" + std::to_string(controls.pcd_laplacian_max_iterations));
        meta.push_back("pcd_kp_restart=" + std::to_string(controls.pcd_laplacian_linear_solver.gmres_restart));
        meta.push_back("pcd_kp_tolerance=1e-8");
        meta.push_back("Mp_nnz=" + std::to_string(pcd_ops.mass.nnz()));
        meta.push_back("Kp_nnz=" + std::to_string(pcd_ops.laplacian.nnz()));
        meta.push_back("Fp_nnz=" + std::to_string(pcd_ops.convection_diffusion.nnz()));
        const auto hex = [](std::uint64_t h) {
            char buf[32];
            std::snprintf(buf, sizeof(buf), "%016llx", static_cast<unsigned long long>(h));
            return std::string(buf);
        };
        meta.push_back("fnv1a_A_iter1_values=" + hex(matrix_checksum("A_iter1_values", audit1.A)));
        meta.push_back("fnv1a_A_iter2_values=" + hex(matrix_checksum("A_iter2_values", audit2.A)));
        meta.push_back("fnv1a_Mp_values=" + hex(matrix_checksum("Mp_values", pcd_ops.mass)));
        meta.push_back("fnv1a_Kp_values=" + hex(matrix_checksum("Kp_values", pcd_ops.laplacian)));
        meta.push_back("fnv1a_Fp_values=" + hex(matrix_checksum("Fp_values", pcd_ops.convection_diffusion)));
        meta.push_back("fnv1a_b_iter1=" + hex(checksum("b_iter1", b1.data(), b1.size())));
        meta.push_back("fnv1a_x0_iter1=" + hex(checksum("x0_iter1", x01.data(), x01.size())));
        meta.push_back("fnv1a_b_iter2=" + hex(checksum("b_iter2", b2.data(), b2.size())));
        meta.push_back("fnv1a_x0_iter2=" + hex(checksum("x0_iter2", x02.data(), x02.size())));
        meta.push_back("fnv1a_schur_diag_iter2=" + hex(checksum("schur_diag_iter2", schur2.data(), schur2.size())));
        meta.push_back("format=cf_dx_csr_chunk_v1 cfdx_vec_v1 cfdx_field_v1 values_percent_a");
        meta.push_back("validation_overall=" + std::string(all_ok ? "PASS" : "FAIL"));
        meta.push_back("elapsed_s=" + std::to_string(elapsed()));

        std::FILE* mf = audit_open(out_dir + "/meta.txt");
        for (const auto& line : meta) std::fprintf(mf, "%s\n", line.c_str());
        std::fclose(mf);

        std::cerr << "AUDIT_DONE validation_overall=" << (all_ok ? "PASS" : "FAIL")
                  << " elapsed_s=" << elapsed() << "\n";
        return all_ok ? 0 : 2;
    } catch (const std::exception& e) {
        std::cerr << "AUDIT_FAIL: " << e.what() << "\n";
        return 1;
    }
}
